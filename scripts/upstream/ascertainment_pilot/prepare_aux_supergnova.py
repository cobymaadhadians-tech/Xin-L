#!/usr/bin/env python3
"""Prepare auxiliary GWAS files for SUPERGNOVA from audited SNP/Z tables.

The SUPERGNOVA implementation uses one scalar N when --N1/--N2 are supplied.
This script therefore records the source N distribution and writes the audited
scalar N explicitly to every output row.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path


ALLELES = {"A", "C", "G", "T"}
COMP = str.maketrans("ACGT", "TGCA")
AMBIGUOUS = ({"A", "T"}, {"C", "G"})


def opener(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open("rt", encoding="utf-8", newline="")


def number(value: str) -> float | None:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def allele_mode(a1: str, a2: str, r1: str, r2: str) -> str | None:
    a1, a2 = a1.upper(), a2.upper()
    if (a1, a2) == (r1, r2):
        return "direct"
    if (a1, a2) == (r2, r1):
        return "swapped"
    c1, c2 = a1.translate(COMP), a2.translate(COMP)
    if (c1, c2) == (r1, r2):
        return "complement"
    if (c1, c2) == (r2, r1):
        return "complement_swapped"
    return None


def load_reference(pattern: str, rsid_map_path: Path | None) -> tuple[dict[str, tuple[str, int, str, str]], dict[str, tuple[int, str, str]], dict[str, str]]:
    by_id: dict[str, tuple[str, int, str, str]] = {}
    by_coord: dict[str, tuple[int, str, str]] = {}
    for chrom in range(1, 23):
        path = Path(pattern.replace("@", str(chrom)))
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 6:
                    raise ValueError(f"Malformed BIM row in {path}: {line.rstrip()}")
                c, snp, bp = int(fields[0]), fields[1], int(fields[3])
                r1, r2 = fields[4].upper(), fields[5].upper()
                if snp in by_id:
                    raise ValueError(f"Duplicate reference SNP ID: {snp}")
                key = f"{c}:{bp}"
                if key in by_coord:
                    raise ValueError(f"Duplicate reference coordinate: {key}")
                by_id[snp] = (key, bp, r1, r2)
                by_coord[key] = (c, r1, r2)
    by_rsid: dict[str, str] = {}
    if rsid_map_path is not None:
        with rsid_map_path.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 4:
                    continue
                try:
                    if ":" in fields[0]:
                        c_text, bp_text = fields[0].removeprefix("chr").split(":", 1)
                        c, bp = int(c_text), int(bp_text)
                    else:
                        c, bp = int(fields[0]), int(fields[3])
                except ValueError:
                    continue
                key = f"{c}:{bp}"
                if key in by_coord:
                    by_rsid.setdefault(fields[1], key)
    return by_id, by_coord, by_rsid


def source_coordinate(snp: str, by_id: dict[str, tuple[str, int, str, str]], by_coord: dict[str, tuple[int, str, str]], by_rsid: dict[str, str]) -> str | None:
    if snp in by_id:
        return by_id[snp][0]
    if snp in by_rsid:
        return by_rsid[snp]
    raw = snp.removeprefix("chr")
    parts = raw.split(":")
    if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
        key = f"{int(parts[0])}:{int(parts[1])}"
        return key if key in by_coord else None
    return None


def write_audit(path: Path, metadata: dict[str, object], counts: Counter) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in metadata.items():
            handle.write(f"{key}\t{value}\n")
        for key in sorted(counts):
            handle.write(f"{key}\t{counts[key]}\n")


def build(args: argparse.Namespace) -> None:
    by_id, by_coord, by_rsid = load_reference(args.bim_pattern, args.rsid_map)
    counts = Counter(reference_snps=len(by_id))
    candidates: dict[str, tuple[str, str, float]] = {}
    duplicated: set[str] = set()
    source_ns: list[float] = []
    with opener(args.source) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SNP", "N", "Z", "A1", "A2"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{args.source} header lacks {sorted(missing)}")
        for row in reader:
            counts["rows_read"] += 1
            snp = (row.get("SNP") or "").strip()
            src_n, z = number(row.get("N")), number(row.get("Z"))
            if src_n is not None and src_n > 0:
                source_ns.append(src_n)
            if z is None:
                counts["invalid_z"] += 1
                continue
            a1, a2 = (row.get("A1") or "").upper(), (row.get("A2") or "").upper()
            if a1 not in ALLELES or a2 not in ALLELES or a1 == a2:
                counts["invalid_allele"] += 1
                continue
            if {a1, a2} in AMBIGUOUS:
                counts["strand_ambiguous"] += 1
                continue
            key = source_coordinate(snp, by_id, by_coord, by_rsid)
            if key is None:
                counts["snp_not_in_reference"] += 1
                continue
            counts["coordinate_in_reference"] += 1
            _, r1, r2 = by_coord[key]
            mode = allele_mode(a1, a2, r1, r2)
            if mode is None:
                counts["allele_mismatch"] += 1
                continue
            counts[f"allele_{mode}"] += 1
            if key in duplicated:
                continue
            if key in candidates:
                candidates.pop(key)
                duplicated.add(key)
                continue
            candidates[key] = (a1, a2, z)

    counts["duplicate_coordinate_excluded"] = len(duplicated)
    counts["rows_written"] = len(candidates)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        handle.write("SNP\tA1\tA2\tZ\tN\n")
        for key in sorted(candidates, key=lambda x: tuple(map(int, x.split(":")))):
            a1, a2, z = candidates[key]
            handle.write(f"{key}\t{a1}\t{a2}\t{z:.12g}\t{args.scalar_n}\n")
    source_n_median = "NA"
    if source_ns:
        ordered = sorted(source_ns)
        mid = len(ordered) // 2
        source_n_median = ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2
    metadata = {
        "trait": args.trait,
        "source": args.source,
        "coordinate_build": "GRCh37 via 1000G EUR BIM SNP-ID/coordinate map",
        "effect_allele": "A1",
        "source_N_median": source_n_median,
        "scalar_effective_N": args.scalar_n,
        "scalar_N_rule": args.scalar_rule,
        "reference_bim_pattern": args.bim_pattern,
        "reference_rsid_map": args.rsid_map,
        "output": args.output,
    }
    write_audit(Path(args.audit), metadata, counts)
    print(f"AUX_PREP_COMPLETE\t{args.trait}\trows={counts['rows_written']}\tscalar_N={args.scalar_n}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--bim-pattern", required=True)
    parser.add_argument("--rsid-map", type=Path, required=True)
    parser.add_argument("--scalar-n", type=int, required=True)
    parser.add_argument("--scalar-rule", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    build(args)


if __name__ == "__main__":
    main()
