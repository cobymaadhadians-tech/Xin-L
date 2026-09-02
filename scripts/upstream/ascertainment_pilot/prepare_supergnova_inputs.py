#!/usr/bin/env python3
"""Prepare coordinate-keyed SUPERGNOVA inputs and auditable FinnGen liftover BEDs."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path

from pyliftover import LiftOver


COMP = str.maketrans("ACGT", "TGCA")


def opener(path: Path, mode: str = "rt"):
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def finite(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


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


def load_reference(bim_pattern: str) -> dict[str, tuple[str, str]]:
    reference: dict[str, tuple[str, str]] = {}
    for chrom in range(1, 23):
        path = Path(bim_pattern.replace("@", str(chrom)))
        with path.open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 6:
                    raise ValueError(f"Malformed BIM row in {path}: {line.rstrip()}")
                key = fields[1]
                if key in reference:
                    raise ValueError(f"Duplicate reference SNP: {key}")
                reference[key] = (fields[4].upper(), fields[5].upper())
    return reference


def add_candidate(
    candidates: dict[str, tuple[str, str, float]],
    duplicates: set[str],
    key: str,
    a1: str,
    a2: str,
    z: float,
) -> None:
    if key in duplicates:
        return
    if key in candidates:
        candidates.pop(key)
        duplicates.add(key)
        return
    candidates[key] = (a1.upper(), a2.upper(), z)


def make_finngen_beds(source: Path, out_dir: Path, audit_path: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    handles = {
        chrom: (out_dir / f"finngen_chr{chrom}.hg38.bed").open("w")
        for chrom in range(1, 23)
    }
    counts = Counter()
    try:
        with opener(source) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"#chrom", "pos", "ref", "alt", "beta", "sebeta"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"FinnGen header lacks {sorted(required)}")
            for row in reader:
                counts["rows_read"] += 1
                try:
                    chrom = int(row["#chrom"])
                    pos = int(row["pos"])
                except ValueError:
                    counts["invalid_coordinate"] += 1
                    continue
                if chrom not in handles:
                    counts["non_autosomal"] += 1
                    continue
                a1, a2 = row["alt"].upper(), row["ref"].upper()
                if len(a1) != 1 or len(a2) != 1 or a1 not in "ACGT" or a2 not in "ACGT":
                    counts["non_snp"] += 1
                    continue
                beta, se = finite(row["beta"]), finite(row["sebeta"])
                if beta is None or se is None or se <= 0:
                    counts["invalid_z"] += 1
                    continue
                z = beta / se
                name = f"{a1}|{a2}|{z:.12g}"
                handles[chrom].write(f"chr{chrom}\t{pos - 1}\t{pos}\t{name}\t0\t+\n")
                counts["bed_rows_written"] += 1
    finally:
        for handle in handles.values():
            handle.close()
    write_audit(audit_path, counts, {"source": str(source), "coordinate_build": "GRCh38"})


def liftover_bed(source: Path, chain: Path, output: Path, audit_path: Path) -> None:
    converter = LiftOver(str(chain))
    counts = Counter()
    output.parent.mkdir(parents=True, exist_ok=True)
    with source.open() as source_handle, output.open("w") as output_handle:
        for line in source_handle:
            counts["rows_read"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                counts["malformed_bed_row"] += 1
                continue
            chrom, start, name, score, strand = fields[0], int(fields[1]), fields[3], fields[4], fields[5]
            mappings = converter.convert_coordinate(chrom, start, strand)
            if not mappings:
                counts["unmapped"] += 1
                continue
            if len(mappings) != 1:
                counts["multiple_mappings_excluded"] += 1
                continue
            target_chrom, target_start, target_strand, _ = mappings[0]
            target_start = int(target_start)
            output_handle.write(
                f"{target_chrom}\t{target_start}\t{target_start + 1}\t{name}\t{score}\t{target_strand}\n"
            )
            counts["rows_written"] += 1
            if target_strand == "-":
                counts["reverse_strand"] += 1
    write_audit(audit_path, counts, {
        "source": str(source),
        "chain": str(chain),
        "output": str(output),
        "coordinate_conversion": "GRCh38_to_GRCh37_0_based_point",
    })


def iter_pgc(path: Path, counts: Counter):
    with opener(path) as handle:
        while True:
            line = handle.readline()
            if not line:
                raise ValueError("PGC source has no data header")
            if not line.startswith("##"):
                header = line.rstrip("\n").split("\t")
                break
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for row in reader:
            counts["rows_read"] += 1
            try:
                chrom, bp = int(row["CHROM"]), int(row["POS"])
            except ValueError:
                counts["invalid_coordinate"] += 1
                continue
            beta, se = finite(row["BETA"]), finite(row["SE"])
            if beta is None or se is None or se <= 0:
                counts["invalid_z"] += 1
                continue
            yield f"{chrom}:{bp}", row["A1"], row["A2"], beta / se


def iter_adhd(path: Path, counts: Counter):
    with opener(path) as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        for row in reader:
            counts["rows_read"] += 1
            try:
                chrom, bp = int(row["CHR"]), int(row["BP"])
            except ValueError:
                counts["invalid_coordinate"] += 1
                continue
            odds_ratio, se = finite(row["OR"]), finite(row["SE"])
            if odds_ratio is None or odds_ratio <= 0 or se is None or se <= 0:
                counts["invalid_z"] += 1
                continue
            yield f"{chrom}:{bp}", row["A1"], row["A2"], math.log(odds_ratio) / se


def iter_lifted_finngen(pattern: str, counts: Counter):
    for chrom in range(1, 23):
        path = Path(pattern.replace("@", str(chrom)))
        with path.open() as handle:
            for line in handle:
                counts["lifted_rows_read"] += 1
                fields = line.split()
                if len(fields) < 6:
                    counts["malformed_lifted_row"] += 1
                    continue
                target = fields[0].removeprefix("chr")
                if not target.isdigit() or not 1 <= int(target) <= 22:
                    counts["non_autosomal_lift"] += 1
                    continue
                start, end = int(fields[1]), int(fields[2])
                if end - start != 1:
                    counts["non_single_base_lift"] += 1
                    continue
                name = fields[3].split("|")
                if len(name) != 3:
                    counts["malformed_lifted_name"] += 1
                    continue
                z = finite(name[2])
                if z is None:
                    counts["invalid_z"] += 1
                    continue
                yield f"{int(target)}:{end}", name[0], name[1], z


def write_audit(path: Path, counts: Counter, metadata: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("metric\tvalue\n")
        for key, value in metadata.items():
            handle.write(f"{key}\t{value}\n")
        for key in sorted(counts):
            handle.write(f"{key}\t{counts[key]}\n")


def build_sumstats(
    trait: str,
    source: Path | None,
    lifted_pattern: str | None,
    bim_pattern: str,
    scalar_n: int,
    output: Path,
    audit_path: Path,
) -> None:
    reference = load_reference(bim_pattern)
    counts = Counter(reference_snps=len(reference))
    candidates: dict[str, tuple[str, str, float]] = {}
    duplicates: set[str] = set()
    if trait == "PGC_SCZ":
        iterator = iter_pgc(source, counts)
        build = "GRCh37"
    elif trait == "ADHD":
        iterator = iter_adhd(source, counts)
        build = "GRCh37"
    elif trait == "FINNGEN_SCZ":
        iterator = iter_lifted_finngen(lifted_pattern, counts)
        build = "GRCh38_lifted_to_GRCh37"
    else:
        raise ValueError(trait)
    for key, a1, a2, z in iterator:
        if key not in reference:
            counts["coordinate_not_in_reference"] += 1
            continue
        counts["coordinate_in_reference"] += 1
        mode = allele_mode(a1, a2, *reference[key])
        if mode is None:
            counts["allele_mismatch"] += 1
            continue
        counts[f"allele_{mode}"] += 1
        add_candidate(candidates, duplicates, key, a1, a2, z)
    counts["duplicate_coordinate_excluded"] = len(duplicates)
    counts["rows_written"] = len(candidates)
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt") as handle:
        handle.write("SNP\tA1\tA2\tZ\tN\n")
        for key in sorted(candidates, key=lambda x: tuple(map(int, x.split(":")))):
            a1, a2, z = candidates[key]
            handle.write(f"{key}\t{a1}\t{a2}\t{z:.12g}\t{scalar_n}\n")
    metadata = {
        "trait": trait,
        "source": str(source) if source else lifted_pattern or "",
        "coordinate_build": build,
        "effect_allele": "A1",
        "scalar_effective_N": str(scalar_n),
        "output": str(output),
    }
    write_audit(audit_path, counts, metadata)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    bed = sub.add_parser("make-finngen-beds")
    bed.add_argument("--source", type=Path, required=True)
    bed.add_argument("--out-dir", type=Path, required=True)
    bed.add_argument("--audit", type=Path, required=True)
    lift = sub.add_parser("liftover-bed")
    lift.add_argument("--source", type=Path, required=True)
    lift.add_argument("--chain", type=Path, required=True)
    lift.add_argument("--output", type=Path, required=True)
    lift.add_argument("--audit", type=Path, required=True)
    build = sub.add_parser("build-sumstats")
    build.add_argument("--trait", choices=["PGC_SCZ", "FINNGEN_SCZ", "ADHD"], required=True)
    build.add_argument("--source", type=Path)
    build.add_argument("--lifted-pattern")
    build.add_argument("--bim-pattern", required=True)
    build.add_argument("--scalar-n", type=int, required=True)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "make-finngen-beds":
        make_finngen_beds(args.source, args.out_dir, args.audit)
    elif args.command == "liftover-bed":
        liftover_bed(args.source, args.chain, args.output, args.audit)
    else:
        if args.trait == "FINNGEN_SCZ" and not args.lifted_pattern:
            parser.error("FINNGEN_SCZ requires --lifted-pattern")
        if args.trait != "FINNGEN_SCZ" and not args.source:
            parser.error(f"{args.trait} requires --source")
        build_sumstats(args.trait, args.source, args.lifted_pattern, args.bim_pattern,
                       args.scalar_n, args.output, args.audit)


if __name__ == "__main__":
    main()
