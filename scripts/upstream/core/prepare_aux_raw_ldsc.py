#!/usr/bin/env python3
"""Convert the full raw ASD, AN and OCD summary statistics to LDSC format."""

import argparse
import csv
import gzip
import itertools
import math
from pathlib import Path


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.endswith(".gz") else open(path, "rt", encoding="utf-8", newline="")


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", choices=["ASD", "AN", "OCD"], required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--hm3", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()

    hm3 = set()
    with open(args.hm3, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().split()
            if fields and fields[0].upper() != "SNP":
                hm3.add(fields[0])

    if args.trait == "AN":
        header_starts = "CHROM\t"
        required = {"ID", "REF", "ALT", "BETA", "SE", "PVAL", "IMPINFO"}
        fixed_n = 52041.9101727871
    else:
        header_starts = "CHR\t"
        required = {"SNP", "A1", "A2", "INFO", "OR", "SE", "P"}
        fixed_n = 44368.0747340942 if args.trait == "ASD" else None

    stats = {"input_rows": 0, "hm3_rows": 0, "info_below_0.9": 0, "invalid_allele": 0, "strand_ambiguous": 0, "invalid_effect": 0, "duplicate_snp": 0, "written_rows": 0}
    best = {}
    with open_text(args.input) as handle:
        header_line = None
        for line in handle:
            if line.startswith(header_starts):
                header_line = line
                break
        if header_line is None:
            raise ValueError(f"summary-statistics header not found for {args.trait}")
        reader = csv.DictReader(itertools.chain([header_line], handle), delimiter="\t")
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing columns for {args.trait}: {sorted(missing)}")
        for row in reader:
            stats["input_rows"] += 1
            if args.trait == "AN":
                snp = (row.get("ID") or "").strip()
                a1 = (row.get("ALT") or "").strip().upper()
                a2 = (row.get("REF") or "").strip().upper()
                beta = parse_float(row.get("BETA"))
                pval = parse_float(row.get("PVAL"))
                n = fixed_n
                info = parse_float(row.get("IMPINFO"))
            else:
                snp = (row.get("SNP") or "").strip()
                a1 = (row.get("A1") or "").strip().upper()
                a2 = (row.get("A2") or "").strip().upper()
                odds_ratio = parse_float(row.get("OR"))
                beta = math.log(odds_ratio) if math.isfinite(odds_ratio) and odds_ratio > 0 else math.nan
                pval = parse_float(row.get("P"))
                info = parse_float(row.get("INFO"))
                n = fixed_n if fixed_n is not None else parse_float(row.get("Neff_half"))
            se = parse_float(row.get("SE"))

            if snp not in hm3:
                continue
            stats["hm3_rows"] += 1
            if not math.isfinite(info) or info < 0.9:
                stats["info_below_0.9"] += 1
                continue
            if a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"} or a1 == a2:
                stats["invalid_allele"] += 1
                continue
            if {a1, a2} in ({"A", "T"}, {"C", "G"}):
                stats["strand_ambiguous"] += 1
                continue
            if not all(math.isfinite(x) for x in (beta, se, pval, n)) or se <= 0 or n <= 0 or pval <= 0 or pval > 1:
                stats["invalid_effect"] += 1
                continue
            z = beta / se
            if not math.isfinite(z):
                stats["invalid_effect"] += 1
                continue
            candidate = (abs(z), n, z, a1, a2)
            if snp in best:
                stats["duplicate_snp"] += 1
                if candidate[0] <= best[snp][0]:
                    continue
            best[snp] = candidate

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        handle.write("SNP\tN\tZ\tA1\tA2\n")
        for snp in sorted(best):
            _, n, z, a1, a2 = best[snp]
            handle.write(f"{snp}\t{n:.10g}\t{z:.10g}\t{a1}\t{a2}\n")
            stats["written_rows"] += 1

    audit = Path(args.audit)
    audit.parent.mkdir(parents=True, exist_ok=True)
    with audit.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in stats.items():
            handle.write(f"{key}\t{value}\n")
    print(f"AUX_RAW_LDSC_PREP_COMPLETE\t{args.trait}\t{stats['written_rows']}")


if __name__ == "__main__":
    main()
