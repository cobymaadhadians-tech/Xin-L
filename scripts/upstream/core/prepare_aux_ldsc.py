#!/usr/bin/env python3
"""Convert project-normalized auxiliary traits to LDSC sumstats format."""

import argparse
import csv
import gzip
import math
from pathlib import Path


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.endswith(".gz") else open(path, "rt", encoding="utf-8", newline="")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--hm3", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()

    hm3 = set()
    with open(args.hm3, "rt", encoding="utf-8") as handle:
        for line in handle:
            token = line.strip().split()[0] if line.strip() else ""
            if token and token.upper() != "SNP":
                hm3.add(token)

    stats = {
        "input_rows": 0,
        "hm3_rows": 0,
        "written_rows": 0,
        "duplicate_snp": 0,
        "invalid_snp": 0,
        "invalid_allele": 0,
        "strand_ambiguous": 0,
        "invalid_effect": 0,
        "invalid_n": 0,
    }
    best = {}
    with open_text(args.input) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SNP", "A1", "A2", "BETA", "SE", "P", "N_RSS"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing columns: {sorted(missing)}")
        for row in reader:
            stats["input_rows"] += 1
            snp = (row.get("SNP") or "").strip().strip('"')
            if snp not in hm3:
                continue
            stats["hm3_rows"] += 1
            a1 = (row.get("A1") or "").strip().strip('"').upper()
            a2 = (row.get("A2") or "").strip().strip('"').upper()
            if len(snp) < 2 or not snp.startswith("rs"):
                stats["invalid_snp"] += 1
                continue
            if a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"} or a1 == a2:
                stats["invalid_allele"] += 1
                continue
            if {a1, a2} in ({"A", "T"}, {"C", "G"}):
                stats["strand_ambiguous"] += 1
                continue
            try:
                beta = float(row["BETA"])
                se = float(row["SE"])
                pval = float(row["P"])
                n = float(row["N_RSS"])
            except (TypeError, ValueError):
                stats["invalid_effect"] += 1
                continue
            if not all(math.isfinite(x) for x in (beta, se, pval, n)) or se <= 0 or n <= 0 or pval <= 0 or pval > 1:
                stats["invalid_effect"] += 1
                continue
            z = beta / se
            if not math.isfinite(z):
                stats["invalid_effect"] += 1
                continue
            candidate = (abs(z), snp, n, z, a1, a2)
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
            _, _, n, z, a1, a2 = best[snp]
            handle.write(f"{snp}\t{n:.10g}\t{z:.10g}\t{a1}\t{a2}\n")
            stats["written_rows"] += 1

    audit = Path(args.audit)
    audit.parent.mkdir(parents=True, exist_ok=True)
    with audit.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in stats.items():
            handle.write(f"{key}\t{value}\n")
    print("\t".join(["AUX_LDSC_PREP_COMPLETE", str(output), str(stats["written_rows"])]))


if __name__ == "__main__":
    main()
