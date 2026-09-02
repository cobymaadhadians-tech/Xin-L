#!/usr/bin/env python3
"""Standardize the revised OCD 2025 and CUD 2023 EUR auxiliaries for LDSC."""

import argparse
import csv
import gzip
import math
from pathlib import Path


ALLELES = {"A", "C", "G", "T"}
AMBIGUOUS = ({"A", "T"}, {"C", "G"})


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def load_hm3_snps(path: str) -> set[str]:
    result = set()
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if fields and fields[0].upper() != "SNP":
                result.add(fields[0])
    return result


def load_position_map(prefix: str) -> tuple[dict[tuple[int, int], str], int]:
    position_map = {}
    duplicate_positions = 0
    for chromosome in range(1, 23):
        path = f"{prefix}.{chromosome}.l2.ldscore.gz"
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                key = (int(row["CHR"]), int(row["BP"]))
                if key in position_map and position_map[key] != row["SNP"]:
                    duplicate_positions += 1
                    position_map[key] = ""
                else:
                    position_map[key] = row["SNP"]
    return position_map, duplicate_positions


def valid_alleles(a1: str, a2: str) -> bool:
    return a1 in ALLELES and a2 in ALLELES and a1 != a2


def keep_best(best, snp, n, z, a1, a2, stats):
    candidate = (abs(z), n, z, a1, a2)
    if snp in best:
        stats["duplicate_snp"] += 1
        if candidate[0] <= best[snp][0]:
            return
    best[snp] = candidate


def process_ocd(args, stats, best):
    hm3 = load_hm3_snps(args.hm3)
    with gzip.open(args.input, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SNP", "A1", "A2", "INFO", "OR", "SE", "P", "Neff_half"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"OCD input is missing columns: {sorted(missing)}")
        for row in reader:
            stats["input_rows"] += 1
            snp = (row.get("SNP") or "").strip()
            if snp not in hm3:
                continue
            stats["hm3_rows"] += 1
            a1 = (row.get("A1") or "").upper()
            a2 = (row.get("A2") or "").upper()
            info = parse_float(row.get("INFO"))
            odds_ratio = parse_float(row.get("OR"))
            se = parse_float(row.get("SE"))
            pval = parse_float(row.get("P"))
            n = parse_float(row.get("Neff_half"))
            beta = math.log(odds_ratio) if math.isfinite(odds_ratio) and odds_ratio > 0 else math.nan
            if not math.isfinite(info) or info < 0.9:
                stats["info_below_0.9"] += 1
                continue
            if not valid_alleles(a1, a2):
                stats["invalid_allele"] += 1
                continue
            if {a1, a2} in AMBIGUOUS:
                stats["strand_ambiguous"] += 1
                continue
            if not all(math.isfinite(x) for x in (beta, se, pval, n)) or se <= 0 or n <= 0 or pval <= 0 or pval > 1:
                stats["invalid_effect"] += 1
                continue
            keep_best(best, snp, n, beta / se, a1, a2, stats)


def process_cud(args, stats, best):
    position_map, duplicate_positions = load_position_map(args.weights_prefix)
    stats["reference_duplicate_position"] = duplicate_positions
    with gzip.open(args.input, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        required = {"CHR", "BP", "A1", "A2", "BETA", "SE", "P.value"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CUD input is missing columns: {sorted(missing)}")
        for row in reader:
            stats["input_rows"] += 1
            try:
                key = (int(row["CHR"]), int(row["BP"]))
            except (TypeError, ValueError):
                stats["invalid_position"] += 1
                continue
            snp = position_map.get(key, "")
            if not snp:
                continue
            stats["hm3_rows"] += 1
            a1 = (row.get("A1") or "").upper()
            a2 = (row.get("A2") or "").upper()
            beta = parse_float(row.get("BETA"))
            se = parse_float(row.get("SE"))
            pval = parse_float(row.get("P.value"))
            if not valid_alleles(a1, a2):
                stats["invalid_allele"] += 1
                continue
            if {a1, a2} in AMBIGUOUS:
                stats["strand_ambiguous"] += 1
                continue
            if not all(math.isfinite(x) for x in (beta, se, pval)) or se <= 0 or pval <= 0 or pval > 1:
                stats["invalid_effect"] += 1
                continue
            keep_best(best, snp, args.cud_neff, beta / se, a1, a2, stats)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", choices=["OCD_2025", "CUD_2023_EUR"], required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--hm3", required=True)
    parser.add_argument("--weights-prefix", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--cud-neff", type=float, default=161053.42428938235)
    args = parser.parse_args()

    stats = {
        "input_rows": 0,
        "hm3_rows": 0,
        "written_rows": 0,
        "duplicate_snp": 0,
        "invalid_position": 0,
        "invalid_allele": 0,
        "strand_ambiguous": 0,
        "invalid_effect": 0,
        "info_below_0.9": 0,
    }
    best = {}
    if args.trait == "OCD_2025":
        process_ocd(args, stats, best)
    else:
        process_cud(args, stats, best)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        handle.write("SNP\tN\tZ\tA1\tA2\n")
        for snp in sorted(best):
            _, n, z, a1, a2 = best[snp]
            handle.write(f"{snp}\t{n:.12g}\t{z:.12g}\t{a1}\t{a2}\n")
            stats["written_rows"] += 1

    audit = Path(args.audit)
    audit.parent.mkdir(parents=True, exist_ok=True)
    with audit.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in stats.items():
            handle.write(f"{key}\t{value}\n")
        handle.write(f"input\t{args.input}\n")
        handle.write(f"output\t{args.output}\n")
        if args.trait == "CUD_2023_EUR":
            handle.write(f"fixed_neff\t{args.cud_neff}\n")
    print(f"P25_AUX_PREP_COMPLETE\t{args.trait}\t{stats['written_rows']}")


if __name__ == "__main__":
    main()
