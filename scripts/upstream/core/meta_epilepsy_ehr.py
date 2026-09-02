#!/usr/bin/env python3
"""Fixed-effect inverse-variance meta-analysis for the epilepsy EHR sources.

Inputs are the already standardized HM3/no-MHC tables for FinnGen R13,
MVP EUR, and UKB X345.  Alleles are harmonized per SNP, reverse-complement
matches are accepted, and beta signs are changed when the effect allele is
reversed.  The primary output keeps SNPs present in at least two cohorts;
the sensitivity output keeps SNPs present in all three cohorts.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
from collections import Counter
from typing import Dict, Iterable, List, Optional, Tuple


COHORTS = (
    ("FINNGEN", 62405.792541),
    ("MVP", 24332.120012),
    ("UKB", 20089.415662409814),
)


def comp(allele: str) -> str:
    return allele.translate(str.maketrans("ACGT", "TGCA"))


def open_gzip(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")


def parse_number(raw: str) -> Optional[float]:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def load_source(path: str, cohort: str, counters: Counter) -> Dict[str, Tuple[str, str, str, str, float, float, float, float]]:
    data: Dict[str, Tuple[str, str, str, str, float, float, float, float]] = {}
    with open_gzip(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        expected = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]
        if header[:9] != expected:
            raise ValueError(f"{cohort}: unexpected header {header}")
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                counters[f"{cohort}_malformed"] += 1
                continue
            snp, chrom, bp, a1, a2 = fields[:5]
            beta, se, p, n = (parse_number(x) for x in fields[5:9])
            if not snp or beta is None or se is None or se <= 0 or p is None or not 0 < p <= 1 or n is None or n <= 0:
                counters[f"{cohort}_invalid_core"] += 1
                continue
            if len(a1) != 1 or len(a2) != 1 or a1 not in "ACGT" or a2 not in "ACGT" or a1 == a2:
                counters[f"{cohort}_invalid_allele"] += 1
                continue
            snp = snp.lower()
            if snp in data:
                counters[f"{cohort}_duplicate_snp"] += 1
                continue
            data[snp] = (chrom, bp, a1, a2, beta, se, p, n)
            counters[f"{cohort}_rows"] += 1
    counters[f"{cohort}_unique_snp"] = len(data)
    return data


def harmonize(anchor: Tuple[str, str, str, str, float, float, float, float], candidate: Tuple[str, str, str, str, float, float, float, float]) -> Optional[Tuple[Tuple[str, str, str, str, float, float, float, float], int, str]]:
    """Return candidate with anchor alleles and beta sign, or None."""
    chrom, bp, a1, a2, beta, se, p, n = candidate
    _, _, aa1, aa2, _, _, _, _ = anchor
    pair = (a1, a2)
    target = (aa1, aa2)
    if pair == target:
        return (chrom, bp, aa1, aa2, beta, se, p, n), 1, "same"
    if pair == (aa2, aa1):
        return (chrom, bp, aa1, aa2, -beta, se, p, n), -1, "reverse"
    cp = (comp(a1), comp(a2))
    if cp == target:
        return (chrom, bp, aa1, aa2, beta, se, p, n), 1, "complement"
    if cp == (aa2, aa1):
        return (chrom, bp, aa1, aa2, -beta, se, p, n), -1, "reverse_complement"
    return None


def fmt(value: float) -> str:
    return f"{value:.12g}"


def write_meta(
    path: str,
    sources: List[Dict[str, Tuple[str, str, str, str, float, float, float, float]]],
    counters: Counter,
    min_cohorts: int,
    label: str,
) -> None:
    union = set().union(*(source.keys() for source in sources))
    output_part = path + ".part"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    rows = 0
    with gzip.open(output_part, "wt", encoding="utf-8", newline="") as out:
        out.write("SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tN\tN_cohort\tN_finngen\tN_mvp\tN_ukb\tdirection\n")
        for snp in sorted(union):
            present = [(i, source[snp]) for i, source in enumerate(sources) if snp in source]
            if len(present) < min_cohorts:
                continue
            anchor = present[0][1]
            harmonized: List[Tuple[int, Tuple[str, str, str, str, float, float, float, float], int, str]] = []
            for cohort_index, candidate in present:
                result = harmonize(anchor, candidate)
                if result is None:
                    counters[f"{label}_allele_mismatch"] += 1
                    continue
                harmonized.append((cohort_index, result[0], result[1], result[2]))
                counters[f"{label}_allele_{result[2]}"] += 1
            if len(harmonized) < min_cohorts:
                counters[f"{label}_below_min_after_harmonization"] += 1
                continue
            weights = [1.0 / (item[1][5] ** 2) for item in harmonized]
            beta_meta = sum(weight * item[1][4] for weight, item in zip(weights, harmonized)) / sum(weights)
            se_meta = math.sqrt(1.0 / sum(weights))
            z_meta = beta_meta / se_meta
            p_meta = math.erfc(abs(z_meta) / math.sqrt(2.0))
            n_meta = sum(COHORTS[item[0]][1] for item in harmonized)
            direction_by_cohort = ["." for _ in COHORTS]
            for cohort_index, record, _, _ in harmonized:
                direction_by_cohort[cohort_index] = "+" if record[4] > 0 else "-" if record[4] < 0 else "0"
            direction = "".join(direction_by_cohort)
            chrom, bp, a1, a2 = anchor[:4]
            out.write("\t".join((snp, chrom, bp, a1, a2, fmt(beta_meta), fmt(se_meta), fmt(p_meta), fmt(n_meta), str(len(harmonized)), fmt(COHORTS[0][1] if any(item[0] == 0 for item in harmonized) else 0.0), fmt(COHORTS[1][1] if any(item[0] == 1 for item in harmonized) else 0.0), fmt(COHORTS[2][1] if any(item[0] == 2 for item in harmonized) else 0.0), direction)) + "\n")
            rows += 1
    os.replace(output_part, path)
    counters[f"{label}_rows_written"] = rows
    counters[f"{label}_min_cohorts"] = min_cohorts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--finngen", required=True)
    parser.add_argument("--mvp", required=True)
    parser.add_argument("--ukb", required=True)
    parser.add_argument("--primary-output", required=True)
    parser.add_argument("--three-output", required=True)
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()

    counters = Counter()
    sources = [
        load_source(args.finngen, "FINNGEN", counters),
        load_source(args.mvp, "MVP", counters),
        load_source(args.ukb, "UKB", counters),
    ]
    write_meta(args.primary_output, sources, counters, 2, "META_2OF3")
    write_meta(args.three_output, sources, counters, 3, "META_3OF3")
    os.makedirs(os.path.dirname(args.audit), exist_ok=True)
    with open(args.audit, "w", encoding="utf-8", newline="") as handle:
        handle.write("metric\tvalue\n")
        for key, value in sorted(counters.items()):
            handle.write(f"{key}\t{value}\n")
        handle.write(f"input_FINNGEN\t{args.finngen}\ninput_MVP\t{args.mvp}\ninput_UKB\t{args.ukb}\nprimary_output\t{args.primary_output}\nthree_output\t{args.three_output}\nN_eff_FINNGEN\t{COHORTS[0][1]}\nN_eff_MVP\t{COHORTS[1][1]}\nN_eff_UKB\t{COHORTS[2][1]}\nmethod\tfixed-effect inverse-variance; no genomic control; primary requires at least 2 of 3 cohorts; sensitivity requires 3 of 3 cohorts\n")
    print(f"EHR_META_COMPLETE\tprimary={counters['META_2OF3_rows_written']}\tthree={counters['META_3OF3_rows_written']}", flush=True)


if __name__ == "__main__":
    main()
