#!/usr/bin/env python3
"""Check the files and columns required by the figure renderers."""

from __future__ import annotations

import csv
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

REQUIRED_COLUMNS = {
    "data/derived/phbc/integrated_cross_definition_rg.tsv": {
        "comparison",
        "rg",
        "rg_se",
    },
    "data/derived/phbc/integrated_target_phbc.tsv": {
        "section",
        "target",
        "phbc",
        "official_se",
    },
    "data/derived/phbc/integrated_paired_delta_phbc.tsv": {
        "comparison",
        "delta_pp",
        "ci95_low_pp",
        "ci95_high_pp",
    },
    "data/derived/phbc/leave_one_auxiliary_out_paired_phbc.tsv": {
        "comparison",
        "condition",
        "delta_pp",
    },
    "data/derived/phbc/leave_category_out_paired_did.tsv": {
        "disease",
        "domain",
        "left_reduction_pp",
        "right_reduction_pp",
        "difference_in_reduction_pp",
    },
    "data/derived/phbc/broad_scz_figure1.tsv": {
        "label",
        "disease",
        "rg",
        "delta_pp",
    },
    "data/derived/supergnova/pgc_finngen_scz_local_covariance.tsv": {
        "chr",
        "start",
        "end",
        "rho",
        "bonferroni",
    },
    "data/derived/supergnova/scz_adhd_local_covariance.tsv": {
        "chr",
        "start",
        "end",
        "rho_PGC_ADHD",
        "rho_FG_ADHD",
    },
    "data/derived/supergnova/scz_adhd_profile_comparison.tsv": {
        "pearson_r",
        "informative_sign_concordance",
        "top5pct_n_each",
        "top5pct_intersection_n",
    },
    "data/derived/supergnova/six_trait/sixtrait_profile_summary.tsv": {
        "auxiliary",
        "pearson_r",
    },
    "data/derived/supergnova/six_trait/panel_standardized_profile_similarity.tsv": {
        "panel_standardized_profile_pearson_r",
    },
}


def columns(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return set(next(csv.reader(handle, delimiter="\t")))


def main() -> int:
    missing: list[str] = []
    invalid: list[str] = []
    for relative_path, expected in REQUIRED_COLUMNS.items():
        path = REPOSITORY_ROOT / relative_path
        if not path.is_file():
            missing.append(relative_path)
            continue
        actual = columns(path)
        absent = sorted(expected - actual)
        if absent:
            invalid.append(f"{relative_path}: missing columns {', '.join(absent)}")

    if missing or invalid:
        for path in missing:
            print(f"MISSING\t{path}")
        for message in invalid:
            print(f"INVALID\t{message}")
        return 1

    print(f"OK\t{len(REQUIRED_COLUMNS)} required TSV files and column sets found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
