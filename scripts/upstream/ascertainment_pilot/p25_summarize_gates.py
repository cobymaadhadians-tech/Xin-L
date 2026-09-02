#!/usr/bin/env python3
"""Summarize P2.5 h2 and target-auxiliary rg eligibility results."""

import csv
import math
import os
import re
from pathlib import Path


ROOT = Path(os.environ.get("ANALYSIS_ROOT", Path.cwd()))
H2_SOURCES = {
    "AUD": ROOT / "results/p3/aux_h2/AUD.log",
    "OCD_2025": ROOT / "results/p25/aux_h2/OCD_2025.log",
    "CUD_2023_EUR": ROOT / "results/p25/aux_h2/CUD_2023_EUR.log",
}


def h2_from_log(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Total Observed scale h2:\s+([-.0-9eE]+)\s+\(([-.0-9eE]+)\)", text)
    if not match:
        raise ValueError(f"h2 result missing in {path}")
    h2, se = map(float, match.groups())
    return h2, se, h2 / se


def rg_from_log(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if "Summary of Genetic Correlation Results" in line:
            header = lines[index + 1].split()
            values = lines[index + 2].split()
            row = dict(zip(header, values))
            return float(row["rg"]), float(row["se"]), float(row["z"]), float(row["p"])
    raise ValueError(f"rg result missing in {path}")


def main():
    out = ROOT / "results/p25/gates"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "updated_auxiliary_h2.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["trait", "h2_observed", "se", "z_h2", "eligible_z_gt_6"])
        for trait, path in H2_SOURCES.items():
            h2, se, z = h2_from_log(path)
            writer.writerow([trait, h2, se, z, z > 6])

    rows = []
    for path in sorted((ROOT / "results/p25/target_aux_rg").glob("*.log")):
        target, aux = path.stem.split("__", 1)
        rg, se, z, p = rg_from_log(path)
        rows.append((target, aux, rg, se, z, p, rg * rg < 0.5))
    with (out / "target_aux_rg_grid.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target", "auxiliary", "rg", "se", "z", "p", "eligible_rg2_lt_0.5"])
        writer.writerows(rows)

    families = {
        "MDD": ["MDD_CLIN_PGC2025", "MDD_EHR_PGC2025", "MDD_QUEST_PGC2025", "MDD_FINNGEN_R13"],
        "BD": ["BD_CLIN_PGC4", "BD_FINNGEN_R13"],
        "SCZ": ["SCZ_PGC2022", "SCZ_FINNGEN_R13"],
    }
    lookup = {(row[0], row[1]): row for row in rows}
    with (out / "family_common_eligibility.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["family", "auxiliary", "definitions_tested", "max_abs_rg", "all_rg2_lt_0.5"])
        for family, targets in families.items():
            auxiliaries = sorted({row[1] for row in rows if row[0] in targets})
            for aux in auxiliaries:
                selected = [lookup[(target, aux)] for target in targets]
                max_abs_rg = max(abs(row[2]) for row in selected)
                writer.writerow([family, aux, len(selected), max_abs_rg, all(row[6] for row in selected)])
    print("P25_GATE_SUMMARY_COMPLETE")


if __name__ == "__main__":
    main()
