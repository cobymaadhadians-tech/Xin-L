#!/usr/bin/env python3
"""Build the referee-facing LDSC intercept and overlap audit table.

The script reads every ``*.log`` file in the LDSC audit log directory and
parses only values printed by LDSC. It never infers an intercept from a
genetic correlation or from an h2/rg result table. Missing logs remain NA and
are labelled as unavailable in the audit columns.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "data" / "supplementary_tables"
OUTPUTS = [
    TABLE_DIR / "Supplementary_Table_22_LDSC_diagnostic_and_overlap_audit.tsv",
    ROOT.parent
    / "submission_package"
    / "04_supplementary_information"
    / "tables"
    / "Supplementary_Table_22_LDSC_diagnostic_and_overlap_audit.tsv",
]
DEFAULT_LOG_DIR = ROOT / "results" / "ldsc_audit" / "logs"

AUDIT_COLUMNS = [
    "comparison_type",
    "section",
    "disease",
    "left_trait",
    "right_trait",
    "rg",
    "rg_se",
    "left_h2",
    "left_h2_se",
    "right_h2",
    "right_h2_se",
    "left_source_url",
    "right_source_identifier",
    "right_source_url",
    "intercept_gencov",
    "intercept_gencov_se",
    "intercept_left_h2",
    "intercept_left_h2_se",
    "intercept_right_h2",
    "intercept_right_h2_se",
    "intercept_constrained",
    "overlap_status",
    "archive_status",
    "notes",
    "ldsc_version",
    "intercept_command_mode",
    "log_file",
    "run_status",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value: str | None) -> str:
    return "NA" if value in (None, "") else value


def source_url_map(rows: list[dict[str, str]]) -> dict[str, str]:
    return {row["target"]: fmt(row.get("source_url")) for row in rows}


def overlap_map(rows: list[dict[str, str]]) -> dict[str, str]:
    return {row["target"]: fmt(row.get("overlap_group")) for row in rows}


def h2_map(rows: list[dict[str, str]]) -> dict[str, tuple[str, str]]:
    return {row["analysis_key"]: (fmt(row.get("h2")), fmt(row.get("se"))) for row in rows}


def pair_overlap(left: str, right: str, groups: dict[str, str]) -> str:
    left_group = groups.get(left, "NA")
    right_group = groups.get(right, "NA")
    if left_group != "NA" and left_group == right_group:
        return f"shared_source_group:{left_group}; participant_overlap_unknown"
    return "participant_overlap_unknown_or_unreported"


def parse_ldsc_log(path: Path) -> dict[str, str] | None:
    """Parse a pairwise LDSC log into h2 and genetic-covariance intercepts."""
    text = path.read_text(encoding="utf-8", errors="replace")
    versions = re.findall(r"\* Version\s+([^\s*]+)", text)
    intercepts = re.findall(
        r"^Intercept:\s*([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)",
        text,
        flags=re.MULTILINE,
    )
    if len(intercepts) < 3:
        return None
    command_match = re.search(r"(?ms)^Call:\s*(.*?)(?:\n\nBeginning analysis|\nBeginning analysis)", text)
    command = command_match.group(1) if command_match else ""
    constrained = bool(re.search(r"--(?:no-intercept|intercept-h2)", command))
    left, right = path.stem.split("__", 1) if "__" in path.stem else ("NA", "NA")
    return {
        "left": left,
        "right": right,
        "gencov": intercepts[2][0],
        "gencov_se": intercepts[2][1],
        "left_h2": intercepts[0][0],
        "left_h2_se": intercepts[0][1],
        "right_h2": intercepts[1][0],
        "right_h2_se": intercepts[1][1],
        "version": versions[0] if versions else "NA",
        "constrained": "TRUE" if constrained else "FALSE",
        "command_mode": "constrained" if constrained else "freely_estimated_default",
        "log_file": path.name,
    }


def read_log_index(log_dir: Path) -> dict[frozenset[str], dict[str, str]]:
    index: dict[frozenset[str], dict[str, str]] = {}
    if not log_dir.exists():
        return index
    for path in sorted(log_dir.glob("*.log")):
        parsed = parse_ldsc_log(path)
        if parsed is not None:
            index[frozenset({parsed["left"], parsed["right"]})] = parsed
    return index


def base_row(
    *,
    comparison_type: str,
    section: str,
    disease: str,
    left: str,
    right: str,
    rg: str,
    rg_se: str,
    h2: dict[str, tuple[str, str]],
    urls: dict[str, str],
    source_ids: dict[str, str],
    groups: dict[str, str],
    logs: dict[frozenset[str], dict[str, str]],
) -> dict[str, str]:
    left_h2, left_h2_se = h2.get(left, ("NA", "NA"))
    right_h2, right_h2_se = h2.get(right, ("NA", "NA"))
    audit = {
        "comparison_type": comparison_type,
        "section": section,
        "disease": disease,
        "left_trait": left,
        "right_trait": right,
        "rg": fmt(rg),
        "rg_se": fmt(rg_se),
        "left_h2": left_h2,
        "left_h2_se": left_h2_se,
        "right_h2": right_h2,
        "right_h2_se": right_h2_se,
        "left_source_url": urls.get(left, "NA"),
        "right_source_identifier": source_ids.get(right, "NA"),
        "right_source_url": urls.get(right, "NA"),
        "intercept_gencov": "NA",
        "intercept_gencov_se": "NA",
        "intercept_left_h2": "NA",
        "intercept_left_h2_se": "NA",
        "intercept_right_h2": "NA",
        "intercept_right_h2_se": "NA",
        "intercept_constrained": "FALSE",
        "overlap_status": pair_overlap(left, right, groups),
        "archive_status": "not_archived_in_current_project",
        "notes": "intercept not recoverable from current archive; LDSC command mode recorded as freely estimated",
        "ldsc_version": "NA",
        "intercept_command_mode": "freely_estimated_default",
        "log_file": "NA",
        "run_status": "log_unavailable",
    }
    parsed = logs.get(frozenset({left, right}))
    if parsed is not None:
        audit.update(
            {
                "intercept_gencov": parsed["gencov"],
                "intercept_gencov_se": parsed["gencov_se"],
                "intercept_left_h2": parsed["left_h2"] if parsed["left"] == left else parsed["right_h2"],
                "intercept_left_h2_se": parsed["left_h2_se"] if parsed["left"] == left else parsed["right_h2_se"],
                "intercept_right_h2": parsed["right_h2"] if parsed["left"] == left else parsed["left_h2"],
                "intercept_right_h2_se": parsed["right_h2_se"] if parsed["left"] == left else parsed["left_h2_se"],
                "intercept_constrained": parsed["constrained"],
                "archive_status": "archived_ldsc_log",
                "notes": "parsed from archived LDSC pairwise log; intercepts freely estimated unless command indicates a constraint",
                "ldsc_version": parsed["version"],
                "intercept_command_mode": parsed["command_mode"],
                "log_file": parsed["log_file"],
                "run_status": "log_parsed",
            }
        )
    return audit


def build_rows(log_dir: Path) -> list[dict[str, str]]:
    provenance = read_tsv(TABLE_DIR / "Supplementary_Table_1_target_GWAS_provenance.tsv")
    target_h2 = read_tsv(TABLE_DIR / "Supplementary_Table_3_target_heritability_gate.tsv")
    auxiliary_h2 = read_tsv(TABLE_DIR / "Supplementary_Table_2_psychiatric_auxiliary_provenance_and_h2.tsv")
    target_aux = read_tsv(TABLE_DIR / "Supplementary_Table_4_psychiatric_target_auxiliary_eligibility.tsv")
    target_pairs = read_tsv(TABLE_DIR / "Supplementary_Table_6_cross_definition_genetic_correlations.tsv")

    urls = source_url_map(provenance)
    groups = overlap_map(provenance)
    h2 = h2_map(target_h2)
    source_ids = {row["trait"]: row["source_identifier"] for row in auxiliary_h2}
    for row in auxiliary_h2:
        h2[row["trait"]] = (fmt(row.get("h2")), fmt(row.get("se")))

    logs = read_log_index(log_dir)
    rows: list[dict[str, str]] = []
    for row in target_aux:
        audit = base_row(
            comparison_type="target_auxiliary",
            section="psychiatric_target_auxiliary",
            disease="NA",
            left=row["target"],
            right=row["auxiliary"],
            rg=row["rg"],
            rg_se=row["se"],
            h2=h2,
            urls=urls,
            source_ids=source_ids,
            groups=groups,
            logs=logs,
        )
        audit["overlap_status"] = "participant_overlap_unknown_or_unreported"
        audit["notes"] = "target-auxiliary overlap was not resolved from the available public source records; LDSC command mode recorded as freely estimated"
        rows.append(audit)

    for row in target_pairs:
        audit = base_row(
            comparison_type="target_target",
            section=row["section"],
            disease=row["disease"],
            left=row["left_target"],
            right=row["right_target"],
            rg=row["rg"],
            rg_se=row["rg_se"],
            h2=h2,
            urls=urls,
            source_ids=source_ids,
            groups=groups,
            logs=logs,
        )
        parsed = logs.get(frozenset({row["left_target"], row["right_target"]}))
        if parsed is not None:
            audit["overlap_status"] = (
                "documented_shared_FinnGen_cohorts_and_cases"
                if frozenset({row["left_target"], row["right_target"]})
                == frozenset({"SCZ_FINNGEN_R13_EXMORE", "SCZ_FINNGEN_R13"})
                else "participant_overlap_not_reported_in_log; source records kept separate"
            )
        rows.append(audit)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR)
    args = parser.parse_args()
    rows = build_rows(args.log_dir)
    for output in OUTPUTS:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=AUDIT_COLUMNS, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {len(OUTPUTS)} outputs")


if __name__ == "__main__":
    main()
