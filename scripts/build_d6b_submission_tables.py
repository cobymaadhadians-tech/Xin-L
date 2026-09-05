#!/usr/bin/env python3
"""Materialize referee-facing D6B supplementary tables from the archived run."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "results/d6b_published_pi_braincell_stratified_gsem_20260905"
TABLE_DIR = ROOT / "data/supplementary_tables"
OUTPUT = TABLE_DIR / "Supplementary_Table_30_SCZ_stratified_GSEM_PI_brain_cell_provenance_QC.tsv"
FULL_OUTPUT = TABLE_DIR / "Supplementary_Table_28_SCZ_stratified_GSEM_PI_brain_cell_full.tsv"
CUSTOM_OUTPUT = TABLE_DIR / "Supplementary_Table_29_SCZ_stratified_GSEM_PI_brain_cell_custom.tsv"


def _public_annotation_value(value: object) -> object:
    """Remove private reference-root prefixes while retaining annotation IDs."""
    if isinstance(value, str) and value.startswith("/public/home/"):
        return value.rstrip("/").rsplit("/", 1)[-1]
    return value


def _materialize_public_enrichment(source: Path, output: Path) -> None:
    frame = pd.read_csv(source, sep="\t")
    for column in frame.columns:
        if pd.api.types.is_object_dtype(frame[column]) or pd.api.types.is_string_dtype(frame[column]):
            frame[column] = frame[column].map(_public_annotation_value)
    frame.to_csv(output, sep="\t", index=False)


def main() -> None:
    long = pd.read_csv(ARCHIVE / "D6B_enrichment_long.tsv", sep="\t")
    manifest = pd.read_csv(ARCHIVE / "D6B_annotation_manifest.tsv", sep="\t")
    tested = long.loc[long["is_nonbase"].eq(True)].copy()

    _materialize_public_enrichment(ARCHIVE / "D6B_enrichment_long.tsv", FULL_OUTPUT)
    _materialize_public_enrichment(ARCHIVE / "D6B_custom_enrichment.tsv", CUSTOM_OUTPUT)

    rows: list[dict[str, object]] = []

    def add(**kwargs: object) -> None:
        rows.append(
            {
                "record_type": kwargs.get("record_type", ""),
                "record_id": kwargs.get("record_id", ""),
                "annotation_class": kwargs.get("annotation_class", ""),
                "published_label": kwargs.get("published_label", ""),
                "source_label": kwargs.get("source_label", ""),
                "source_protocol": kwargs.get("source_protocol", ""),
                "n_input_annotations": kwargs.get("n_input_annotations", ""),
                "n_s_ldsc_matrices": kwargs.get("n_s_ldsc_matrices", ""),
                "n_tested_rows": kwargs.get("n_tested_rows", ""),
                "n_finite_p": kwargs.get("n_finite_p", ""),
                "n_significant_168": kwargs.get("n_significant_168", ""),
                "n_warning_rows": kwargs.get("n_warning_rows", ""),
                "n_error_rows": kwargs.get("n_error_rows", ""),
                "notes": kwargs.get("notes", ""),
            }
        )

    add(
        record_type="run_summary",
        record_id="D6B",
        source_protocol="Grotzinger et al. 2022 PI/brain-cell annotation protocol",
        n_input_annotations=210,
        n_s_ldsc_matrices=197,
        n_tested_rows=len(tested),
        n_finite_p=int(tested["Enrichment_p_value"].notna().sum()),
        n_significant_168=int(tested["Bonferroni_significant_168"].eq(True).sum()),
        n_warning_rows=int(tested["Annotation_model_warning"].eq(True).sum()),
        n_error_rows=int(tested["Annotation_model_error"].eq(True).sum()),
        notes="200 LDSC blocks; published 0.05/168 and covariance-protocol 0.05/155 thresholds gave identical calls; baseline continuous/flanking annotations were retained for conditioning and excluded from enrichment output.",
    )

    class_inputs = {"baseline": 97, "expression": 24, "chromatin": 60, "custom": 29}
    for annotation_class, group in tested.groupby("annotation_class", sort=False):
        add(
            record_type="class_summary",
            record_id=annotation_class,
            annotation_class=annotation_class,
            n_input_annotations=class_inputs.get(annotation_class, ""),
            n_tested_rows=len(group),
            n_finite_p=int(group["Enrichment_p_value"].notna().sum()),
            n_significant_168=int(group["Bonferroni_significant_168"].eq(True).sum()),
            n_warning_rows=int(group["Annotation_model_warning"].eq(True).sum()),
            n_error_rows=int(group["Annotation_model_error"].eq(True).sum()),
            notes="Annotation family summary after official enrich() exclusions.",
        )

    for target, group in tested.groupby("target_parameter", sort=False):
        add(
            record_type="target_summary",
            record_id=target,
            n_tested_rows=len(group),
            n_finite_p=int(group["Enrichment_p_value"].notna().sum()),
            n_significant_168=int(group["Bonferroni_significant_168"].eq(True).sum()),
            n_warning_rows=int(group["Annotation_model_warning"].eq(True).sum()),
            n_error_rows=int(group["Annotation_model_error"].eq(True).sum()),
            source_protocol="PGC SCZ" if "SCZ_PGC" in target else "FinnGen SCZ",
            notes="Six annotation-dependent COMP, INT and SUD covariance parameters were evaluated.",
        )

    custom = manifest.loc[manifest["annotation_class"].eq("custom")].copy()
    custom = custom.drop_duplicates(subset=["published_label"], keep="first")
    for _, row in custom.iterrows():
        add(
            record_type="annotation_provenance",
            record_id=str(row["published_label"]),
            annotation_class="custom",
            published_label=str(row["published_label"]),
            source_label=str(row["source_label"]),
            source_protocol="Grotzinger et al. 2022, Supplementary Table 10",
            n_input_annotations=1,
            notes="Published PI, brain-cell or PI × brain-cell annotation; custom annotation reconstructed with the reference protocol.",
        )

    result = pd.DataFrame(rows)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, sep="\t", index=False)
    print(OUTPUT)


if __name__ == "__main__":
    main()
