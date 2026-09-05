#!/usr/bin/env python3
"""Render the revised Figure 4b inferential summary panel."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CUSTOM_TABLE = ROOT / "data/supplementary_tables/Supplementary_Table_29_SCZ_stratified_GSEM_PI_brain_cell_custom.tsv"
DIRECT_TABLE = ROOT / "results/d6b_direct_contrast_20260905/D6B_direct_annotation_contrast.tsv"
RAW_TABLE = ROOT / "results/d6b_direct_contrast_20260905/D6B_raw_annotation_covariance_contrast.tsv"


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def add_box(ax: plt.Axes, y: float, title: str, lines: list[str], color: str) -> None:
    ax.text(0.04, y, title, transform=ax.transAxes, fontsize=27, fontweight="bold", color=color, va="top")
    ax.text(0.04, y - 0.065, "\n".join(lines), transform=ax.transAxes, fontsize=22, color="#222222", va="top", linespacing=1.35)


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def make_panel() -> plt.Figure:
    data = pd.read_csv(CUSTOM_TABLE, sep="\t")
    sig = as_bool(data["Bonferroni_significant_984_warning_free"])
    sig_data = data.loc[sig].copy()
    domain_counts = sig_data["target_parameter"].str.extract(r"F_(COMP|INT|SUD)", expand=False).value_counts().reindex(["COMP", "INT", "SUD"], fill_value=0)
    direct_text = "Direct PGC–FinnGen contrast: no result passed 0.05/492"
    if DIRECT_TABLE.exists():
        direct = pd.read_csv(DIRECT_TABLE, sep="\t")
        if "Bonferroni_significant_492_warning_free" in direct:
            n_direct = int(as_bool(direct["Bonferroni_significant_492_warning_free"]).sum())
            finite = int(pd.to_numeric(direct["P_two_sided"], errors="coerce").notna().sum())
            direct_text = f"Direct PGC–FinnGen contrast: {n_direct} of 492 tests passed 0.05/492 ({finite} finite)"
    raw_text = "Targeted raw covariance pilot: no QC-eligible result passed 0.05/29."
    raw_flag_text = ""
    if RAW_TABLE.exists():
        raw = pd.read_csv(RAW_TABLE, sep="\t")
        sig = as_bool(raw["Bonferroni_significant_domain_29"])
        primary = as_bool(raw["Bonferroni_significant_domain_29_primary"])
        eligible = raw.loc[primary, "Annotation"].drop_duplicates().tolist()
        flagged = raw.loc[sig & ~primary, "Annotation"].drop_duplicates().tolist()
        if eligible:
            labels = " and ".join(["exCA1", "exDG"] if set(eligible) >= {"exCA1L2", "exDGL2"} else eligible)
            raw_text = f"Targeted raw covariance pilot: {labels} SUD passed 0.05/29 after QC."
            if flagged:
                raw_flag_text = "PI Genes was flagged by model QC."

    fig = plt.figure(figsize=(7.2, 16.5), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.text(0.03, 0.97, "(b)", fontsize=31, fontweight="bold", color="#111111", va="top")
    ax.text(0.13, 0.97, "Inferential summary under the 984-test family", fontsize=29, fontweight="bold", color="#111111", va="top")
    ax.text(0.13, 0.925, "Primary calls require P < 0.05/984 and an annotation-model warning-free record.", fontsize=21, color="#444444", va="top")

    add_box(ax, 0.77, "Primary family", [
        "164 analyzable annotations × 6 SCZ–domain covariance parameters = 984 tests",
        "The 168- and 155-annotation thresholds are protocol sensitivities.",
    ], "#285B82")
    add_box(ax, 0.49, "Custom PI/brain-cell results", [
        "17 custom rows pass 0.05/984; 9 are warning-free.",
        f"Warning-free PGC rows by domain: COMP {int(domain_counts['COMP'])}, INT {int(domain_counts['INT'])}, SUD {int(domain_counts['SUD'])}.",
        "No custom warning-free global-family call was observed for narrow FinnGen SCZ.",
    ], "#2A9994")
    add_box(ax, 0.22, "Interpretation", [
        direct_text,
        raw_text,
        raw_flag_text,
        "The annotation analysis identifies candidate genomic contexts within each representation.",
        "It does not establish cell-specific effects, causal mechanisms or mediation of the PHBC contrast.",
    ], "#7B8F4E")
    return fig


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    configure_style()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    make_panel().savefig(args.output, format="svg", bbox_inches="tight")
    plt.close("all")
    print(args.output)


if __name__ == "__main__":
    main()
