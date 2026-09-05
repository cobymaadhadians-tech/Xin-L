#!/usr/bin/env python3
"""Render standalone Figure 4a for the D6B PI/brain-cell analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TABLE = REPOSITORY_ROOT / "data/supplementary_tables/Supplementary_Table_29_SCZ_stratified_GSEM_PI_brain_cell_custom.tsv"
DEFAULT_OUTPUT = REPOSITORY_ROOT / "outputs/Figure_4a_SCZ_STRATIFIED_PI_BRAIN_CELL.tif"

TEXT = "#111111"
GRID = "#D9D9D9"
WARNING = "#C57A00"
NEGATIVE = "#3B1F6F"
POSITIVE = "#B43E3E"
NEUTRAL = "#F7F7F7"
GROUP_BRAIN = "#2A6FAD"
GROUP_PI = "#3B1F6F"
GROUP_INTERSECTION = "#3F8A62"

TARGET_ORDER = [
    "F_COMP~~SCZ_PGC",
    "F_INT~~SCZ_PGC",
    "F_SUD~~SCZ_PGC",
    "F_COMP~~SCZ_FG",
    "F_INT~~SCZ_FG",
    "F_SUD~~SCZ_FG",
]
TARGET_LABELS = ["COMP", "INT", "SUD", "COMP", "INT", "SUD"]
ANNOTATION_ORDER = [
    "ASC1",
    "ASC2",
    "Endothelial Cells",
    "exCA1",
    "exCA3",
    "exDG",
    "exPFC1",
    "exPFC2",
    "GABA1",
    "GABA2",
    "Microglia",
    "Neuronal Stem Cells",
    "Oligodendrocytes",
    "Oligodendrocyte Precursor Cells",
    "PI Genes",
    "PI x ASC1",
    "PI x ASC2",
    "PI x Endothelial Cells",
    "PI x exCA1",
    "PI x exCA3",
    "PI x exDG",
    "PI x exPFC1",
    "PI x exPFC2",
    "PI x GABA1",
    "PI x GABA2",
    "PI x Microglia",
    "PI x Neuronal Stem Cells",
    "PI x Oligodendrocytes",
    "PI x Oligodendrocyte Precursor Cells",
]
DISPLAY_LABELS = {
    "Endothelial Cells": "Endothelial",
    "Neuronal Stem Cells": "NSC",
    "Oligodendrocytes": "Oligodendrocytes",
    "Oligodendrocyte Precursor Cells": "OPC",
    "PI Genes": "PI genes",
    "PI x Endothelial Cells": "PI × Endothelial",
    "PI x Neuronal Stem Cells": "PI × NSC",
    "PI x Oligodendrocytes": "PI × Oligodendrocytes",
    "PI x Oligodendrocyte Precursor Cells": "PI × OPC",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 8,
            "axes.titleweight": "bold",
            "axes.linewidth": 0.8,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.0,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def read_table() -> pd.DataFrame:
    data = pd.read_csv(TABLE, sep="\t")
    required = {
        "published_label",
        "target_parameter",
        "Enrichment",
        "Enrichment_p_value",
        "Annotation_model_warning",
        "Annotation_model_error",
        "Bonferroni_significant_168",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if len(data) != 174 or data["published_label"].nunique() != 29:
        raise ValueError("Expected 174 rows for 29 custom annotations")
    if set(data["target_parameter"]) != set(TARGET_ORDER):
        raise ValueError("Unexpected target-parameter set")
    if data.duplicated(["published_label", "target_parameter"]).any():
        raise ValueError("Duplicate annotation-target row")
    if data["Annotation_model_error"].map(as_bool).any():
        raise ValueError("Annotation-level errors are present in the custom table")
    return data


def format_annotation(label: str) -> str:
    return DISPLAY_LABELS.get(label, label.replace("PI x ", "PI × "))


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.tick_params(colors=TEXT, pad=2)
    ax.xaxis.label.set_color(TEXT)
    ax.yaxis.label.set_color(TEXT)
    ax.set_axisbelow(True)


def plot_matrix(ax: plt.Axes, data: pd.DataFrame) -> None:
    cmap = LinearSegmentedColormap.from_list("signed_enrichment", [NEGATIVE, NEUTRAL, POSITIVE])
    norm = TwoSlopeNorm(vmin=-30, vcenter=0, vmax=30)
    lookup = data.set_index(["published_label", "target_parameter"])
    x_values = np.arange(len(TARGET_ORDER))
    y_values = np.arange(len(ANNOTATION_ORDER))

    ax.set_xlim(-0.55, 5.55)
    ax.set_ylim(28.7, -0.7)
    ax.set_xticks(x_values, TARGET_LABELS)
    ax.set_yticks(y_values, [format_annotation(a) for a in ANNOTATION_ORDER])
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, lw=0.6)
    for y in y_values:
        ax.axhline(y, color="#EFEFEF", lw=0.45, zorder=0)
    ax.axvline(2.5, color="#777777", lw=0.9)
    ax.axhline(13.5, color="#777777", lw=0.8)
    ax.axhline(14.5, color="#777777", lw=0.8)

    for y, annotation in enumerate(ANNOTATION_ORDER):
        for x, target in enumerate(TARGET_ORDER):
            row = lookup.loc[(annotation, target)]
            enrichment = pd.to_numeric(row["Enrichment"], errors="coerce")
            p_value = pd.to_numeric(row["Enrichment_p_value"], errors="coerce")
            warning = as_bool(row["Annotation_model_warning"])
            significant = as_bool(row["Bonferroni_significant_168"])
            if not np.isfinite(enrichment):
                ax.scatter(x, y, marker="x", s=24, color="#8E8E8E", linewidths=0.85, zorder=4)
                continue
            # Strictly positive P values are required before applying log10.
            if not np.isfinite(p_value) or p_value <= 0:
                size = 18
            else:
                neg_log_p = max(float(-np.log10(p_value)), 0.0)
                size = 18 + 27 * min(neg_log_p, 6) / 6
            marker = "s" if warning else "o"
            edge = TEXT if significant else (WARNING if warning else "#8C8C8C")
            edge_width = 1.05 if significant else (0.9 if warning else 0.35)
            ax.scatter(
                x,
                y,
                s=size,
                marker=marker,
                c=[float(enrichment)],
                cmap=cmap,
                norm=norm,
                edgecolors=edge,
                linewidths=edge_width,
                alpha=0.96,
                zorder=4,
            )

    ax.text(1.0, 1.10, "PGC SCZ", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.text(4.0, 1.10, "FinnGen SCZ", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.text(-0.02, 1.015, "Enrichment", transform=ax.transAxes, ha="left", va="bottom", fontsize=7, color="#444444")

    ax.text(-0.23, 6.5, "Brain-cell\nannotations", transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=7, fontweight="bold", color=GROUP_BRAIN, clip_on=False)
    ax.text(-0.23, 14.0, "PI\ngenes", transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=7, fontweight="bold", color=GROUP_PI, clip_on=False)
    ax.text(-0.23, 21.5, "PI × cell\nannotations", transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=7, fontweight="bold", color=GROUP_INTERSECTION, clip_on=False)

    clean_axis(ax)


def make_figure(data: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(7.0866, 5.5118))  # 180 mm × 140 mm
    ax = fig.add_axes([0.30, 0.23, 0.67, 0.61])
    plot_matrix(ax, data)
    fig.text(0.075, 0.94, "a", fontsize=10, fontweight="bold", color=TEXT, va="bottom")
    fig.text(
        0.11,
        0.94,
        "Stratified Genomic SEM enrichment across PI/brain-cell annotations",
        fontsize=9,
        fontweight="bold",
        color=TEXT,
        va="bottom",
    )

    size_handles = [
        Line2D([0], [0], marker="o", color="#8C8C8C", mfc="#D9D9D9", lw=0, ms=3.7, label="−log10(P) = 1"),
        Line2D([0], [0], marker="o", color="#8C8C8C", mfc="#D9D9D9", lw=0, ms=5.2, label="3"),
        Line2D([0], [0], marker="o", color="#8C8C8C", mfc="#D9D9D9", lw=0, ms=6.7, label="6"),
    ]
    sig_handles = [
        Line2D([0], [0], marker="o", color=TEXT, mfc="#D9D9D9", lw=0, ms=5.0, label="Bonferroni-significant"),
        Line2D([0], [0], marker="s", color=WARNING, mfc="#D9D9D9", lw=0, ms=5.0, label="annotation warning"),
        Line2D([0], [0], marker="x", color="#8E8E8E", lw=0, ms=5.0, label="unavailable"),
    ]
    fig.legend(handles=size_handles, loc="lower left", bbox_to_anchor=(0.30, 0.075), ncol=3, frameon=False, fontsize=6.2, title="Marker size: −log10(enrichment P)", title_fontsize=6.5, handletextpad=0.35, columnspacing=1.0)
    fig.legend(handles=sig_handles, loc="lower right", bbox_to_anchor=(0.985, 0.075), ncol=3, frameon=False, fontsize=6.2, handletextpad=0.35, columnspacing=1.0)

    cax = fig.add_axes([0.11, 0.115, 0.13, 0.018])
    sm = mpl.cm.ScalarMappable(norm=TwoSlopeNorm(vmin=-30, vcenter=0, vmax=30), cmap=LinearSegmentedColormap.from_list("signed_enrichment_legend", [NEGATIVE, NEUTRAL, POSITIVE]))
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_ticks([-30, 0, 30])
    cbar.set_ticklabels(["−30", "0", "+30"])
    cbar.ax.tick_params(labelsize=5.8, length=2)
    cbar.set_label("Signed enrichment (clipped at ±30)", fontsize=6.2, labelpad=2)
    return fig


def save_figure(fig: plt.Figure, output: Path, pdf_output: Path | None, svg_output: Path | None) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.stem}.rgba.tif")
    fig.savefig(temporary, format="tiff", dpi=600, pil_kwargs={"compression": "raw"})
    with Image.open(temporary) as rendered:
        rendered.convert("RGB").save(output, format="TIFF", compression="raw", dpi=(600, 600))
    temporary.unlink()
    if pdf_output is not None:
        pdf_output.parent.mkdir(parents=True, exist_ok=True)
        pdf_path = pdf_output or output.with_name(f".{output.stem}.pdf")
        fig.savefig(pdf_path, format="pdf")
    if svg_output is not None:
        svg_output.parent.mkdir(parents=True, exist_ok=True)
        svg_path = svg_output or output.with_name(f".{output.stem}.svg")
        fig.savefig(svg_path, format="svg")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--pdf-output", type=Path, default=None)
    parser.add_argument("--svg-output", type=Path, default=None)
    args = parser.parse_args()
    configure_style()
    save_figure(make_figure(read_table()), args.output, args.pdf_output, args.svg_output)
    print(args.output)
    if args.pdf_output is not None:
        print(args.pdf_output)
    if args.svg_output is not None:
        print(args.svg_output)


if __name__ == "__main__":
    main()
