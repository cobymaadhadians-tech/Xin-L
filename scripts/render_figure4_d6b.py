#!/usr/bin/env python3
"""Render the main Figure 4 as a Fig. 3-style custom-annotation profile."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image
import seaborn as sns


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TABLE = REPOSITORY_ROOT / "data/supplementary_tables/Supplementary_Table_29_SCZ_stratified_GSEM_PI_brain_cell_custom.tsv"
FIGURE_STEM = "Figure_4_SCZ_STRATIFIED_PI_BRAIN_CELL"
MM = 1 / 25.4
FIG_WIDTH_MM = 180
FIG_HEIGHT_MM = 125

TEXT = "#111111"
PGC = "#5B7891"
TEAL = "#2A9994"
CORAL = "#C15B5A"

DOMAIN_ORDER = ["COMP", "INT", "SUD"]
ANNOTATION_ORDER = [
    "PI Genes",
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
    "PI Genes": "PI genes",
    "Endothelial Cells": "Endothelial",
    "Neuronal Stem Cells": "NSC",
    "Oligodendrocytes": "Oligodendrocytes",
    "Oligodendrocyte Precursor Cells": "OPC",
    "PI x Endothelial Cells": "PI × Endothelial",
    "PI x Neuronal Stem Cells": "PI × NSC",
    "PI x Oligodendrocytes": "PI × Oligodendrocytes",
    "PI x Oligodendrocyte Precursor Cells": "PI × OPC",
}


def configure_style() -> None:
    sns.set_theme(style="white", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.labelsize": 6.8,
            "axes.titlesize": 7.2,
            "axes.linewidth": 1.0,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 6.2,
            "ytick.labelsize": 6.2,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.tick_params(colors=TEXT, pad=2)
    ax.set_axisbelow(True)


def format_annotation(label: str) -> str:
    return DISPLAY_LABELS.get(label, label.replace("PI x ", "PI × "))


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def read_table() -> pd.DataFrame:
    data = pd.read_csv(TABLE, sep="\t")
    required = {
        "published_label",
        "target_parameter",
        "Enrichment",
        "Enrichment_SE",
        "Enrichment_p_value",
        "Annotation_model_warning",
        "Annotation_model_error",
        "Bonferroni_significant_984_warning_free",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if len(data) != 174:
        raise ValueError(f"Expected 174 custom rows, found {len(data)}")
    if data["published_label"].nunique() != 29:
        raise ValueError("Expected 29 unique custom annotations")
    if data.duplicated(["published_label", "target_parameter"]).any():
        raise ValueError("Duplicate annotation-target row")
    for column in [
        "Annotation_model_warning",
        "Annotation_model_error",
        "Bonferroni_significant_984_warning_free",
    ]:
        data[column] = data[column].map(as_bool)
    return data


def significant_annotation_order(data: pd.DataFrame) -> list[str]:
    significant = set(data.loc[data["Bonferroni_significant_984_warning_free"], "published_label"])
    selected = [annotation for annotation in ANNOTATION_ORDER if annotation in significant]
    if len(selected) != 9:
        raise ValueError(f"Expected 9 warning-free global-family significant annotations, found {len(selected)}")
    return selected


def target_parameter(domain: str, representation: str) -> str:
    return f"F_{domain}~~SCZ_{representation}"


def draw_estimate(
    ax: plt.Axes,
    row: pd.Series,
    y: float,
    color: str,
    significant: bool,
    warning: bool,
) -> None:
    estimate = pd.to_numeric(row["Enrichment"], errors="coerce")
    se = pd.to_numeric(row["Enrichment_SE"], errors="coerce")
    if not np.isfinite(estimate):
        ax.text(0, y, "×", ha="center", va="center", fontsize=7, color="#8A8A8A", zorder=4)
        return

    point_clipped = bool(estimate < -30 or estimate > 30)
    display_estimate = float(np.clip(estimate, -29.4, 29.4)) if point_clipped else float(estimate)
    marker = "|" if point_clipped else "o"
    if np.isfinite(se) and se >= 0:
        lo = float(np.clip(estimate - se, -30, 30))
        hi = float(np.clip(estimate + se, -30, 30))
        xerr = np.array([[max(display_estimate - lo, 0)], [max(hi - display_estimate, 0)]])
    else:
        xerr = None

    clipped = point_clipped
    if np.isfinite(se) and se >= 0:
        clipped = clipped or bool(estimate - se < -30 or estimate + se > 30)
    alpha = 1.0 if significant else 0.52
    marker_size = 4.8 if significant else 3.2
    line_width = 0.9 if significant else 0.55
    cap_size = 1.5 if significant else 1.0
    if clipped:
        alpha *= 0.65
        marker_size = min(marker_size, 3.0)
        line_width = 0.50
        cap_size = 1.0
    ax.errorbar(
        display_estimate,
        y,
        xerr=xerr,
        fmt=marker,
        ms=marker_size,
        mfc=color if significant else "white",
        mec=color,
        mew=0.80 if significant else 0.55,
        ecolor=color,
        elinewidth=line_width,
        capsize=cap_size,
        alpha=alpha,
        zorder=3,
    )
    if warning:
        warning_x = float(np.clip(display_estimate + (0.85 if display_estimate < 28 else -1.0), -29.2, 29.2))
        ax.text(warning_x, y + 0.06, "†", ha="center", va="bottom", fontsize=5.2, color=CORAL, alpha=0.68, zorder=5)


def draw_profile(ax: plt.Axes, data: pd.DataFrame, domain: str, annotations: list[str], show_y: bool) -> None:
    lookup = data.set_index(["published_label", "target_parameter"])
    positions = np.arange(len(annotations))[::-1]
    offsets = {"PGC": 0.105, "FG": -0.105}
    colors = {"PGC": PGC, "FG": TEAL}

    ax.axvline(0, color="#A6A6A6", lw=0.9, ls=(0, (3, 3)), zorder=0)

    base_count = sum(not annotation.startswith("PI x ") for annotation in annotations)
    if 0 < base_count < len(annotations):
        boundary = (positions[base_count - 1] + positions[base_count]) / 2
        ax.axhline(boundary, color="#AFAFAF", lw=0.8, zorder=1)

    for annotation, y in zip(annotations, positions):
        for representation in ["PGC", "FG"]:
            row = lookup.loc[(annotation, target_parameter(domain, representation))]
            draw_estimate(
                ax=ax,
                row=row,
                y=y + offsets[representation],
                color=colors[representation],
                significant=as_bool(row["Bonferroni_significant_984_warning_free"]),
                warning=as_bool(row["Annotation_model_warning"]),
            )

    ax.set_title(domain, fontsize=7.0, fontweight="bold", pad=3)
    ax.set_xlim(-30, 30)
    ax.set_xticks([-30, -15, 0, 15, 30])
    ax.set_xlabel("Signed enrichment")
    ax.set_ylim(-0.55, len(annotations) - 0.45)
    ax.set_yticks(positions, [format_annotation(annotation) for annotation in annotations] if show_y else [])
    ax.tick_params(axis="x", length=2.5)
    ax.tick_params(axis="y", length=0, labelsize=5.8)
    if not show_y:
        ax.tick_params(axis="y", labelleft=False)
    clean_axis(ax)


def make_figure(data: pd.DataFrame) -> plt.Figure:
    annotations = significant_annotation_order(data)
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, FIG_HEIGHT_MM * MM), facecolor="white")
    gs = fig.add_gridspec(1, 3, left=0.275, right=0.985, top=0.88, bottom=0.20, wspace=0.30)
    axes = []
    for i in range(3):
        axes.append(fig.add_subplot(gs[0, i]))
    for i, (ax, domain) in enumerate(zip(axes, DOMAIN_ORDER)):
        draw_profile(ax, data, domain, annotations, show_y=i == 0)

    axes[0].text(-0.39, 1.08, "(a)", transform=axes[0].transAxes, ha="left", va="bottom", fontsize=8, fontweight="bold")
    axes[0].text(-0.26, 1.08, "Warning-free global-family profile", transform=axes[0].transAxes, ha="left", va="bottom", fontsize=7, fontweight="bold")

    legend_handles = [
        Line2D([0], [0], marker="o", color=PGC, markerfacecolor=PGC, markeredgecolor=PGC, lw=0, ms=4.8, label="PGC SCZ"),
        Line2D([0], [0], marker="o", color=TEAL, markerfacecolor=TEAL, markeredgecolor=TEAL, lw=0, ms=4.8, label="narrow FinnGen SCZ"),
        Line2D([0], [0], marker="o", color="#555555", markerfacecolor="white", markeredgecolor="#555555", lw=0, ms=4.0, label="open: non-significant"),
        Line2D([0], [0], marker="$†$", color=CORAL, lw=0, ms=6.5, label="annotation-model warning"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.57, 0.045),
        ncol=4,
        frameon=False,
        fontsize=5.7,
        handlelength=1.0,
        handletextpad=0.35,
        columnspacing=1.0,
    )
    return fig


def save_figure(fig: plt.Figure, outdir: Path, pdf: bool, emit_vector: bool) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    if pdf:
        output = outdir / f"{FIGURE_STEM}.pdf"
        fig.savefig(output, format="pdf", bbox_inches="tight")
    else:
        temporary = outdir / f".{FIGURE_STEM}.rgba.tif"
        output = outdir / f"{FIGURE_STEM}.tif"
        fig.savefig(temporary, format="tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "raw"})
        with Image.open(temporary) as rendered:
            rendered.convert("RGB").save(output, format="TIFF", compression="raw", dpi=(600, 600))
        temporary.unlink()
    if emit_vector:
        fig.savefig(outdir / f".{FIGURE_STEM}.svg", format="svg", bbox_inches="tight")
        if not pdf:
            fig.savefig(outdir / f".{FIGURE_STEM}.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("submission_package/03_figures/main"))
    parser.add_argument("--pdf", action="store_true", help="Write the user-facing figure as a vector PDF")
    parser.add_argument("--emit-vector", action="store_true", help="Write hidden SVG/PDF copies for QA")
    args = parser.parse_args()
    configure_style()
    data = read_table()
    figure = make_figure(data)
    output = save_figure(figure, args.outdir, args.pdf, args.emit_vector)
    print(output)


if __name__ == "__main__":
    main()
