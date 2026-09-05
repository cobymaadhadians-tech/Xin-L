#!/usr/bin/env python3
"""Render Supplementary Figure 4 from the archived SCZ covariance tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = REPOSITORY_ROOT / "data/supplementary_tables"
FIG_WIDTH_MM = 180
FIG_WIDTH_IN = 7.0866  # 180 mm
FIG_HEIGHT_IN = 5.5118  # 140 mm
FIGURE_STEM = "Supplementary_Figure_4_SCZ_R2_SHAPLEY"

PGC = "#5B7891"
NARROW = "#2A9994"
BROAD = "#D47F6F"
PRIMARY = "#5B7891"
SENSITIVITY = "#8E9BAA"
PTSD = "#C15B5A"
GRID = "#D9D9D9"
TEXT = "#111111"


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "axes.titleweight": "bold",
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


def read_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    full = pd.read_csv(TABLE_DIR / "Supplementary_Table_23_SCZ_downstream_R2_and_Shapley_full.tsv", sep="\t")
    total = pd.read_csv(TABLE_DIR / "Supplementary_Table_24_SCZ_downstream_total_R2_contrasts.tsv", sep="\t")
    omnibus = pd.read_csv(TABLE_DIR / "Supplementary_Table_25_SCZ_downstream_shapley_omnibus.tsv", sep="\t")
    component = pd.read_csv(TABLE_DIR / "Supplementary_Table_26_SCZ_downstream_component_contrasts.tsv", sep="\t")
    return full, total, omnibus, component


def validate_inputs(full: pd.DataFrame, total: pd.DataFrame, omnibus: pd.DataFrame, component: pd.DataFrame) -> None:
    expected_sources = {"PGC", "narrow", "broad"}
    observed_sources = set(full.loc[full["model"].eq("primary_six"), "source"])
    if observed_sources != expected_sources:
        raise ValueError(f"Unexpected primary six-predictor sources: {observed_sources}")

    r2 = (
        full.loc[full["model"].eq("primary_six")]
        .groupby("source", as_index=True)["full_R2"]
        .nunique()
    )
    if (r2 != 1).any():
        raise ValueError("Each primary source must have one full-model R2")

    broad_narrow_total = total[
        total["model"].eq("primary_six") & total["comparison"].eq("broad_minus_narrow")
    ]
    if len(broad_narrow_total) != 1:
        raise ValueError("Missing primary broad-minus-narrow total-R2 contrast")

    broad_narrow_component = component[
        component["comparison"].eq("broad_minus_narrow")
        & component["model"].isin(["primary_six", "sensitivity_all_ten"])
    ]
    if set(broad_narrow_component["model"]) != {"primary_six", "sensitivity_all_ten"}:
        raise ValueError("Missing broad-minus-narrow component contrasts")

    broad_narrow_omnibus = omnibus[
        omnibus["comparison"].eq("broad_minus_narrow")
        & omnibus["test"].eq("raw_phi_full_rank")
        & omnibus["model"].isin(["primary_six", "sensitivity_all_ten"])
    ]
    if len(broad_narrow_omnibus) != 2:
        raise ValueError("Missing broad-minus-narrow absolute Shapley omnibus tests")

    q_values = omnibus.loc[
        omnibus["test"].eq("raw_phi_full_rank")
        & omnibus["model"].isin(["primary_six", "sensitivity_all_ten"]),
        "BH_q_within_model_test",
    ].to_numpy(dtype=float)
    if not np.isfinite(q_values).all() or not (q_values > 0).all() or not (q_values <= 1).all():
        raise ValueError("All q-values used for the strictly positive -log10 transformation must be finite and within (0, 1]")


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.tick_params(colors=TEXT, pad=2)
    ax.xaxis.label.set_color(TEXT)
    ax.yaxis.label.set_color(TEXT)
    ax.set_axisbelow(True)


def add_panel_label(fig: plt.Figure, ax: plt.Axes, label: str, title: str) -> None:
    box = ax.get_position()
    fig.text(box.x0, box.y1 + 0.060, f"({label})", fontsize=8, fontweight="bold", color=TEXT, va="bottom")
    fig.text(box.x0 + 0.035, box.y1 + 0.060, title, fontsize=7, fontweight="bold", color=TEXT, va="bottom")


def format_q(q: float) -> str:
    if q < 0.001:
        return f"{q:.2e}"
    return f"{q:.5f}"


def make_figure(full: pd.DataFrame, total: pd.DataFrame, omnibus: pd.DataFrame, component: pd.DataFrame) -> plt.Figure:
    primary_full = full[full["model"].eq("primary_six")].drop_duplicates("source").copy()
    source_order = ["PGC", "narrow", "broad"]
    source_labels = ["PGC\nSCZ", "narrow\nFinnGen\nSCZ", "broad\nFinnGen\nSCZ"]
    source_colors = [PGC, NARROW, BROAD]
    primary_full["source"] = pd.Categorical(primary_full["source"], categories=source_order, ordered=True)
    primary_full = primary_full.sort_values("source")

    six = component[
        component["model"].eq("primary_six") & component["comparison"].eq("broad_minus_narrow")
    ].copy()
    ten = component[
        component["model"].eq("sensitivity_all_ten") & component["comparison"].eq("broad_minus_narrow")
    ].copy()
    six = six.sort_values("delta_phi_left_minus_right", ascending=True)
    ten = ten.sort_values("delta_phi_left_minus_right", ascending=True)
    predictor_labels = {
        "MDD_CLIN_PGC2025": "MDD",
        "BD_PGC2021": "BD",
        "PTSD": "PTSD",
        "ANX": "ANX",
        "ADHD": "ADHD",
        "ASD": "ASD",
        "OCD_2025": "OCD",
        "AN": "AN",
        "AUD": "AUD",
        "CUD_2023_EUR": "CUD",
    }
    six["predictor"] = six["predictor"].map(predictor_labels).fillna(six["predictor"])
    ten["predictor"] = ten["predictor"].map(predictor_labels).fillna(ten["predictor"])

    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
    gs = fig.add_gridspec(
        1,
        3,
        left=0.08,
        right=0.985,
        top=0.90,
        bottom=0.22,
        width_ratios=[0.92, 1.02, 2.05],
        wspace=0.50,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    cgs = gs[0, 2].subgridspec(2, 1, height_ratios=[6, 10], hspace=0.28)
    ax_c_primary = fig.add_subplot(cgs[0, 0])
    ax_c_sensitivity = fig.add_subplot(cgs[1, 0], sharex=ax_c_primary)

    # Panel a: full-matrix estimates from the prespecified six-predictor model.
    x = np.arange(len(primary_full))
    values = 100 * primary_full["full_R2"].to_numpy()
    ax_a.scatter(x, values, s=34, c=source_colors, edgecolors="white", linewidths=0.7, zorder=3)
    for xx, value in zip(x, values):
        label_x = xx + 0.08 if xx == 0 else xx
        ax_a.text(
            label_x,
            value + 1.3,
            f"{value:.2f}%",
            ha="left" if xx == 0 else "center",
            va="bottom",
            fontsize=6.5,
            color=TEXT,
        )
    ax_a.set_xticks(x, source_labels)
    ax_a.set_ylim(0, 72)
    ax_a.set_yticks([0, 20, 40, 60])
    ax_a.set_ylabel("Joint genetic R² (%)")
    ax_a.grid(axis="y", color=GRID, lw=1.0, zorder=0)
    ax_a.text(
        0.02,
        0.99,
        "Six traits",
        transform=ax_a.transAxes,
        ha="left",
        va="top",
        fontsize=6.3,
        color="#444444",
    )
    clean_axis(ax_a)

    omnibus_plot = omnibus[
        omnibus["test"].eq("raw_phi_full_rank")
        & omnibus["model"].isin(["primary_six", "sensitivity_all_ten"])
    ].copy()
    comparison_order = ["narrow_minus_PGC", "broad_minus_PGC", "broad_minus_narrow"]
    comparison_labels = ["narrow − PGC", "broad − PGC", "broad − narrow"]
    comparison_y = np.arange(len(comparison_order))
    for offset, model, color, label in [
        (-0.13, "primary_six", PRIMARY, "six traits"),
        (0.13, "sensitivity_all_ten", SENSITIVITY, "ten traits"),
    ]:
        subset = omnibus_plot[omnibus_plot["model"].eq(model)].set_index("comparison").reindex(comparison_order)
        q_values = subset["BH_q_within_model_test"].to_numpy(dtype=float)
        x_values = -np.log10(q_values)
        ax_b.scatter(x_values, comparison_y + offset, s=26, color=color, edgecolors="white", linewidths=0.6, zorder=3, label=label)
        for xv, yv, qv in zip(x_values, comparison_y + offset, q_values):
            ax_b.text(xv + 0.14, yv, f"{format_q(qv)}", fontsize=5.2, va="center", color="#444444")
    ax_b.axvline(-np.log10(0.05), color="#777777", lw=1.0, ls=(0, (3, 3)), zorder=0)
    ax_b.set_yticks(comparison_y, comparison_labels)
    ax_b.set_xlim(0, max(-np.log10(omnibus_plot["BH_q_within_model_test"])) + 1.0)
    ax_b.set_xlabel("−log10(q)")
    ax_b.text(1.0, 1.02, "dashed line: q = 0.05", transform=ax_b.transAxes, ha="right", va="bottom", fontsize=5.4, color="#444444")
    legend_b_handles, legend_b_labels = ax_b.get_legend_handles_labels()
    clean_axis(ax_b)

    def plot_component_axis(ax: plt.Axes, data: pd.DataFrame, model_label: str, omnibus_q: float) -> None:
        y = np.arange(len(data))
        values = 100 * data["delta_phi_left_minus_right"].to_numpy()
        lows = 100 * data["CI95_low"].to_numpy()
        highs = 100 * data["CI95_high"].to_numpy()
        colors = [PTSD if p == "PTSD" else (PRIMARY if model_label.startswith("Primary") else SENSITIVITY) for p in data["predictor"]]
        filled = data["BH_q_within_comparison_model"].to_numpy() < 0.05
        ax.axvline(0, color="#777777", lw=1.0, ls=(0, (3, 3)), zorder=0)
        for yi, value, low, high, color, is_filled in zip(y, values, lows, highs, colors, filled):
            ax.errorbar(
                value,
                yi,
                xerr=np.array([[value - low], [high - value]]),
                fmt="o",
                ms=4.8 if color == PTSD else 4.2,
                mfc=color if is_filled or color == PTSD else "white",
                mec=color,
                mew=1.0,
                ecolor=color,
                elinewidth=1.0,
                capsize=0,
                zorder=3,
            )
        ax.set_yticks(y, data["predictor"].tolist())
        for tick, predictor in zip(ax.get_yticklabels(), data["predictor"]):
            if predictor == "PTSD":
                tick.set_fontweight("bold")
                tick.set_color(PTSD)
        ax.set_xlim(-8, 22)
        ax.set_xticks([-5, 0, 5, 10, 15, 20])
        ax.text(
            0.0,
            1.02,
            model_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=6.5,
            fontweight="bold",
            color=TEXT,
        )
        ax.text(
            1.0,
            1.02,
            f"absolute-vector omnibus q = {format_q(omnibus_q)}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=5.6,
            color="#444444",
        )
        clean_axis(ax)

    six_q = float(
        omnibus.loc[
            omnibus["model"].eq("primary_six")
            & omnibus["comparison"].eq("broad_minus_narrow")
            & omnibus["test"].eq("raw_phi_full_rank"),
            "BH_q_within_model_test",
        ].iloc[0]
    )
    ten_q = float(
        omnibus.loc[
            omnibus["model"].eq("sensitivity_all_ten")
            & omnibus["comparison"].eq("broad_minus_narrow")
            & omnibus["test"].eq("raw_phi_full_rank"),
            "BH_q_within_model_test",
        ].iloc[0]
    )
    plot_component_axis(ax_c_primary, six, "Six traits", six_q)
    plot_component_axis(ax_c_sensitivity, ten, "Ten traits", ten_q)
    ax_c_sensitivity.set_xlabel("Broad − narrow Shapley difference (percentage points of R²)")
    ax_c_primary.tick_params(labelbottom=False)
    ax_c_primary.set_xlabel("")

    handles = [
        Line2D([0], [0], marker="o", color="#444444", mfc="#444444", lw=0, ms=4.5, label="BH-FDR < 0.05"),
        Line2D([0], [0], marker="o", color="#444444", mfc="white", lw=0, ms=4.5, label="BH-FDR ≥ 0.05"),
        Line2D([0], [0], marker="o", color=PTSD, mfc=PTSD, lw=0, ms=4.8, label="PTSD"),
    ]
    box_b = ax_b.get_position()
    box_c = ax_c_sensitivity.get_position()
    legend_y = 0.035
    callout_fontsize = 6.2
    fig.legend(
        handles=legend_b_handles,
        labels=legend_b_labels,
        loc="lower center",
        bbox_to_anchor=((box_b.x0 + box_b.x1) / 2, legend_y),
        ncol=2,
        handlelength=0.8,
        columnspacing=1.0,
        fontsize=callout_fontsize,
        frameon=False,
    )
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=((box_c.x0 + box_c.x1) / 2, legend_y),
        ncol=3,
        handlelength=0.8,
        columnspacing=1.0,
        fontsize=callout_fontsize,
        frameon=False,
    )
    add_panel_label(fig, ax_a, "a", "Joint genetic R²")
    add_panel_label(fig, ax_b, "b", "Shapley omnibus q values")
    add_panel_label(fig, ax_c_primary, "c", "Absolute Shapley contrasts")
    return fig


def save_tiff(fig: plt.Figure, outdir: Path, emit_vector: bool = False) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    temporary = outdir / f".{FIGURE_STEM}.rgba.tif"
    final = outdir / f"{FIGURE_STEM}.tif"
    fig.savefig(temporary, format="tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "raw"})
    with Image.open(temporary) as rendered:
        rendered.convert("RGB").save(final, format="TIFF", compression="raw", dpi=(600, 600))
    temporary.unlink()
    if emit_vector:
        fig.savefig(outdir / f".{FIGURE_STEM}.svg", format="svg", bbox_inches="tight")
        fig.savefig(outdir / f".{FIGURE_STEM}.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)
    return final


def save_pdf(fig: plt.Figure, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    final = outdir / f"{FIGURE_STEM}.pdf"
    fig.savefig(final, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return final


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("submission_package/03_figures/supplementary"))
    parser.add_argument("--emit-vector", action="store_true", help="Write hidden SVG/PDF copies for QA")
    parser.add_argument("--pdf", action="store_true", help="Write the user-facing figure as a vector PDF")
    args = parser.parse_args()
    configure_style()
    full, total, omnibus, component = read_tables()
    validate_inputs(full, total, omnibus, component)
    figure = make_figure(full, total, omnibus, component)
    output = save_pdf(figure, args.outdir) if args.pdf else save_tiff(figure, args.outdir, emit_vector=args.emit_vector)
    print(output)


if __name__ == "__main__":
    main()
