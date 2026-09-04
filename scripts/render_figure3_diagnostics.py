#!/usr/bin/env python3
"""Render the SCZ local-covariance diagnostic figure."""

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
DERIVED = REPOSITORY_ROOT / "data/derived"
MM = 1 / 25.4
FIG_WIDTH_MM = 180
FIG_HEIGHT_MM = 150
SCZ = "#3B1F6F"
TEAL = "#2A9994"
CORAL = "#C15B5A"
ORANGE = "#D47F6F"
INDIGO = "#3B1F6F"
GRID = "#D9D9D9"


def style() -> None:
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
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def clean(ax: plt.Axes) -> None:
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.tick_params(colors="#111111", pad=2)
    ax.set_axisbelow(True)


def label_panel(fig: plt.Figure, ax: plt.Axes, letter: str, title: str) -> None:
    ax.text(-0.12, 1.14, f"({letter})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8, fontweight="bold")
    ax.text(-0.01, 1.14, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=7, fontweight="bold")


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.svg", format="svg", bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.pdf", format="pdf", bbox_inches="tight")
    rgba = output_dir / f".{stem}.rgba.tif"
    fig.savefig(rgba, dpi=600, format="tiff", pil_kwargs={"compression": "raw"}, bbox_inches="tight")
    with Image.open(rgba) as im:
        im.convert("RGB").save(output_dir / f"{stem}.tif", format="TIFF", compression="raw", dpi=(600, 600))
    rgba.unlink()
    plt.close(fig)


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, float, float]:
    supergnova = DERIVED / "supergnova"
    p1 = pd.read_csv(supergnova / "pgc_finngen_scz_local_covariance.tsv", sep="\t")
    prof = pd.read_csv(supergnova / "six_trait/sixtrait_profile_summary.tsv", sep="\t")
    panel = pd.read_csv(supergnova / "six_trait/panel_standardized_profile_similarity.tsv", sep="\t")
    adhd = pd.read_csv(supergnova / "scz_adhd_local_covariance.tsv", sep="\t")
    adhd_profile = pd.read_csv(supergnova / "scz_adhd_profile_comparison.tsv", sep="\t")
    rg = pd.read_csv(DERIVED / "phbc/integrated_cross_definition_rg.tsv", sep="\t")
    delta = pd.read_csv(DERIVED / "phbc/integrated_paired_delta_phbc.tsv", sep="\t")
    rg_value = float(rg.loc[rg.comparison.eq("SCZ_PGC2022__SCZ_FINNGEN_R13"), "rg"].iloc[0])
    delta_value = float(delta.loc[delta.comparison.eq("SCZ_FINNGEN_minus_PGC"), "delta_pp"].iloc[0])
    return p1, prof, panel, adhd, adhd_profile, rg_value, delta_value


def draw_summary(ax: plt.Axes, rg_value: float, delta_value: float) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rows = [(0.68, r"$r_g$ (PGC–FinnGen SCZ)", rg_value, (0.70, 1.01), [0.8, 0.9, 1.0], "#3B1F6F"),
            (0.29, r"ΔPHBC (FinnGen − PGC; pp)", delta_value, (-40, 5), [-40, -20, 0], SCZ)]
    for y, title, value, limits, ticks, color in rows:
        left, right = 0.24, 0.96
        ax.text(0.02, y + 0.07, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.2, color="#222222")
        x0, x1 = limits
        ax.plot([left, right], [y, y], transform=ax.transAxes, color="#555555", lw=1.0, clip_on=False)
        for tick in ticks:
            xx = left + (tick - x0) / (x1 - x0) * (right - left)
            ax.plot([xx, xx], [y - 0.018, y + 0.018], transform=ax.transAxes, color="#777777", lw=1.0, clip_on=False)
            ax.text(xx, y - 0.055, f"{tick:g}", transform=ax.transAxes, ha="center", va="top", fontsize=5.8, color="#444444")
        xx = left + (value - x0) / (x1 - x0) * (right - left)
        ax.plot(xx, y, "o", transform=ax.transAxes, ms=6.4, color=color, zorder=3, clip_on=False)
        ax.text(
            0.98,
            y + 0.025,
            f"{value:.3f}" if value < 1 else f"{value:.2f}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=6.2,
            color="#444444",
        )
    ax.text(0.02, 0.08, "Pairwise genome-wide similarity does not determine the fixed-panel multivariate contrast.", transform=ax.transAxes, ha="left", va="bottom", fontsize=5.9, color="#444444", wrap=True)


def draw_local_genome(ax: plt.Axes, p1: pd.DataFrame) -> None:
    p1 = p1.sort_values(["chr", "start"]).copy()
    p1["chr"] = p1["chr"].astype(int)
    chromosome_max = p1.groupby("chr")["end"].max().to_dict()
    p1["x"] = p1.apply(
        lambda r: (int(r.chr) - 1) + ((r.start + r.end) / 2) / chromosome_max[int(r.chr)],
        axis=1,
    )
    colors = np.where(p1.rho >= 0, TEAL, CORAL)
    ax.scatter(p1.x, p1.rho * 1000, s=7, c=colors, alpha=0.52, linewidths=0, rasterized=True, zorder=2)
    sig = p1[p1.bonferroni.astype(str).str.lower().eq("true")]
    ax.scatter(sig.x, sig.rho * 1000, s=28, facecolor=ORANGE, edgecolor=INDIGO, linewidth=0.65, zorder=4)
    ax.axhline(0, color="#777777", lw=1.0, zorder=1)
    ax.set_xticks(np.arange(22) + 0.5, [str(i) for i in range(1, 23)], fontsize=5.4)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"Local covariance, $\rho$ ($\times 10^{-3}$)")
    positive_significant = int((sig.rho >= 0).sum())
    significance_text = f"{positive_significant} positive" if positive_significant != len(sig) else "all positive"
    ax.text(0.99, 0.97, f"K = {len(p1):,}; {len(sig)} Bonferroni-positive blocks\n{significance_text}", transform=ax.transAxes, ha="right", va="top", fontsize=5.9, color="#333333")
    clean(ax)
    ax.grid(axis="x", visible=False)


def draw_profile_similarity(ax: plt.Axes, prof: pd.DataFrame, panel: pd.DataFrame) -> None:
    order = ["ADHD", "ASD", "OCD", "AN", "AUD", "CUD"]
    d = prof.set_index("auxiliary").loc[order].reset_index()
    y = np.arange(len(d))
    ax.axvline(0, color="#777777", lw=1.0, ls=(0, (3, 3)), zorder=0)
    ax.hlines(y, 0, d.pearson_r, color="#C7C7C7", lw=1.1, zorder=1)
    ax.scatter(d.pearson_r, y, s=32, color=TEAL, edgecolor="#1B5E5A", linewidth=0.5, zorder=3)
    panel_r = float(panel.panel_standardized_profile_pearson_r.iloc[0])
    ax.axvline(panel_r, color=INDIGO, lw=1.0, ls=(0, (4, 3)), zorder=2)
    common_blocks = int(panel.common_blocks_all_six.iloc[0])
    trait_count = int(panel.traits.iloc[0])
    ax.text(0.99, 1.03, f"panel-wide r = {panel_r:.3f}\n{common_blocks:,} common blocks × {trait_count} traits", transform=ax.transAxes, color="#333333", ha="right", va="bottom", fontsize=5.8)
    ax.set_yticks(y, order)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.42)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_xlabel("Pearson correlation of local covariance profiles")
    ax.text(0.01, -0.23, "Auxiliary-specific values are descriptive; the panel-wide value uses standardized local ρ vectors.", transform=ax.transAxes, ha="left", va="top", fontsize=5.6, color="#444444")
    clean(ax)
    ax.grid(axis="x", visible=False)


def draw_adhd_scatter(ax: plt.Axes, adhd: pd.DataFrame, profile: pd.DataFrame) -> None:
    x = adhd.rho_PGC_ADHD.to_numpy(float) * 1000
    y = adhd.rho_FG_ADHD.to_numpy(float) * 1000
    same = np.sign(x) == np.sign(y)
    ax.scatter(x[same], y[same], s=8, color=TEAL, alpha=0.46, linewidths=0, rasterized=True)
    ax.scatter(x[~same], y[~same], s=8, color=CORAL, alpha=0.46, linewidths=0, rasterized=True)
    lim = max(np.abs(np.r_[x, y]).max(), 0.8) * 1.06
    ax.plot([-lim, lim], [-lim, lim], color="#555555", lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.axhline(0, color=GRID, lw=0.6, zorder=0)
    ax.axvline(0, color=GRID, lw=0.6, zorder=0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"PGC SCZ × ADHD local $\rho$ ($\times 10^{-3}$)")
    ax.set_ylabel(r"FinnGen SCZ × ADHD local $\rho$ ($\times 10^{-3}$)")
    summary = profile.iloc[0]
    ax.text(
        0.03,
        0.97,
        f"Pearson r = {summary.pearson_r:.3f}\n"
        f"Informative sign concordance = {100 * summary.informative_sign_concordance:.1f}%\n"
        f"Top-5% overlap = {int(summary.top5pct_intersection_n)}/{int(summary.top5pct_n_each)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        color="#333333",
    )
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=TEAL, markeredgewidth=0, markersize=4.5, label="same sign"), Line2D([0], [0], marker="o", color="none", markerfacecolor=CORAL, markeredgewidth=0, markersize=4.5, label="opposite sign"), Line2D([0], [0], color="#555555", lw=1.0, ls=(0, (4, 3)), label="$y=x$")], loc="lower right", fontsize=5.7, handlelength=1.4, borderaxespad=0.3)
    clean(ax)
    ax.grid(False)


def render_figure3(output_dir: Path) -> None:
    style()
    p1, prof, panel, adhd, adhd_profile, rg_value, delta_value = load_data()
    fig = plt.figure(figsize=(FIG_WIDTH_MM * MM, FIG_HEIGHT_MM * MM), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.095, right=0.985, top=0.91, bottom=0.12, hspace=0.60, wspace=0.34, height_ratios=[0.92, 1.08])
    a, b, c, d = [fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    draw_summary(a, rg_value, delta_value)
    draw_local_genome(b, p1)
    draw_profile_similarity(c, prof, panel)
    draw_adhd_scatter(d, adhd, adhd_profile)
    label_panel(fig, a, "a", "Global concordance and PHBC contrast")
    label_panel(fig, b, "b", "PGC–FinnGen SCZ local covariance")
    label_panel(fig, c, "c", "Six-trait local-profile similarity")
    label_panel(fig, d, "d", "ADHD local-covariance profile")
    stem = output_dir / "Figure_3_SCZ_DIAGNOSTICS"
    save_figure(fig, output_dir, stem.name)
    print(f"wrote {stem}.tif")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    render_figure3(args.output_dir)


if __name__ == "__main__":
    main()
