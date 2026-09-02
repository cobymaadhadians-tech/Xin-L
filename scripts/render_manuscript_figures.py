#!/usr/bin/env python3
"""Render the manuscript figures from derived TSV files."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DERIVED = REPOSITORY_ROOT / "data/derived"
MM = 1 / 25.4

TEAL = "#2A9994"
CORAL = "#C15B5A"
BLUE = "#648CB5"
ORANGE = "#D47F6F"
INDIGO = "#3B1F6F"

COLORS = {
    "MDD": BLUE,
    "BD": ORANGE,
    "SCZ": INDIGO,
    "AD": TEAL,
    "Epilepsy": CORAL,
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.titlesize": 7,
            "axes.titleweight": "bold",
            "axes.labelsize": 7,
            "axes.linewidth": 1.0,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 6.2,
            "ytick.labelsize": 6.2,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "legend.fontsize": 6.0,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_color("#222222")
    ax.tick_params(colors="#111111", pad=2)
    ax.xaxis.label.set_color("#111111")
    ax.yaxis.label.set_color("#111111")


def panel_label(fig: plt.Figure, ax: plt.Axes, label: str, dx: float = -0.035, dy: float = 0.018) -> None:
    box = ax.get_position()
    fig.text(box.x0 + dx, box.y1 + dy, f"({label})", fontsize=8, fontweight="bold", va="bottom", ha="left")


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    rgba_tiff = out_dir / f".{stem}.rgba.tif"
    fig.savefig(
        rgba_tiff,
        format="tiff",
        dpi=600,
        pil_kwargs={"compression": "raw"},
    )
    with Image.open(rgba_tiff) as rendered:
        rendered.convert("RGB").save(
            out_dir / f"{stem}.tif",
            format="TIFF",
            compression="raw",
            dpi=(600, 600),
        )
    rgba_tiff.unlink()
    plt.close(fig)


def plot_herror(ax: plt.Axes, x, y, low, high, color, marker="o", filled=True, alpha=1.0, ms=5.5, lw=1.0, zorder=3):
    face = color if filled else "white"
    ax.errorbar(
        x,
        y,
        xerr=np.array([[x - low], [high - x]]),
        fmt=marker,
        ms=ms,
        mfc=face,
        mec=color,
        mew=1.0,
        ecolor=color,
        elinewidth=lw,
        capsize=0,
        alpha=alpha,
        zorder=zorder,
    )


def figure1(out_dir: Path, rg: pd.DataFrame, delta: pd.DataFrame) -> None:
    rows = [
        ("MDD EHR - Clinical", "MDD", "MDD_CLIN_PGC2025__MDD_EHR_PGC2025", "MDD_EHR_minus_CLIN"),
        ("MDD FinnGen - Clinical", "MDD", "MDD_CLIN_PGC2025__MDD_FINNGEN_R13", "MDD_FINNGEN_minus_CLIN"),
        ("MDD Questionnaire - Clinical", "MDD", "MDD_CLIN_PGC2025__MDD_QUEST_PGC2025", "MDD_QUEST_minus_CLIN"),
        ("MDD EHR - Questionnaire", "MDD", "MDD_EHR_PGC2025__MDD_QUEST_PGC2025", "MDD_EHR_minus_QUEST"),
        ("MDD EHR - FinnGen", "MDD", "MDD_EHR_PGC2025__MDD_FINNGEN_R13", "MDD_EHR_minus_FINNGEN"),
        ("MDD FinnGen - Questionnaire", "MDD", "MDD_QUEST_PGC2025__MDD_FINNGEN_R13", "MDD_FINNGEN_minus_QUEST"),
        ("BD FinnGen - Clinical", "BD", "BD_CLIN_PGC4__BD_FINNGEN_R13", "BD_FINNGEN_minus_CLIN"),
        ("SCZ narrow FinnGen - PGC", "SCZ", "SCZ_PGC2022__SCZ_FINNGEN_R13", "SCZ_FINNGEN_minus_PGC"),
        ("AD No-proxy - Main", "AD", "AD_GCST90704646_MAIN__AD_GCST90704647_NOPROXY", "AD_NOPROXY_minus_MAIN"),
        ("AD No-biobank - Main", "AD", "AD_GCST90704646_MAIN__AD_GCST90704648_NOBIOBANK", "AD_NOBIOBANK_minus_MAIN"),
        ("AD No-proxy - No-biobank", "AD", "AD_GCST90704647_NOPROXY__AD_GCST90704648_NOBIOBANK", "AD_NOPROXY_minus_NOBIOBANK"),
        ("Epilepsy EHR2 - ILAE", "Epilepsy", "ILAE_vs_EHR_META_2OF3", "EPILEPSY_EHR2_minus_ILAE"),
        ("Epilepsy EHR3 - ILAE", "Epilepsy", "ILAE_vs_EHR_META_3OF3", "EPILEPSY_EHR3_minus_ILAE"),
        ("Epilepsy FinnGen - ILAE", "Epilepsy", "EPILEPSY_ILAE2023_EUR__EPILEPSY_FINNGEN_R13", "EPILEPSY_FINNGEN_minus_ILAE"),
    ]
    fm = pd.DataFrame(rows, columns=["label", "disease", "rg_key", "delta_key"])
    fm = fm.merge(rg[["comparison", "rg", "rg_se"]], left_on="rg_key", right_on="comparison")
    fm = fm.merge(delta[["comparison", "delta_pp", "ci95_low_pp", "ci95_high_pp"]], left_on="delta_key", right_on="comparison")
    broad = read_tsv(DERIVED / "phbc/broad_scz_figure1.tsv")
    fm = pd.concat([fm, broad], ignore_index=True, sort=False)
    core = {
        "MDD EHR - Questionnaire",
        "BD FinnGen - Clinical",
        "SCZ narrow FinnGen - PGC",
        "SCZ broad FinnGen - PGC",
        "SCZ broad - narrow FinnGen",
        "AD No-proxy - Main",
        "Epilepsy FinnGen - ILAE",
    }
    fig, ax = plt.subplots(figsize=(180 * MM, 112 * MM))
    fig.subplots_adjust(left=0.12, right=0.985, top=0.96, bottom=0.25)
    ax.axhline(0, color="#777777", lw=1.0, ls=(0, (4, 4)), zorder=0)
    ax.axvline(1, color="#777777", lw=1.0, ls=(0, (4, 4)), zorder=0)
    for _, r in fm.iterrows():
        is_core = r.label in core
        alpha = 0.9 if is_core else 0.32
        marker = "^" if r.disease in {"AD", "Epilepsy"} else "o"
        ax.errorbar(
            r.rg,
            r.delta_pp,
            xerr=1.96 * r.rg_se,
            yerr=np.array([[r.delta_pp - r.ci95_low_pp], [r.ci95_high_pp - r.delta_pp]]),
            fmt=marker,
            ms=6.5 if is_core else 4.2,
            mfc="white",
            mec=COLORS[r.disease],
            mew=1.0,
            ecolor=COLORS[r.disease],
            elinewidth=1.2 if is_core else 1.0,
            alpha=alpha,
            capsize=0,
            zorder=3 if is_core else 1,
        )
    offsets = {
        "MDD EHR - Questionnaire": (10, 18),
        "BD FinnGen - Clinical": (-28, 10),
        "SCZ narrow FinnGen - PGC": (0, -16),
        "SCZ broad FinnGen - PGC": (-8, -17),
        "SCZ broad - narrow FinnGen": (6, 13),
        "AD No-proxy - Main": (-48, -12),
        "Epilepsy FinnGen - ILAE": (-8, 12),
    }
    aligns = {
        "MDD EHR - Questionnaire": "left",
        "BD FinnGen - Clinical": "right",
        "AD No-proxy - Main": "left",
        "Epilepsy FinnGen - ILAE": "center",
        "SCZ broad FinnGen - PGC": "right",
        "SCZ broad - narrow FinnGen": "left",
    }
    for _, r in fm[fm.label.isin(core)].iterrows():
        ax.annotate(
            r.label,
            (r.rg, r.delta_pp),
            xytext=offsets[r.label],
            textcoords="offset points",
            ha=aligns.get(r.label, "center"),
            va="center",
            fontsize=6.5,
            color="#222222",
            arrowprops=dict(arrowstyle="-", color="#999999", lw=1.0, shrinkA=2, shrinkB=4),
            zorder=5,
        )
    ax.set(xlim=(0.68, 1.12), ylim=(-42, 48), xlabel=r"Cross-definition $r_g$", ylabel="ΔPHBC (percentage points)")
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1])
    ax.set_yticks([-40, -20, 0, 20, 40])
    clean_axis(ax)
    disease_handles = [Line2D([0], [0], marker="o", color=COLORS[k], lw=1.1, ms=5, label=k) for k in ["AD", "BD", "Epilepsy", "MDD", "SCZ"]]
    shape_handles = [
        Line2D([0], [0], marker="^", color="#111111", mfc="white", lw=0, ms=6, label="Neurological"),
        Line2D([0], [0], marker="o", color="#111111", mfc="white", lw=0, ms=6, label="Psychiatric"),
    ]
    leg1 = fig.legend(disease_handles, [h.get_label() for h in disease_handles], loc="lower center", bbox_to_anchor=(0.5, 0.095), ncol=5, handlelength=1.0, columnspacing=0.9)
    fig.add_artist(leg1)
    fig.legend(shape_handles, [h.get_label() for h in shape_handles], loc="lower center", bbox_to_anchor=(0.5, 0.035), ncol=2, handlelength=1.0, columnspacing=1.2)
    save_figure(fig, out_dir, "Figure_1_RG_PHBC_DISSOCIATION")


def psych_label(target: str) -> tuple[str, str, str]:
    mapping = {
        "MDD_CLIN_PGC2025": ("MDD", "Clinical", "MDD  Clinical"),
        "MDD_EHR_PGC2025": ("MDD", "EHR", "MDD  EHR"),
        "MDD_QUEST_PGC2025": ("MDD", "Questionnaire", "MDD  Questionnaire"),
        "MDD_FINNGEN_R13": ("MDD", "FinnGen", "MDD  FinnGen"),
        "BD_CLIN_PGC4": ("BD", "Clinical", "BD  Clinical"),
        "BD_FINNGEN_R13": ("BD", "FinnGen", "BD  FinnGen"),
        "SCZ_PGC2022": ("SCZ", "PGC", "SCZ  PGC"),
        "SCZ_FINNGEN_R13": ("SCZ", "FinnGen", "SCZ  narrow FinnGen"),
    }
    return mapping[target]


def figure2(out_dir: Path, phbc: pd.DataFrame, delta: pd.DataFrame) -> None:
    psych = phbc[phbc.section.eq("psychiatric")].copy()
    psych[["disease", "definition", "display"]] = psych.target.apply(lambda x: pd.Series(psych_label(x)))
    order = ["MDD  Clinical", "MDD  EHR", "MDD  Questionnaire", "MDD  FinnGen", "BD  Clinical", "BD  FinnGen", "SCZ  PGC", "SCZ  narrow FinnGen"]
    labels = {
        "MDD_EHR_minus_CLIN": "MDD EHR - Clinical",
        "MDD_FINNGEN_minus_CLIN": "MDD FinnGen - Clinical",
        "MDD_QUEST_minus_CLIN": "MDD Questionnaire - Clinical",
        "MDD_EHR_minus_QUEST": "MDD EHR - Questionnaire",
        "MDD_EHR_minus_FINNGEN": "MDD EHR - FinnGen",
        "MDD_FINNGEN_minus_QUEST": "MDD FinnGen - Questionnaire",
        "BD_FINNGEN_minus_CLIN": "BD FinnGen - Clinical",
        "SCZ_FINNGEN_minus_PGC": "SCZ narrow FinnGen - PGC",
    }
    d = delta[delta.section.eq("psychiatric")].copy()
    d["display"] = d.comparison.map(labels)
    fig, (a, b) = plt.subplots(1, 2, figsize=(180 * MM, 136 * MM), gridspec_kw={"width_ratios": [1.0, 1.12]})
    fig.subplots_adjust(left=0.14, right=0.98, top=0.89, bottom=0.22, wspace=0.62)
    shape = {"Clinical": "o", "PGC": "s", "EHR": "s", "Questionnaire": "^", "FinnGen": "o"}
    filled = {"Clinical": False, "PGC": False, "EHR": True, "Questionnaire": False, "FinnGen": True}
    for y, name in enumerate(order[::-1]):
        r = psych.loc[psych.display.eq(name)].iloc[0]
        low = 100 * (r.phbc - 1.96 * r.official_se)
        high = 100 * (r.phbc + 1.96 * r.official_se)
        plot_herror(a, r.phbc_pp, y, low, high, COLORS[r.disease], shape[r.definition], filled[r.definition])
    a.set_yticks(range(len(order)), order[::-1])
    a.set(xlim=(-3, 72), xticks=[0, 20, 40, 60], xlabel="PHBC (estimate and 95% CI)", title="Psychiatric PHBC estimates")
    a.set_xticklabels(["0%", "20%", "40%", "60%"])
    for y, name in enumerate(list(labels.values())[::-1]):
        r = d.loc[d.display.eq(name)].iloc[0]
        plot_herror(b, r.delta_pp, y, r.ci95_low_pp, r.ci95_high_pp, COLORS[r.disease], "o", bool(r.fdr_pass))
    b.axvline(0, color="#777777", lw=1.0, ls=(0, (4, 4)), zorder=0)
    b.set_yticks(range(len(labels)), list(labels.values())[::-1])
    b.set(xlim=(-45, 32), xticks=[-40, -20, 0, 20], xlabel="ΔPHBC (percentage points; 95% CI)", title="Complete paired contrasts")
    for ax in (a, b):
        clean_axis(ax)
    panel_label(fig, a, "a")
    panel_label(fig, b, "b")
    disease_handles = [Line2D([0], [0], marker="o", color=COLORS[k], lw=0, ms=5.5, label=k) for k in ["BD", "MDD", "SCZ"]]
    def_handles = [
        Line2D([0], [0], marker="o", color="#111111", mfc="white", lw=0, ms=5.5, label="Clinical"),
        Line2D([0], [0], marker="s", color="#111111", mfc="white", lw=0, ms=5.5, label="PGC"),
        Line2D([0], [0], marker="o", color="#111111", mfc="#111111", lw=0, ms=5.5, label="FinnGen"),
        Line2D([0], [0], marker="s", color="#111111", mfc="#111111", lw=0, ms=5.5, label="EHR"),
        Line2D([0], [0], marker="^", color="#111111", mfc="white", lw=0, ms=6, label="Questionnaire"),
    ]
    sig_handles = [
        Line2D([0], [0], marker="o", color="#111111", mfc="#111111", lw=0, ms=5.5, label="BH-FDR < 0.05"),
        Line2D([0], [0], marker="o", color="#111111", mfc="white", lw=0, ms=5.5, label="BH-FDR ≥ 0.05"),
    ]
    fig.legend(handles=disease_handles, loc="lower center", bbox_to_anchor=(0.5, 0.105), ncol=3, handlelength=0.8, columnspacing=0.9)
    fig.legend(handles=def_handles, loc="lower center", bbox_to_anchor=(0.30, 0.035), ncol=3, handlelength=0.8, columnspacing=0.8)
    fig.legend(handles=sig_handles, loc="lower center", bbox_to_anchor=(0.76, 0.035), ncol=1, handlelength=0.8, columnspacing=0.8)
    save_figure(fig, out_dir, "Figure_2_PSYCHIATRIC_PHBC")


def figure3(out_dir: Path) -> None:
    try:
        from render_figure3_diagnostics import render_figure3
    except ImportError:
        from scripts.render_figure3_diagnostics import render_figure3

    render_figure3(out_dir)


def strip(ax: plt.Axes, text: str) -> None:
    ax.add_patch(Rectangle((0, 1.0), 1, 0.085, transform=ax.transAxes, clip_on=False, facecolor="#F2F2F2", edgecolor="none", zorder=0))
    ax.text(0.5, 1.042, text, transform=ax.transAxes, ha="center", va="center", fontsize=6.5, fontweight="bold")


def supplementary3(out_dir: Path, loo: pd.DataFrame) -> None:
    bs = loo[loo.disease.isin(["BD", "SCZ"])].copy()
    label_common = {"internalizing": "Internalizing", "neurodevelopmental": "Neurodevelopmental", "compulsive_eating": "Compulsive / eating", "substance": "Substance"}
    def domain_label(r):
        if r.domain == "mood_psychotic":
            return "Mood / psychotic (SCZ)" if r.disease == "BD" else "Mood / psychotic (MDD + BD)"
        return label_common[r.domain]
    bs["label"] = bs.apply(domain_label, axis=1)
    fig = plt.figure(figsize=(180 * MM, 140 * MM))
    gs = fig.add_gridspec(2, 2, left=0.19, right=0.98, top=0.83, bottom=0.18, hspace=0.26, wspace=0.78)
    axes = {(d, p): fig.add_subplot(gs[i, j]) for i, d in enumerate(["BD", "SCZ"]) for j, p in enumerate(["raw", "diff"])}
    for disease in ["BD", "SCZ"]:
        data = bs[bs.disease.eq(disease)].copy()
        order = list(data.label)
        ar, ad = axes[(disease, "raw")], axes[(disease, "diff")]
        ar.axvline(0, color="#777777", lw=1.0, ls=(0, (4, 4)), zorder=0)
        ad.axvline(0, color="#777777", lw=1.0, ls=(0, (4, 4)), zorder=0)
        for y, (_, r) in enumerate(data.iloc[::-1].iterrows()):
            ar.plot([r.left_reduction_pp, r.right_reduction_pp], [y + 0.07, y - 0.07], color="#BDBDBD", lw=1.0, zorder=1)
            ar.plot(r.left_reduction_pp, y + 0.07, "D", ms=4.6, mfc="white", mec=COLORS[disease], mew=1.0, zorder=2)
            ar.plot(r.right_reduction_pp, y - 0.07, "o", ms=4.8, mfc="white", mec=COLORS[disease], mew=1.0, zorder=2)
            plot_herror(ad, r.difference_in_reduction_pp, y, r.ci95_low_pp, r.ci95_high_pp, COLORS[disease], filled=bool(r.survives_bh_0_05), ms=4.8)
        ar.set_yticks(range(len(order)), order[::-1])
        ad.set_yticks(range(len(order)), order[::-1])
        ar.set(xlim=(-11, 33), xticks=[-10, 0, 10, 20, 30])
        ad.set(xlim=(-36, 14), xticks=[-30, -20, -10, 0, 10])
        strip(ar, disease)
        strip(ad, f"{disease}: {'FinnGen' if disease == 'BD' else 'narrow FinnGen'} - {'Clinical' if disease == 'BD' else 'PGC'}")
        for ax in (ar, ad):
            clean_axis(ax)
    axes[("SCZ", "raw")].set_xlabel("Full-panel-scaled reduction (percentage points)")
    axes[("SCZ", "diff")].set_xlabel("Difference in reduction (pp; 95% CI)")
    fig.text(0.15, 0.965, "(a)", fontsize=8, fontweight="bold", va="top")
    fig.text(0.19, 0.965, "Category-removal effects", fontsize=7, fontweight="bold", va="top")
    fig.text(0.19, 0.925, "Negative values indicate higher PHBC after category removal", fontsize=6.2, color="#444444", va="top")
    fig.text(0.61, 0.965, "(b)", fontsize=8, fontweight="bold", va="top")
    fig.text(0.65, 0.965, "Paired differences", fontsize=7, fontweight="bold", va="top")
    raw_handles = [Line2D([0], [0], marker="o", color="#111111", mfc="white", lw=0, ms=5.5, label="Consortium (Clinical/PGC)"), Line2D([0], [0], marker="D", color="#111111", mfc="white", lw=0, ms=5.2, label="FinnGen (narrow for SCZ)")]
    sig_handles = [Line2D([0], [0], marker="o", color="#111111", mfc="#111111", lw=0, ms=5.5, label="BH-FDR < 0.05"), Line2D([0], [0], marker="o", color="#111111", mfc="white", lw=0, ms=5.5, label="BH-FDR ≥ 0.05")]
    fig.legend(handles=raw_handles, loc="lower left", bbox_to_anchor=(0.20, 0.045), ncol=2, handlelength=0.8, columnspacing=0.9)
    fig.legend(handles=sig_handles, loc="lower right", bbox_to_anchor=(0.97, 0.045), ncol=2, handlelength=0.8, columnspacing=0.9)
    save_figure(fig, out_dir, "Supplementary_Figure_3_CATEGORY_REMOVAL")


def supplementary1(out_dir: Path, rg: pd.DataFrame, phbc: pd.DataFrame, delta: pd.DataFrame) -> None:
    nrg = rg[rg.section.eq("neurological")].copy()
    rg_labels = {
        "AD_GCST90704646_MAIN__AD_GCST90704647_NOPROXY": "AD Main / No-proxy",
        "AD_GCST90704646_MAIN__AD_GCST90704648_NOBIOBANK": "AD Main / No-biobank",
        "AD_GCST90704647_NOPROXY__AD_GCST90704648_NOBIOBANK": "AD No-proxy / No-biobank",
        "ILAE_vs_EHR_META_2OF3": "Epilepsy ILAE / EHR meta 2/3",
        "ILAE_vs_EHR_META_3OF3": "Epilepsy ILAE / EHR meta 3/3",
        "EPILEPSY_ILAE2023_EUR__EPILEPSY_FINNGEN_R13": "Epilepsy ILAE / FinnGen",
    }
    nrg["display"] = nrg.comparison.map(rg_labels)
    rg_order = list(rg_labels.values())
    npd = phbc[phbc.section.eq("neurological")].copy()
    definitions = {
        "AD_GCST90704646_MAIN": ("AD", "Main", "AD  Main"), "AD_GCST90704647_NOPROXY": ("AD", "No-proxy", "AD  No-proxy"), "AD_GCST90704648_NOBIOBANK": ("AD", "No-biobank", "AD  No-biobank"),
        "EPILEPSY_ILAE2023_EUR": ("Epilepsy", "ILAE", "Epilepsy  ILAE"), "EPILEPSY_FINNGEN_R13": ("Epilepsy", "FinnGen", "Epilepsy  FinnGen"), "EPILEPSY_EHR_META_2OF3": ("Epilepsy", "EHR meta 2/3", "Epilepsy  EHR meta 2/3"), "EPILEPSY_EHR_META_3OF3": ("Epilepsy", "EHR meta 3/3", "Epilepsy  EHR meta 3/3"),
    }
    npd[["disease", "definition", "display"]] = npd.target.apply(lambda x: pd.Series(definitions[x]))
    p_order = [v[2] for v in definitions.values()]
    nd = delta[(delta.section.eq("neurological")) & (~delta.comparison.eq("EPILEPSY_EHR_META_2OF3_minus_3OF3"))].copy()
    dlabels = {"AD_NOPROXY_minus_MAIN": "AD No-proxy - Main", "AD_NOBIOBANK_minus_MAIN": "AD No-biobank - Main", "AD_NOPROXY_minus_NOBIOBANK": "AD No-proxy - No-biobank", "EPILEPSY_EHR2_minus_ILAE": "Epilepsy EHR2 - ILAE", "EPILEPSY_EHR3_minus_ILAE": "Epilepsy EHR3 - ILAE", "EPILEPSY_FINNGEN_minus_ILAE": "Epilepsy FinnGen - ILAE"}
    nd["display"] = nd.comparison.map(dlabels)
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 136 * MM), gridspec_kw={"width_ratios": [1.08, 0.95, 1.08]})
    fig.subplots_adjust(left=0.215, right=0.985, top=0.88, bottom=0.22, wspace=0.90)
    a, b, c = axes
    a.axvline(1, color="#777777", lw=1.0, ls=(0, (4, 4)), zorder=0)
    for y, lab in enumerate(rg_order[::-1]):
        r = nrg.loc[nrg.display.eq(lab)].iloc[0]
        plot_herror(a, r.rg, y, r.rg - 1.96 * r.rg_se, r.rg + 1.96 * r.rg_se, COLORS[r.disease], ms=4.8)
    a.set_yticks(range(len(rg_order)), rg_order[::-1])
    a.set(xlim=(0.52, 1.06), xticks=[0.6, 0.8, 1.0], xlabel=r"$r_g$ (95% CI)", title=r"Cross-definition $r_g$")
    markers = {"Main": ("o", False), "No-proxy": ("^", False), "No-biobank": ("s", False), "ILAE": ("D", False), "FinnGen": ("o", True), "EHR meta 2/3": ("s", True), "EHR meta 3/3": ("^", True)}
    for y, lab in enumerate(p_order[::-1]):
        r = npd.loc[npd.display.eq(lab)].iloc[0]
        m, fill = markers[r.definition]
        plot_herror(b, r.phbc_pp, y, 100 * (r.phbc - 1.96 * r.official_se), 100 * (r.phbc + 1.96 * r.official_se), COLORS[r.disease], m, fill, ms=4.8)
    b.set_yticks(range(len(p_order)), p_order[::-1])
    b.set(xlim=(-1, 16), xticks=[0, 5, 10, 15], xlabel="PHBC (estimate and 95% CI)", title="PHBC estimates")
    b.set_xticklabels(["0%", "5%", "10%", "15%"])
    c.axvline(0, color="#777777", lw=1.0, ls=(0, (4, 4)), zorder=0)
    for y, lab in enumerate(list(dlabels.values())[::-1]):
        r = nd.loc[nd.display.eq(lab)].iloc[0]
        plot_herror(c, r.delta_pp, y, r.ci95_low_pp, r.ci95_high_pp, COLORS[r.disease], ms=4.8)
    c.set_yticks(range(len(dlabels)), list(dlabels.values())[::-1])
    c.set(xlim=(-7, 7.5), xticks=[-5, 0, 5], xlabel="ΔPHBC (pp; 95% CI)", title="Paired ΔPHBC")
    for ax in axes:
        clean_axis(ax)
        ax.tick_params(axis="y", labelsize=6.3)
    for ax, lab in zip(axes, "abc"):
        panel_label(fig, ax, lab)
    disease_handles = [Line2D([0], [0], marker="o", color=COLORS[k], lw=1, ms=4.8, label=k) for k in ["AD", "Epilepsy"]]
    shape_handles = [Line2D([0], [0], marker=m, color="#111111", mfc=("#111111" if f else "white"), lw=0, ms=5.5, label=l) for l, (m, f) in markers.items()]
    fig.legend(handles=disease_handles, loc="lower center", bbox_to_anchor=(0.5, 0.08), ncol=2, handlelength=1, columnspacing=0.8)
    fig.legend(handles=shape_handles, loc="lower center", bbox_to_anchor=(0.5, 0.02), ncol=4, handlelength=0.8, columnspacing=0.7)
    save_figure(fig, out_dir, "Supplementary_Figure_1_NEUROLOGICAL_BOUNDARY")


def supplementary2(out_dir: Path, leave_one: pd.DataFrame) -> None:
    comparisons = [
        ("SCZ_FINNGEN_minus_PGC", "SCZ: narrow FinnGen – PGC", CORAL),
        ("MDD_EHR_minus_QUEST", "MDD: EHR – questionnaire", TEAL),
        ("MDD_FINNGEN_minus_QUEST", "MDD: FinnGen – questionnaire", TEAL),
    ]
    aux_labels = {
        "MDD_CLIN_PGC2025": "MDD", "BD_PGC2021": "BD", "SCZ_PGC2022": "SCZ", "PTSD": "PTSD", "ANX": "ANX",
        "ADHD": "ADHD", "ASD": "ASD", "OCD_2025": "OCD", "AN": "AN", "AUD": "AUD",
        "CUD_2023_EUR": "CUD",
    }
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 120 * MM))
    fig.subplots_adjust(left=0.10, right=0.985, top=0.88, bottom=0.18, wspace=0.62)
    for ax, (comparison, title, color), label in zip(axes, comparisons, "abc"):
        data = leave_one[leave_one.comparison.eq(comparison)].copy()
        full = data[data.condition.eq("full_panel")].iloc[0]
        reduced = data[data.condition.eq("leave_one_auxiliary_out")].copy()
        reduced["label"] = reduced.omitted_auxiliary.map(aux_labels)
        top_to_bottom = (["MDD", "BD", "PTSD", "ANX", "ADHD", "ASD", "OCD", "AN", "AUD", "CUD"] if comparison == "SCZ_FINNGEN_minus_PGC" else ["BD", "SCZ", "ADHD", "ASD", "OCD", "AN", "AUD", "CUD"])
        reduced = reduced.set_index("label").loc[top_to_bottom[::-1]].reset_index()
        ax.axvline(0, color="#777777", lw=1.0, ls=(0, (2, 3)), zorder=0)
        ax.axvline(full.delta_pp, color=color, lw=1.0, ls=(0, (5, 3)), zorder=0)
        for y, (_, row) in enumerate(reduced.iterrows()):
            plot_herror(ax, row.delta_pp, y, row.ci95_low_pp, row.ci95_high_pp, color, ms=4.8)
        full_y = len(reduced)
        plot_herror(ax, full.delta_pp, full_y, full.ci95_low_pp, full.ci95_high_pp, color, marker="D", ms=5.0)
        ax.set_yticks(range(full_y + 1), list(reduced.label) + ["Full panel"])
        ax.set_xlabel("ΔPHBC (pp; 95% CI)")
        ax.set_title(title)
        clean_axis(ax)
        panel_label(fig, ax, label, dx=-0.025)
    save_figure(fig, out_dir, "Supplementary_Figure_2_LEAVE_ONE_AUXILIARY_OUT")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=REPOSITORY_ROOT / "outputs")
    args = parser.parse_args()
    configure_style()
    main_dir = args.output_dir / "main"
    supplementary_dir = args.output_dir / "supplementary"
    rg = read_tsv(DERIVED / "phbc/integrated_cross_definition_rg.tsv")
    phbc = read_tsv(DERIVED / "phbc/integrated_target_phbc.tsv")
    delta = read_tsv(DERIVED / "phbc/integrated_paired_delta_phbc.tsv")
    loo = read_tsv(DERIVED / "phbc/leave_category_out_paired_did.tsv")
    leave_one = read_tsv(DERIVED / "phbc/leave_one_auxiliary_out_paired_phbc.tsv")
    figure1(main_dir, rg, delta)
    figure2(main_dir, phbc, delta)
    figure3(main_dir)
    supplementary1(supplementary_dir, rg, phbc, delta)
    supplementary2(supplementary_dir, leave_one)
    supplementary3(supplementary_dir, loo)


if __name__ == "__main__":
    main()
