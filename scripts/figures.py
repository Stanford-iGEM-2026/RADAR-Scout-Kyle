"""Publication-quality RADAR-Scout figures from a RACS results table.

Nature/Cell-clean: brand palette, no bold, thin spines, no chartjunk. Reads the
parquet written by modal_app build_and_score (+ optional meta.json).

Usage:
    python scripts/figures.py outputs/keloid_skinfibro_racs.parquet \
        --meta outputs/keloid_skinfibro_meta.json --outdir figures
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CRIMSON = "#8e1918"
TEAL = "#1c7170"
GREY = "#9e9e9e"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.weight": "normal",
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
    "figure.dpi": 300,
})


def _italic_genes(ax, axis="y"):
    labels = ax.get_yticklabels() if axis == "y" else ax.get_xticklabels()
    for lab in labels:
        lab.set_fontstyle("italic")


def fig_racs_ranking(df, outdir, n=20):
    d = df.head(n)[::-1]
    fig, ax = plt.subplots(figsize=(4.8, 6.2))
    ax.barh(d["gene"], d["RACS"], color=CRIMSON, height=0.7)
    ax.set_xlabel("RADAR Activation Compatibility Score")
    ax.set_xlim(0, max(0.05, df["RACS"].max() * 1.08))
    ax.tick_params(labelsize=9)
    _italic_genes(ax)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_racs_ranking.png", bbox_inches="tight")
    plt.close(fig)


def fig_knee(df, outdir, label_top=12):
    """Abundance (Feas) vs specificity (Sep): the RADAR trade-off / 'knee'."""
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    sizes = 20 + 240 * df["RACS"].clip(lower=0)
    sc = ax.scatter(df["Feas"], df["Sep"], s=sizes, c=df["RACS"], cmap="RdGy_r",
                    edgecolor="white", linewidth=0.4, alpha=0.9, vmin=0)
    top = df.head(label_top)
    for _, r in top.iterrows():
        ax.annotate(r["gene"], (r["Feas"], r["Sep"]), fontsize=7.5, fontstyle="italic",
                    xytext=(3, 3), textcoords="offset points", color="#333")
    ax.set_xlabel("Feasibility  (P(activation) at reachable threshold)")
    ax.set_ylabel("Specificity  (donor-level AUC, P vs off-target)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0.4, 1.02)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("RACS", fontsize=9)
    cb.outline.set_visible(False)
    ax.text(0.97, 0.46, "RADAR sweet spot →\nhigh abundance + high specificity",
            ha="right", va="bottom", fontsize=8, color=TEAL)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_knee_abundance_specificity.png", bbox_inches="tight")
    plt.close(fig)


def fig_components(df, outdir, n=12):
    d = df.head(n)[::-1]
    comps = [("Sep", CRIMSON), ("Feas", "#c0653f"), ("Repro", TEAL), ("1-OffMax", "#5aa9a3")]
    vals = {"Sep": d["Sep"], "Feas": d["Feas"], "Repro": d["Repro"], "1-OffMax": 1 - d["OffMax"]}
    y = np.arange(len(d))
    h = 0.2
    fig, ax = plt.subplots(figsize=(5.6, 6.0))
    for i, (c, col) in enumerate(comps):
        ax.barh(y + (i - 1.5) * h, vals[c], height=h, color=col, label=c)
    ax.set_yticks(y)
    ax.set_yticklabels(d["gene"])
    ax.set_xlabel("component score")
    ax.set_xlim(0, 1)
    ax.tick_params(labelsize=9)
    _italic_genes(ax)
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="lower right")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_racs_components.png", bbox_inches="tight")
    plt.close(fig)


def fig_window(df, outdir, n=12):
    """On-target vs worst off-target activation — the therapeutic window."""
    d = df.head(n)[::-1]
    off = d[[c for c in ("act_H", "act_R") if c in d.columns]].max(axis=1)
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(5.2, 6.0))
    ax.barh(y + 0.2, d.get("act_P", d["Feas"]), height=0.38, color=CRIMSON, label="pathogenic (on-target)")
    ax.barh(y - 0.2, off, height=0.38, color=GREY, label="worst off-target")
    ax.set_yticks(y)
    ax.set_yticklabels(d["gene"])
    ax.set_xlabel("RADAR activation")
    ax.set_xlim(0, 1)
    ax.tick_params(labelsize=9)
    _italic_genes(ax)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_therapeutic_window.png", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parquet")
    ap.add_argument("--meta", default=None)
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet(args.parquet).sort_values("RACS", ascending=False).reset_index(drop=True)
    df = df[df["RACS"].notna()]

    fig_racs_ranking(df, args.outdir)
    fig_knee(df, args.outdir)
    fig_components(df, args.outdir)
    fig_window(df, args.outdir)

    title = ""
    if args.meta and Path(args.meta).exists():
        m = json.load(open(args.meta))
        title = f"{m.get('disease','')} / {','.join(m.get('pathogenic_cell_types', []))}"
    print(f"wrote 4 figures to {args.outdir}/  [{title}]  (n={len(df)} genes)")


if __name__ == "__main__":
    main()
