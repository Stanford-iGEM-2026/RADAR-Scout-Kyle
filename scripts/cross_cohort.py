"""Cross-cohort validation (spec section 6): rank stability & consensus across
independent cohorts of the same disease.

Usage:
    python scripts/cross_cohort.py CELLxGENE=outputs/keloid_subpop_racs.parquet \
        GEO=outputs/keloid_geo_racs.parquet --rank-by DSS --outdir figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CRIMSON, TEAL = "#8e1918", "#1c7170"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})


def load(spec, rank_by):
    label, path = spec.split("=", 1)
    df = pd.read_parquet(path)
    if rank_by == "DSS" and "DSS" not in df.columns and {"log2FC", "mean_P"}.issubset(df.columns):
        df["DSS"] = df["log2FC"].clip(lower=0) * np.log10(df["mean_P"].clip(lower=0) + 1)
    df = df[df[rank_by].notna()].copy()
    df["rank"] = df[rank_by].rank(ascending=False, method="min")
    df["pct"] = 100 * (1 - df["rank"] / len(df))  # percentile (100 = best)
    return label, df.set_index("gene")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cohorts", nargs="+", help="LABEL=path.parquet")
    ap.add_argument("--rank-by", default="DSS")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    cohorts = dict(load(s, args.rank_by) for s in args.cohorts)
    labels = list(cohorts)
    shared = set.intersection(*[set(df.index) for df in cohorts.values()])
    print(f"shared genes across {len(labels)} cohorts: {len(shared)}")

    pct = pd.DataFrame({lab: cohorts[lab].loc[list(shared), "pct"] for lab in labels})
    pct["consensus"] = pct[labels].mean(axis=1)
    pct["min_pct"] = pct[labels].min(axis=1)  # robust = high in the WORST cohort
    pct = pct.sort_values("consensus", ascending=False)

    # pairwise rank correlation
    print("\n=== cross-cohort Spearman rank correlation (higher = more reproducible) ===")
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            rho = spearmanr(pct[labels[i]], pct[labels[j]]).correlation
            print(f"  {labels[i]} vs {labels[j]}: rho = {rho:.3f}")

    print(f"\n=== top {args.top} CONSENSUS targets (robust across cohorts) ===")
    show = pct.head(args.top).copy()
    show.index.name = "gene"
    print(show.round(1).to_string())

    out_csv = Path(args.outdir).parent / "outputs" / f"cross_cohort_{args.rank_by}.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pct.to_csv(out_csv)

    # scatter of percentiles between the first two cohorts, consensus-colored
    if len(labels) >= 2:
        a, b = labels[0], labels[1]
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(pct[a], pct[b], s=10, c=pct["consensus"], cmap="RdGy_r", linewidth=0, alpha=0.7)
        for g, r in pct.head(12).iterrows():
            ax.annotate(g, (r[a], r[b]), fontsize=7, fontstyle="italic",
                        xytext=(2, 2), textcoords="offset points")
        ax.plot([0, 100], [0, 100], color="#ccc", lw=0.8, ls="--")
        ax.set_xlabel(f"{a} percentile ({args.rank_by})")
        ax.set_ylabel(f"{b} percentile ({args.rank_by})")
        ax.set_title(f"cross-cohort reproducibility (rho={spearmanr(pct[a], pct[b]).correlation:.2f})",
                     fontsize=10)
        fig.tight_layout()
        fig.savefig(Path(args.outdir) / "fig_cross_cohort.png", bbox_inches="tight")
        print(f"\nwrote {out_csv} and {args.outdir}/fig_cross_cohort.png")


if __name__ == "__main__":
    main()
