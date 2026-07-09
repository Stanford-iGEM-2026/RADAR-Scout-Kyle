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


def fig_volcano(de, outdir, fdr_thr=0.05, lfc_thr=1.0, label_top=15):
    """Volcano of donor-level DE (log2FC vs -log10 FDR).

    Falls back to a fold-change-vs-abundance view when donor-level significance is
    underpowered (few donors in a small pathogenic subpopulation → all-NaN p/FDR),
    so the figure is never empty and stays honest about what's estimable.
    """
    d = de[de["log2FC"].notna()].copy()
    ycol = next((c for c in ["FDR", "p_value", "p_mannwhitney"]
                 if c in d and d[c].notna().sum() > 0), None)
    fig, ax = plt.subplots(figsize=(5.0, 4.6))

    if ycol is None:  # underpowered — MA-style fallback
        mcol = "mean_pos" if "mean_pos" in d else ("mean_P" if "mean_P" in d else None)
        x = np.log10((d[mcol] if mcol else 1.0) + 1)
        big = d["log2FC"].abs() > lfc_thr
        ax.scatter(x[~big], d["log2FC"][~big], s=8, color="#cccccc", alpha=0.6, linewidth=0)
        ax.scatter(x[big & (d.log2FC > 0)], d.log2FC[big & (d.log2FC > 0)], s=12, color=CRIMSON, linewidth=0)
        ax.scatter(x[big & (d.log2FC < 0)], d.log2FC[big & (d.log2FC < 0)], s=12, color=TEAL, linewidth=0)
        for _, r in d.reindex(d["log2FC"].abs().sort_values(ascending=False).index).head(label_top).iterrows():
            ax.annotate(r["gene"], (np.log10((r[mcol] if mcol else 1) + 1), r["log2FC"]),
                        fontsize=7, fontstyle="italic", xytext=(2, 2), textcoords="offset points", color="#333")
        ax.axhline(lfc_thr, color="#999", lw=0.8, ls="--")
        ax.axhline(-lfc_thr, color="#999", lw=0.8, ls="--")
        ax.set_xlabel("log10 mean expression (pathogenic)")
        ax.set_ylabel("log2 fold-change vs reference")
        ax.set_title("fold-change vs abundance (donor-level significance underpowered)",
                     fontsize=8, color="#666")
    else:
        d["nlq"] = -np.log10(d[ycol].clip(lower=1e-300))
        d = d[d["nlq"].notna()]
        sig = (d[ycol] < fdr_thr) & (d["log2FC"].abs() > lfc_thr)
        up = sig & (d["log2FC"] > 0)
        ax.scatter(d.loc[~sig, "log2FC"], d.loc[~sig, "nlq"], s=8, color="#cccccc", alpha=0.6, linewidth=0)
        ax.scatter(d.loc[up, "log2FC"], d.loc[up, "nlq"], s=12, color=CRIMSON, linewidth=0)
        ax.scatter(d.loc[sig & ~up, "log2FC"], d.loc[sig & ~up, "nlq"], s=12, color=TEAL, linewidth=0)
        ax.axhline(-np.log10(fdr_thr), color="#999", lw=0.8, ls="--")
        ax.axvline(lfc_thr, color="#999", lw=0.8, ls="--")
        ax.axvline(-lfc_thr, color="#999", lw=0.8, ls="--")
        for _, r in d[up].nlargest(label_top, "nlq").iterrows():
            ax.annotate(r["gene"], (r["log2FC"], r["nlq"]), fontsize=7, fontstyle="italic",
                        xytext=(2, 2), textcoords="offset points", color="#333")
        ax.set_xlabel("log2 fold-change (pathogenic vs reference)")
        ax.set_ylabel(f"-log10 {ycol}")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_volcano.png", bbox_inches="tight")
    plt.close(fig)


def _pops_present(df):
    return [p for p in ("P", "B", "H", "R") if f"act_{p}" in df.columns]


def fig_dotplot(df, outdir, n=20):
    """Dot plot: top genes x populations; dot size = detection %, color = mean expr."""
    pops = _pops_present(df)
    d = df.head(n)[::-1]
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(1.4 + 0.7 * len(pops), 0.32 * len(d) + 1))
    for xi, p in enumerate(pops):
        det = d.get(f"detect_{p}", pd.Series(np.full(len(d), 50.0)))
        mean = d.get(f"mean_{p}", d[f"act_{p}"] * 100)
        c = mean / (mean.max() + 1e-9)
        ax.scatter(np.full(len(d), xi), y, s=6 + 2.2 * det.to_numpy(),
                   c=c, cmap="Reds", edgecolor="#888", linewidth=0.3, vmin=0, vmax=1)
    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels(pops)
    ax.set_yticks(y)
    ax.set_yticklabels(d["gene"])
    _italic_genes(ax)
    ax.set_xlim(-0.5, len(pops) - 0.5)
    ax.set_title("dot size = detection %, color = mean expr", fontsize=8, color="#666")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_dotplot.png", bbox_inches="tight")
    plt.close(fig)


def fig_heatmap(df, outdir, n=25):
    """Heatmap of RADAR activation across populations for the top genes."""
    pops = _pops_present(df)
    d = df.head(n)[::-1]
    mat = d[[f"act_{p}" for p in pops]].to_numpy()
    fig, ax = plt.subplots(figsize=(1.2 + 0.6 * len(pops), 0.28 * len(d) + 1))
    im = ax.imshow(mat, aspect="auto", cmap="Reds", vmin=0, vmax=1)
    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels(pops)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d["gene"])
    _italic_genes(ax)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("RADAR activation", fontsize=8)
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def fig_umap(umap, outdir, gene=None):
    """UMAP panels: population, cell type, and a candidate gene's expression."""
    gene_cols = [c for c in umap.columns if c not in ("UMAP1", "UMAP2", "cell_type", "pop", "donor")]
    gene = gene or (gene_cols[0] if gene_cols else None)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    # (a) population
    palette = {"P": CRIMSON, "B": "#bdbdbd", "H": TEAL, "R": "#e6a23c"}
    for p, sub in umap.groupby("pop"):
        axes[0].scatter(sub["UMAP1"], sub["UMAP2"], s=4, color=palette.get(p, "#999"),
                        label=p, linewidth=0, alpha=0.7)
    axes[0].legend(frameon=False, fontsize=8, markerscale=2)
    axes[0].set_title("population", fontsize=10)
    # (b) cell type
    cts = umap["cell_type"].value_counts().index[:10]
    cmap = plt.get_cmap("tab10")
    for i, ct in enumerate(cts):
        sub = umap[umap["cell_type"] == ct]
        axes[1].scatter(sub["UMAP1"], sub["UMAP2"], s=4, color=cmap(i % 10), label=str(ct)[:18],
                        linewidth=0, alpha=0.7)
    axes[1].legend(frameon=False, fontsize=6, markerscale=2, loc="best")
    axes[1].set_title("cell type", fontsize=10)
    # (c) gene expression
    if gene:
        sc = axes[2].scatter(umap["UMAP1"], umap["UMAP2"], s=4, c=umap[gene], cmap="Reds", linewidth=0)
        fig.colorbar(sc, ax=axes[2], fraction=0.046, pad=0.02).outline.set_visible(False)
        axes[2].set_title(f"{gene} (CP10k)", fontsize=10, fontstyle="italic")
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_umap.png", bbox_inches="tight")
    plt.close(fig)


def fig_violin(umap, outdir, n_genes=4):
    """Violin plots of candidate-gene expression across populations (P/B/H/R)."""
    gene_cols = [c for c in umap.columns if c not in ("UMAP1", "UMAP2", "cell_type", "pop", "donor")][:n_genes]
    pops = [p for p in ["P", "B", "H", "R"] if p in set(umap["pop"])]
    if not gene_cols or not pops:
        return
    palette = {"P": CRIMSON, "B": "#bdbdbd", "H": TEAL, "R": "#e6a23c"}
    fig, axes = plt.subplots(1, len(gene_cols), figsize=(2.4 * len(gene_cols), 3.4))
    axes = np.atleast_1d(axes)
    for ax, g in zip(axes, gene_cols):
        data = [umap.loc[umap["pop"] == p, g].to_numpy() for p in pops]
        parts = ax.violinplot(data, showmeans=False, showextrema=False, widths=0.85)
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(palette.get(pops[i], "#999"))
            pc.set_alpha(0.8)
            pc.set_edgecolor("white")
        ax.set_xticks(range(1, len(pops) + 1))
        ax.set_xticklabels(pops)
        ax.set_title(g, fontsize=10, fontstyle="italic")
        ax.set_ylabel("expression (CP10k)" if ax is axes[0] else "")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_violin.png", bbox_inches="tight")
    plt.close(fig)


def fig_forest(umap, outdir, gene=None):
    """Per-donor mean expression (pathogenic vs off-target donors) with 95% CI —
    a donor-level effect-size forest for one candidate gene."""
    gene_cols = [c for c in umap.columns if c not in ("UMAP1", "UMAP2", "cell_type", "pop", "donor")]
    gene = gene or (gene_cols[0] if gene_cols else None)
    if gene is None:
        return
    d = umap[["donor", "pop", gene]].copy()
    d["grp"] = np.where(d["pop"] == "P", "pathogenic", "off-target")
    per = d.groupby(["grp", "donor"])[gene].mean().reset_index()
    fig, ax = plt.subplots(figsize=(4.8, 0.28 * len(per) + 1.2))
    y, yticks, ylabels = 0, [], []
    for grp, color in [("pathogenic", CRIMSON), ("off-target", "#9e9e9e")]:
        sub = per[per["grp"] == grp].sort_values(gene)
        for _, r in sub.iterrows():
            ax.plot(r[gene], y, "o", color=color, ms=4, alpha=0.7)
            yticks.append(y)
            ylabels.append(str(r["donor"])[:12])
            y += 1
        v = sub[gene].to_numpy()
        if len(v) >= 2:
            m, ci = v.mean(), 1.96 * v.std(ddof=1) / np.sqrt(len(v))
            ax.plot([m - ci, m + ci], [y, y], color=color, lw=2.5)
            ax.plot(m, y, "D", color=color, ms=7)
            yticks.append(y)
            ylabels.append(f"{grp} mean ± CI")
            y += 1
        y += 0.6
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel(f"{gene} per-donor mean (CP10k)")
    ax.set_title(f"{gene} — donor-level effect", fontsize=10, fontstyle="italic")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_forest.png", bbox_inches="tight")
    plt.close(fig)


def fig_paga(paga, outdir):
    """PAGA graph: cluster nodes at their UMAP centroid, sized by cell count, colored
    by dominant population; edge width = connectivity (differentiation relatedness)."""
    nodes, edges = paga.get("nodes", []), paga.get("edges", [])
    if not nodes:
        return
    pos = {n["cluster"]: (n["x"], n["y"]) for n in nodes}
    palette = {"P": CRIMSON, "B": "#bdbdbd", "H": TEAL, "R": "#e6a23c"}
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for e in edges:
        if e["a"] in pos and e["b"] in pos:
            (x1, y1), (x2, y2) = pos[e["a"]], pos[e["b"]]
            ax.plot([x1, x2], [y1, y2], color="#cccccc", lw=0.5 + 4 * e["w"], alpha=0.6, zorder=1)
    for n in nodes:
        dom = max(n["pop"], key=n["pop"].get) if n["pop"] else "?"
        ax.scatter(n["x"], n["y"], s=40 + n["n"] / 15, color=palette.get(dom, "#999"),
                   edgecolor="white", linewidth=1.2, zorder=2)
        ax.annotate(n["cluster"], (n["x"], n["y"]), fontsize=7, ha="center", va="center", zorder=3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("PAGA — cell-state connectivity (node color = dominant population)", fontsize=9)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "fig_paga.png", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parquet")
    ap.add_argument("--meta", default=None)
    ap.add_argument("--de", default=None)
    ap.add_argument("--umap", default=None)
    ap.add_argument("--paga", default=None)
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet(args.parquet).sort_values("RACS", ascending=False).reset_index(drop=True)
    df = df[df["RACS"].notna()]

    n = 4
    fig_racs_ranking(df, args.outdir)
    fig_knee(df, args.outdir)
    fig_components(df, args.outdir)
    fig_window(df, args.outdir)
    fig_dotplot(df, args.outdir)
    fig_heatmap(df, args.outdir)
    n += 2
    if args.de and Path(args.de).exists():
        fig_volcano(pd.read_parquet(args.de), args.outdir); n += 1
    if args.umap and Path(args.umap).exists():
        um = pd.read_parquet(args.umap)
        fig_umap(um, args.outdir)
        fig_violin(um, args.outdir)
        fig_forest(um, args.outdir, gene=str(df.iloc[0]["gene"]) if len(df) else None)
        n += 3
    if args.paga and Path(args.paga).exists():
        fig_paga(json.load(open(args.paga)), args.outdir); n += 1

    title = ""
    if args.meta and Path(args.meta).exists():
        m = json.load(open(args.meta))
        title = f"{m.get('disease','')} / {','.join(m.get('pathogenic_cell_types', []))}"
    print(f"wrote {n} figures to {args.outdir}/  [{title}]  (n={len(df)} genes)")


if __name__ == "__main__":
    main()
