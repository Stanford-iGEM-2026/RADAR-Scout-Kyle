"""Reproduce the iGEM keloid sensor nomination through RADAR-Scout.

Uses the Deng 2021 fibroblasts with the AUTHORS' Seurat cluster labels (cluster 3 =
pathological MFB) and compares keloid MFB vs normal-scar fibroblasts, ranking by the
detection-difference specificity score. This is the comparison that surfaces the
near-binary keloid markers ADAM12 (92% vs 3%) and POSTN (100% vs 21%).

Usage: python scripts/reproduce_deng_mfb.py [path/to/deng_fibroblasts.h5ad]
"""

from __future__ import annotations

import sys

import anndata as ad
import numpy as np
import scipy.sparse as sp

from radar_scout.scoring import score_matrix
from radar_scout.genesets import filter_technical

H5 = sys.argv[1] if len(sys.argv) > 1 else "/Users/rohk/iGEM/data/processed/scRNA/deng_fibroblasts.h5ad"
MFB_CLUSTER = "3"

a = ad.read_h5ad(H5)
cond = a.obs["condition"].astype(str).to_numpy()
clus = a.obs["seurat_clusters"].astype(str).to_numpy()

# P = keloid pathological MFB (author cluster 3); H = normal-scar fibroblasts (any cluster)
pop = np.full(a.n_obs, "", dtype=object)
pop[(cond == "Keloid") & (clus == MFB_CLUSTER)] = "P"
pop[cond == "Normal scar"] = "H"
sel = pop != ""
a, pop = a[sel].copy(), pop[sel]
donor = a.obs["sample"].astype(str).to_numpy()

# X is log1p(CP10k) -> back to linear CP10k for the scorer
X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
X = np.expm1(X)
genes = a.var_names.to_numpy().astype(str)

# candidate genes: detected in >=10% of pathogenic MFB cells, minus technical genes
detect = (X[pop == "P"] > 0).mean(axis=0)
keep = (detect >= 0.10) & filter_technical(genes)
Xc, g = X[:, keep], genes[keep]

df = score_matrix(Xc, g, donor, pop, pos_label="P")
df = df[df["spec_score"].notna()].sort_values("spec_score", ascending=False).reset_index(drop=True)

print(f"P (keloid MFB C{MFB_CLUSTER}): {int((pop=='P').sum())} cells / {len(set(donor[pop=='P']))} donors | "
      f"H (normal scar fib): {int((pop=='H').sum())} cells / {len(set(donor[pop=='H']))} donors | genes: {len(g)}")
print("\n=== TOP 20 by detection-specificity (RADAR-Scout, Deng MFB vs normal scar) ===")
cols = ["gene", "spec_score", "log2FC", "detect_P", "detect_H", "delta_detect", "RACS", "DSS"]
print(df.head(20)[[c for c in cols if c in df.columns]].round(2).to_string(index=False))
print("\n=== canonical sensor candidates ===")
for m in ["ADAM12", "POSTN", "ASPN", "COL11A1", "SDC1", "LOXL2", "TGFBI", "TNC"]:
    r = df.index[df["gene"] == m]
    print(f"  {m:8s} rank {int(r[0]) + 1}" if len(r) else f"  {m:8s} (filtered/not detected)")

df.to_parquet("outputs/keloid_deng_mfb_racs.parquet")
print("\nwrote outputs/keloid_deng_mfb_racs.parquet")
