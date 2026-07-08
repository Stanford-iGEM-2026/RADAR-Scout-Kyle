"""Donor-aware pseudobulk aggregation.

Sums raw counts within each (donor, population) group — the standard, correct
way to prepare single-cell data for donor-level differential expression and
reproducibility analysis without pseudoreplication. The per-group profiles feed
mixed-effects / pseudobulk DE (statsmodels, or an R limma/DESeq2 bridge).
"""

from __future__ import annotations

import numpy as np


def aggregate_pseudobulk(counts, donor, pop, gene_names, min_cells=10):
    """Sum raw counts per (donor, population) group.

    Parameters
    ----------
    counts : (n_cells, n_genes) array-like of RAW counts.
    donor, pop : (n_cells,) group labels.
    gene_names : (n_genes,) names.
    min_cells : drop groups with fewer than this many cells.

    Returns
    -------
    pandas.DataFrame with a MultiIndex (donor, pop), columns = genes, plus an
    ``n_cells`` column. Values are summed raw counts.
    """
    import pandas as pd

    counts = np.asarray(counts, dtype=float)
    donor = np.asarray(donor)
    pop = np.asarray(pop)
    keys = list(zip(donor.tolist(), pop.tolist()))
    uniq = sorted(set(keys))

    records, index, ncells = [], [], []
    for d, p in uniq:
        mask = (donor == d) & (pop == p)
        if mask.sum() < min_cells:
            continue
        records.append(counts[mask].sum(axis=0))
        index.append((d, p))
        ncells.append(int(mask.sum()))

    if not records:
        return pd.DataFrame(columns=list(gene_names) + ["n_cells"])

    df = pd.DataFrame(records, columns=list(gene_names),
                      index=pd.MultiIndex.from_tuples(index, names=["donor", "pop"]))
    df["n_cells"] = ncells
    return df


def cpm_normalize(pseudobulk_df, log=True, pseudocount=1.0):
    """Counts-per-million normalize a pseudobulk DataFrame (excludes ``n_cells``).

    Returns a new DataFrame; optionally log1p-transformed.
    """
    import pandas as pd

    gene_cols = [c for c in pseudobulk_df.columns if c != "n_cells"]
    mat = pseudobulk_df[gene_cols].to_numpy(dtype=float)
    lib = mat.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    cpm = mat / lib * 1e6
    if log:
        cpm = np.log1p(cpm) if pseudocount == 1.0 else np.log(cpm + pseudocount)
    out = pd.DataFrame(cpm, columns=gene_cols, index=pseudobulk_df.index)
    if "n_cells" in pseudobulk_df.columns:
        out["n_cells"] = pseudobulk_df["n_cells"].to_numpy()
    return out
