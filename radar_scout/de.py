"""Donor-aware differential expression.

Single-cell target prioritization must use the DONOR (patient), not the cell, as
the experimental unit: the thousands of cells from one individual are technical
pseudoreplicates, not independent biological replicates, and treating them as
such inflates significance by orders of magnitude (see docs/RACS_framework.md and
Squair et al. 2021, "Confronting false discoveries in single-cell differential
expression"). This module offers two donor-aware DE routes:

* ``pseudobulk_de`` — the recommended default. Collapse each donor to one
  pseudobulk profile (sum raw counts, then CPM+log), then run an ordinary
  donor-level two-sample test per gene. n = number of donors.
* ``mixedlm_de`` — a cell-level linear mixed model with a per-donor random
  intercept, for a single gene, as a sensitivity check on the pseudobulk call.

``forest_data`` prepares the per-donor points and a bootstrap CI behind a forest
plot for one gene.

Everything is array-based (per-cell ``counts``/``expr``, ``donor``, ``condition``
labels) so the statistics are testable without single-cell dependencies, matching
``scoring.py`` and ``pseudobulk.py``.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu

from .pseudobulk import aggregate_pseudobulk, cpm_normalize

# Small constant added to linear CPM means before the log2 ratio, so a gene that
# is zero in one group gives a large-but-finite fold change rather than +/-inf.
_EPS = 1.0


def _neg_labels(neg_labels, condition, pos_label):
    """Off-target labels: explicit ``neg_labels`` or every condition but ``pos_label``."""
    if neg_labels is not None:
        return list(neg_labels)
    return [c for c in np.unique(condition) if c != pos_label]


def _bh_fdr(pvals):
    """Benjamini-Hochberg FDR, NaN-safe. NaN p-values map to NaN q-values.

    Ranks only the finite p-values, applies the step-up q = p * m / rank, and
    enforces monotonicity. ``m`` is the number of finite p-values.
    """
    pvals = np.asarray(pvals, dtype=float)
    q = np.full(pvals.shape, np.nan)
    finite = np.where(np.isfinite(pvals))[0]
    m = finite.size
    if m == 0:
        return q
    p = pvals[finite]
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adj = ranked * m / (np.arange(1, m + 1))
    # enforce monotone non-decreasing q from the largest p downward, then clip
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    q_finite = np.empty(m)
    q_finite[order] = adj
    q[finite] = q_finite
    return q


def _pseudobulk_matrices(counts, donor, condition, gene_names, pos_label,
                         neg_labels=None, min_cells=10):
    """Aggregate to per-donor log-CPM pseudobulk and split donors by condition.

    Because donors are nested within condition (each donor is one condition), the
    (donor, condition) pseudobulk index has exactly one row per donor. Returns
    ``(gene_names, pos_log, neg_log, pos_lin, neg_lin, pos_donors, neg_donors)``
    where the ``*_log`` arrays are (n_donors_group, n_genes) log-CPM and ``*_lin``
    are the corresponding linear CPM (for fold-change on the natural scale).
    """
    condition = np.asarray(condition)
    neg = _neg_labels(neg_labels, condition, pos_label)

    pb = aggregate_pseudobulk(counts, donor, condition, gene_names, min_cells=min_cells)
    gene_cols = [c for c in pb.columns if c != "n_cells"]
    log_cpm = cpm_normalize(pb, log=True)
    lin_cpm = cpm_normalize(pb, log=False)

    # index level 1 is the condition ("pop") passed to aggregate_pseudobulk
    conds = pb.index.get_level_values(1).to_numpy()
    donors = pb.index.get_level_values(0).to_numpy()
    pos_mask = conds == pos_label
    neg_mask = np.isin(conds, neg)

    log_mat = log_cpm[gene_cols].to_numpy(dtype=float)
    lin_mat = lin_cpm[gene_cols].to_numpy(dtype=float)
    return (
        gene_cols,
        log_mat[pos_mask], log_mat[neg_mask],
        lin_mat[pos_mask], lin_mat[neg_mask],
        donors[pos_mask], donors[neg_mask],
    )


def pseudobulk_de(counts, donor, condition, gene_names, pos_label,
                  neg_labels=None, min_cells=10):
    """Donor-level pseudobulk differential expression, pos condition vs pooled off-target.

    Each donor contributes ONE value per gene (its CPM+log pseudobulk), so the
    per-gene tests have n = number of donors — the correct, pseudoreplication-free
    unit. For every gene we run a Welch t-test (unequal variance) and a
    Mann-Whitney U on the donor-level log-CPM, and report the log2 fold change on
    the linear-CPM group means.

    Parameters
    ----------
    counts : (n_cells, n_genes) RAW counts.
    donor : (n_cells,) donor id per cell.
    condition : (n_cells,) condition label per cell (e.g. "P", "H", "R"). Used as
        the ``pop`` for pseudobulk aggregation; donors are nested within condition.
    gene_names : (n_genes,) names.
    pos_label : the on-target condition.
    neg_labels : off-target condition labels; default = all conditions but ``pos_label``.
    min_cells : passed to ``aggregate_pseudobulk`` (drop tiny donor groups).

    Returns
    -------
    pandas.DataFrame with columns
    ``[gene, log2FC, mean_pos, mean_neg, t_stat, p_value, p_mannwhitney, FDR,
    n_pos_donors, n_neg_donors]``, sorted by FDR ascending then |log2FC| descending.
    ``log2FC`` and ``mean_pos``/``mean_neg`` are on linear CPM; the tests are on
    log-CPM. FDR is Benjamini-Hochberg across genes on ``p_value``.
    """
    import pandas as pd

    (genes, pos_log, neg_log, pos_lin, neg_lin,
     pos_donors, neg_donors) = _pseudobulk_matrices(
        counts, donor, condition, gene_names, pos_label, neg_labels, min_cells)

    n_pos = int(pos_log.shape[0])
    n_neg = int(neg_log.shape[0])

    mean_pos_lin = pos_lin.mean(axis=0) if n_pos else np.full(len(genes), np.nan)
    mean_neg_lin = neg_lin.mean(axis=0) if n_neg else np.full(len(genes), np.nan)
    log2fc = np.log2((mean_pos_lin + _EPS) / (mean_neg_lin + _EPS))

    t_stat = np.full(len(genes), np.nan)
    p_value = np.full(len(genes), np.nan)
    p_mw = np.full(len(genes), np.nan)
    # Welch t-test / Mann-Whitney need >= 2 donors per side (and some variance).
    if n_pos >= 2 and n_neg >= 2:
        with np.errstate(invalid="ignore"):
            t_res = ttest_ind(pos_log, neg_log, axis=0, equal_var=False)
        t_stat = np.asarray(t_res.statistic, dtype=float)
        p_value = np.asarray(t_res.pvalue, dtype=float)
        for j in range(len(genes)):
            a, b = pos_log[:, j], neg_log[:, j]
            # mannwhitneyu raises when one side is entirely constant *and* equal
            # to the other; guard so a single degenerate gene doesn't kill the run.
            try:
                p_mw[j] = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
            except ValueError:
                p_mw[j] = np.nan

    fdr = _bh_fdr(p_value)

    df = pd.DataFrame({
        "gene": list(genes),
        "log2FC": log2fc,
        "mean_pos": mean_pos_lin,
        "mean_neg": mean_neg_lin,
        "t_stat": t_stat,
        "p_value": p_value,
        "p_mannwhitney": p_mw,
        "FDR": fdr,
        "n_pos_donors": n_pos,
        "n_neg_donors": n_neg,
    })
    df["_absfc"] = df["log2FC"].abs()
    df = df.sort_values(
        ["FDR", "_absfc"], ascending=[True, False], na_position="last"
    ).drop(columns="_absfc").reset_index(drop=True)
    return df


def mixedlm_de(expr_cells, donor, condition, pos_label, neg_labels=None):
    """Cell-level linear mixed model for one gene: expression ~ is_pos + (1 | donor).

    A sensitivity alternative to pseudobulk that keeps every cell but adds a
    per-donor random intercept, so the fixed ``is_pos`` effect is estimated
    against between-donor (not between-cell) variance. Off-target cells are pooled
    as the reference; cells whose condition is neither ``pos_label`` nor in
    ``neg_labels`` are dropped.

    Parameters
    ----------
    expr_cells : (n_cells,) per-cell expression for a single gene (any scale; a
        normalized/log scale is typical).
    donor : (n_cells,) donor id (the grouping / random-effect variable).
    condition : (n_cells,) condition label.
    pos_label, neg_labels : on- vs off-target labels (default off-target = rest).

    Returns
    -------
    dict(effect, se, p_value, converged). On any fit failure or non-convergence
    every numeric field is NaN and ``converged`` is False.
    """
    import pandas as pd
    import statsmodels.formula.api as smf

    expr_cells = np.asarray(expr_cells, dtype=float)
    donor = np.asarray(donor)
    condition = np.asarray(condition)
    neg = _neg_labels(neg_labels, condition, pos_label)

    is_pos = condition == pos_label
    is_neg = np.isin(condition, neg)
    keep = is_pos | is_neg

    fail = dict(effect=np.nan, se=np.nan, p_value=np.nan, converged=False)
    if keep.sum() == 0:
        return fail

    data = pd.DataFrame({
        "expression": expr_cells[keep],
        "is_pos": is_pos[keep].astype(float),
        "donor": donor[keep].astype(str),
    })

    try:
        model = smf.mixedlm("expression ~ is_pos", data, groups=data["donor"])
        res = model.fit(reml=True, method="lbfgs")
        converged = bool(getattr(res, "converged", True))
        effect = float(res.params["is_pos"])
        se = float(res.bse["is_pos"])
        p_value = float(res.pvalues["is_pos"])
        if not converged or not np.isfinite([effect, se, p_value]).all():
            return fail
        return dict(effect=effect, se=se, p_value=p_value, converged=True)
    except Exception:
        return fail


def _resolve_gene_index(gene_index_or_name, gene_names):
    """Accept either an integer column index or a gene name; return the index."""
    if isinstance(gene_index_or_name, (int, np.integer)):
        return int(gene_index_or_name)
    names = list(gene_names)
    return names.index(gene_index_or_name)


def forest_data(counts, donor, condition, gene_index_or_name, gene_names,
                pos_label, neg_labels=None, min_cells=10, n_boot=2000):
    """Per-donor log-CPM values and a bootstrap CI of the group difference, one gene.

    Data prep for a forest plot: each donor is one point (its log-CPM pseudobulk
    for the chosen gene), split into pos vs pooled off-target. The mean difference
    (pos - neg) gets a percentile bootstrap 95% CI resampling donors within each
    group. The RNG is ``numpy.random.default_rng(0)`` so the CI is deterministic.

    Parameters
    ----------
    gene_index_or_name : int column index into ``counts`` or a name in ``gene_names``.
    n_boot : bootstrap resamples (default 2000).

    Returns
    -------
    dict with keys ``pos_values, neg_values, pos_donors, neg_donors, pos_mean,
    neg_mean, diff, ci_low, ci_high``. ``diff = pos_mean - neg_mean`` on log-CPM.
    CI bounds are NaN if either group is empty.
    """
    (genes, pos_log, neg_log, _pos_lin, _neg_lin,
     pos_donors, neg_donors) = _pseudobulk_matrices(
        counts, donor, condition, gene_names, pos_label, neg_labels, min_cells)

    j = _resolve_gene_index(gene_index_or_name, gene_names)
    # map the caller's gene index/name onto the surviving pseudobulk columns
    target_name = gene_names[j] if isinstance(gene_index_or_name, (int, np.integer)) \
        else gene_index_or_name
    col = list(genes).index(target_name)

    pos_values = pos_log[:, col].astype(float)
    neg_values = neg_log[:, col].astype(float)
    pos_mean = float(pos_values.mean()) if pos_values.size else np.nan
    neg_mean = float(neg_values.mean()) if neg_values.size else np.nan
    diff = pos_mean - neg_mean

    ci_low = ci_high = np.nan
    if pos_values.size and neg_values.size:
        rng = np.random.default_rng(0)
        boot = np.empty(n_boot)
        for b in range(n_boot):
            pb = rng.choice(pos_values, size=pos_values.size, replace=True)
            nb = rng.choice(neg_values, size=neg_values.size, replace=True)
            boot[b] = pb.mean() - nb.mean()
        ci_low, ci_high = np.percentile(boot, [2.5, 97.5])

    return dict(
        pos_values=pos_values,
        neg_values=neg_values,
        pos_donors=list(pos_donors),
        neg_donors=list(neg_donors),
        pos_mean=pos_mean,
        neg_mean=neg_mean,
        diff=float(diff),
        ci_low=float(ci_low),
        ci_high=float(ci_high),
    )
