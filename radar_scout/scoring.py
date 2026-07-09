"""Donor-aware RACS scoring.

Implements the RADAR Activation Compatibility Score of docs/RACS_framework.md.
Everything here is array-based (per gene: expression, donor id, population label)
so it is testable without single-cell dependencies. The anti-pseudoreplication
principle is enforced structurally: every quantity is computed *within a donor*
and then averaged over donors; cells are never pooled across donors as if
independent (the only pooling is for picking a single operating threshold).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import rankdata, mannwhitneyu

from .hill import HillParams, DEFAULT_HILL, hill_activation, reachable_band


# --------------------------------------------------------------------------- #
# low-level helpers
# --------------------------------------------------------------------------- #
def _auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUC = P(pos > neg), tie-corrected. NaN if either side empty."""
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(np.concatenate([pos, neg]))
    r1 = r[:n1].sum()
    return float((r1 - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def _neg_labels(neg_labels, pop, pos_label):
    if neg_labels is not None:
        return list(neg_labels)
    return [p for p in np.unique(pop) if p != pos_label]


def _youden_threshold(pos_pool: np.ndarray, neg_pool: np.ndarray, band) -> tuple[float, float]:
    """Hard-threshold Youden-optimal operating point, clamped to the reachable band.

    O(n log n) via a single sort (the naive per-candidate broadcast is O(n^2) and
    blows up at ~1e4+ cells). Pooling across donors is acceptable *here* — it only
    selects one operating threshold; all scoring evaluations remain donor-aware.
    Returns (k_op, youden_J).
    """
    lo, hi = band
    n_pos, n_neg = pos_pool.size, neg_pool.size
    if n_pos == 0 or n_neg == 0:
        return float(lo), np.nan
    vals = np.concatenate([pos_pool, neg_pool])
    is_pos = np.concatenate([np.ones(n_pos, bool), np.zeros(n_neg, bool)])
    order = np.argsort(vals, kind="mergesort")[::-1]  # descending value
    vals_s = vals[order]
    pos_s = is_pos[order]
    tp = np.cumsum(pos_s)      # positives with value >= vals_s[k]  (TPR = tp/n_pos)
    fp = np.cumsum(~pos_s)     # negatives with value >= vals_s[k]  (FPR = fp/n_neg)
    j = tp / n_pos - fp / n_neg
    k = int(np.argmax(j))
    return float(np.clip(vals_s[k], lo, hi)), float(j[k])


def _per_donor_mean_activation(expr, donor, mask, params, K, min_cells):
    """dict donor -> mean RADAR activation over cells in ``mask`` for that donor.

    Donors with fewer than ``min_cells`` cells in the mask are omitted.
    """
    out: dict = {}
    if mask.sum() == 0:
        return out
    for d in np.unique(donor[mask]):
        cells = expr[(donor == d) & mask]
        if cells.size >= min_cells:
            out[d] = float(np.mean(hill_activation(cells, params, K=K)))
    return out


# --------------------------------------------------------------------------- #
# RACS component scores (donor-aware)
# --------------------------------------------------------------------------- #
def separability(expr, donor, pop, pos_label="P", neg_labels=None, min_cells=10) -> float:
    """Sep(g): DONOR-LEVEL AUC discriminating pathogenic vs off-target donors.

    Each donor is summarized by its mean expression (the donor is the experimental
    unit), then AUC = P(pathogenic-donor summary > off-target-donor summary). This
    is the correct anti-pseudoreplication specificity for cross-condition designs,
    where pathogenic and healthy cells come from *different individuals* (so a
    within-donor P-vs-O comparison is impossible). A donor contributing cells to
    both sides (mixed designs) contributes a summary to both sets.

    Threshold-free; captures cell-type + disease specificity + detection frequency.
    """
    expr = np.asarray(expr, float)
    donor = np.asarray(donor)
    pop = np.asarray(pop)
    neg = _neg_labels(neg_labels, pop, pos_label)
    is_pos = pop == pos_label
    is_neg = np.isin(pop, neg)
    pos_summ, neg_summ = [], []
    for d in np.unique(donor):
        dm = donor == d
        pc = expr[dm & is_pos]
        if pc.size >= min_cells:
            pos_summ.append(float(pc.mean()))
        nc = expr[dm & is_neg]
        if nc.size >= min_cells:
            neg_summ.append(float(nc.mean()))
    if not pos_summ or not neg_summ:
        return np.nan
    return _auc(np.asarray(pos_summ), np.asarray(neg_summ))


def feasibility(expr, donor, pop, pos_label="P", params: HillParams = DEFAULT_HILL,
                K=None, min_cells=10) -> float:
    """Feas(g): donor-mean on-target RADAR activation at the operating threshold.

    Encodes 'pathogenic abundance clears the activation threshold' (folds in
    detection frequency). ``K`` defaults to the nominal ``params.K``; the
    orchestrator passes the ROC-optimal clamped threshold.
    """
    expr = np.asarray(expr, float)
    donor = np.asarray(donor)
    pop = np.asarray(pop)
    per_donor = _per_donor_mean_activation(expr, donor, pop == pos_label, params, K, min_cells)
    return float(np.mean(list(per_donor.values()))) if per_donor else np.nan


def off_target_max(expr, donor, pop, pos_label="P", neg_labels=None,
                   params: HillParams = DEFAULT_HILL, K=None, min_cells=10):
    """OffMax(g): worst-case off-target activation across off-target populations.

    Returns (offmax, per_population_activation dict). Each population's activation
    is a donor-mean; OffMax is the max over populations.
    """
    expr = np.asarray(expr, float)
    donor = np.asarray(donor)
    pop = np.asarray(pop)
    neg = _neg_labels(neg_labels, pop, pos_label)
    per_pop = {}
    for o in neg:
        pd = _per_donor_mean_activation(expr, donor, pop == o, params, K, min_cells)
        if pd:
            per_pop[o] = float(np.mean(list(pd.values())))
    if not per_pop:
        return np.nan, {}
    return float(max(per_pop.values())), per_pop


def reproducibility(expr, donor, pop, pos_label="P",
                    params: HillParams = DEFAULT_HILL, K=None, min_cells=10) -> float:
    """Repro(g): consistency of on-target activation across PATHOGENIC donors.

    1 - CV over pathogenic donors of their mean RADAR activation, clipped to [0,1].
    A target that fires in every pathogenic donor scores high; one driven by a
    single donor scores low. NaN if < 2 pathogenic donors (absence of evidence,
    handled as neutral by ``racs``). This is the donor/cohort reproducibility term
    and the structural guard against a result resting on a single patient.
    """
    expr = np.asarray(expr, float)
    donor = np.asarray(donor)
    pop = np.asarray(pop)
    per_donor = _per_donor_mean_activation(expr, donor, pop == pos_label, params, K, min_cells)
    vals = np.asarray(list(per_donor.values()))
    if vals.size < 2:
        return np.nan
    cv = vals.std(ddof=1) / (vals.mean() + 1e-9)
    return float(np.clip(1.0 - cv, 0.0, 1.0))


def racs(sep, feas, repro, offmax, weights=(1.0, 1.0, 1.0, 1.0)) -> float:
    """Combine components into RACS in [0,1].

    Core factors (sep, feas, offmax) must be finite or RACS is NaN. A NaN repro
    (too few donors to assess) is treated as neutral (1.0) but should be surfaced
    to the user via ``GeneScore.n_donors``.
    """
    a, b, c, d = weights
    if any(v is None or np.isnan(v) for v in (sep, feas, offmax)):
        return np.nan
    r = 1.0 if (repro is None or np.isnan(repro)) else max(repro, 0.0)
    return float(max(sep, 0.0) ** a * max(feas, 0.0) ** b * r ** c * max(1.0 - offmax, 0.0) ** d)


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #
@dataclass
class GeneScore:
    gene: str
    racs: float
    sep: float
    feas: float
    repro: float
    offmax: float
    k_op: float
    youden_j: float
    n_donors: int
    per_pop_activation: dict = field(default_factory=dict)

    def as_row(self) -> dict:
        row = {
            "gene": self.gene, "RACS": self.racs, "Sep": self.sep, "Feas": self.feas,
            "Repro": self.repro, "OffMax": self.offmax, "k_op": self.k_op,
            "Youden_J": self.youden_j, "n_donors": self.n_donors,
        }
        for pop, val in self.per_pop_activation.items():
            row[f"act_{pop}"] = val
        return row


def score_gene(gene, expr, donor, pop, pos_label="P", neg_labels=None,
               params: HillParams = DEFAULT_HILL, weights=(1.0, 1.0, 1.0, 1.0),
               min_cells=10) -> GeneScore:
    """Compute all RACS components for one gene at the ROC-optimal operating point."""
    expr = np.asarray(expr, float)
    donor = np.asarray(donor)
    pop = np.asarray(pop)
    neg = _neg_labels(neg_labels, pop, pos_label)
    band = reachable_band(params)

    is_pos = pop == pos_label
    is_neg = np.isin(pop, neg)
    k_op, jstat = _youden_threshold(expr[is_pos], expr[is_neg], band)

    sep = separability(expr, donor, pop, pos_label, neg, min_cells)
    feas = feasibility(expr, donor, pop, pos_label, params, K=k_op, min_cells=min_cells)
    offmax, per_pop = off_target_max(expr, donor, pop, pos_label, neg, params, K=k_op, min_cells=min_cells)
    repro = reproducibility(expr, donor, pop, pos_label, params, K=k_op, min_cells=min_cells)
    score = racs(sep, feas, repro, offmax, weights)

    n_donors = int(sum(
        1 for d in np.unique(donor) if (expr[(donor == d) & is_pos]).size >= min_cells
    ))
    # include on-target activation in the per-population map for the dashboard
    pd_pos = _per_donor_mean_activation(expr, donor, is_pos, params, k_op, min_cells)
    if pd_pos:
        per_pop = {pos_label: float(np.mean(list(pd_pos.values()))), **per_pop}

    return GeneScore(
        gene=str(gene), racs=score, sep=sep, feas=feas, repro=repro, offmax=offmax,
        k_op=k_op, youden_j=jstat, n_donors=n_donors, per_pop_activation=per_pop,
    )


def _group_stats(x, groups):
    """Per-group (per-donor) mean expression and detection fraction for one gene."""
    means = np.array([x[idx].mean() for idx in groups])
    detect = np.array([(x[idx] > 0).mean() for idx in groups])
    return means, detect


def _bh_fdr(p):
    """Benjamini-Hochberg FDR, NaN-safe."""
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return out
    pv = p[ok]
    order = np.argsort(pv)
    n = pv.size
    q = pv[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    filled = np.empty(n)
    filled[order] = np.clip(q, 0, 1)
    out[ok] = filled
    return out


def score_matrix(expr2d, gene_names, donor, pop, pos_label="P", neg_labels=None,
                 params: HillParams = DEFAULT_HILL, weights=(1.0, 1.0, 1.0, 1.0),
                 min_cells=10):
    """Score every gene (columns of ``expr2d``) and return a ranked pandas DataFrame.

    Equivalent to looping ``score_gene`` but precomputes the (pop, donor) cell-index
    groups ONCE, so per-gene work touches only the relevant groups instead of
    rescanning all cells/donors. Scales to tens of thousands of genes x cells; run
    on Modal for the largest matrices (embarrassingly parallel over gene blocks).

    ``expr2d`` is (n_cells, n_genes) linear normalized abundance (CP10k).
    """
    import pandas as pd

    expr2d = np.asarray(expr2d, float)
    donor = np.asarray(donor)
    pop = np.asarray(pop)
    neg = _neg_labels(neg_labels, pop, pos_label)
    band = reachable_band(params)

    # precompute (pop, donor) -> cell index array, once for all genes
    gi = pd.DataFrame({"pop": pop, "donor": donor}).groupby(["pop", "donor"], sort=False).indices
    pos_groups = [idx for (p, _), idx in gi.items() if p == pos_label and idx.size >= min_cells]
    neg_groups_by_pop: dict = {}
    for (p, _), idx in gi.items():
        if p in neg and idx.size >= min_cells:
            neg_groups_by_pop.setdefault(p, []).append(idx)
    neg_groups_all = [idx for lst in neg_groups_by_pop.values() for idx in lst]

    pos_cells = np.concatenate(pos_groups) if pos_groups else np.array([], int)
    neg_cells = np.concatenate(neg_groups_all) if neg_groups_all else np.array([], int)
    n_pos_donors = len(pos_groups)

    eps = 1.0  # CP10k pseudocount for log2 fold-changes
    rows = []
    for j, g in enumerate(gene_names):
        x = expr2d[:, j]
        k_op, jstat = _youden_threshold(x[pos_cells], x[neg_cells], band)

        # --- pathogenic (P) donor-level abundance + detection + activation ---
        if pos_groups:
            pmeans, pdet = _group_stats(x, pos_groups)
            pos_act = np.array([hill_activation(x[idx], params, K=k_op).mean() for idx in pos_groups])
            feas = float(pos_act.mean())
            repro = (float(np.clip(1.0 - pos_act.std(ddof=1) / (pos_act.mean() + 1e-9), 0.0, 1.0))
                     if pos_act.size >= 2 else np.nan)
            mean_P, median_P = float(pmeans.mean()), float(np.median(pmeans))
            cv_P = float(pmeans.std(ddof=1) / (pmeans.mean() + 1e-9)) if pmeans.size >= 2 else np.nan
            dynrange = float(np.log2((pmeans.max() + eps) / (pmeans.min() + eps)))
            detect_P, pct_don_P = float(100 * pdet.mean()), float(100 * (pdet >= 0.1).mean())
        else:
            pmeans = np.array([])
            feas = repro = mean_P = median_P = cv_P = dynrange = detect_P = pct_don_P = np.nan

        # --- off-target populations (H healthy, B bystander cell types, R related disease) ---
        per_pop_stats, neg_all_means, offvals = {}, [], []
        for p, idxs in neg_groups_by_pop.items():
            m, d = _group_stats(x, idxs)
            a = np.array([hill_activation(x[idx], params, K=k_op).mean() for idx in idxs])
            per_pop_stats[p] = (m, d, a)
            neg_all_means.append(m)
            offvals.append(float(a.mean()))
        neg_all = np.concatenate(neg_all_means) if neg_all_means else np.array([])

        sep = _auc(pmeans, neg_all) if pmeans.size and neg_all.size else np.nan
        offmax = float(max(offvals)) if offvals else np.nan

        # donor-level significance (P vs healthy, else vs pooled off-target)
        ref = per_pop_stats["H"][0] if "H" in per_pop_stats else neg_all
        pval = np.nan
        if pmeans.size >= 2 and ref.size >= 2 and not np.allclose(np.r_[pmeans, ref], pmeans[0]):
            try:
                pval = float(mannwhitneyu(pmeans, ref, alternative="greater").pvalue)
            except Exception:
                pval = np.nan
        log2fc = float(np.log2((mean_P + eps) / (ref.mean() + eps))) if pmeans.size and ref.size else np.nan

        row = {"gene": str(g), "RACS": racs(sep, feas, repro, offmax, weights),
               "Sep": sep, "Feas": feas, "Repro": repro, "OffMax": offmax,
               "k_op": k_op, "Youden_J": jstat, "n_donors": n_pos_donors,
               "act_P": feas, "mean_P": mean_P, "median_P": median_P, "detect_P": detect_P,
               "cv_P": cv_P, "dynrange": dynrange, "pct_don_P": pct_don_P,
               "log2FC": log2fc, "p_value": pval,
               # Disease Specificity Score (reference-style: high transcription x high
               # fold-change) — surfaces disease-associated markers of the pathogenic pop.
               "DSS": (float(max(log2fc, 0.0) * np.log10(mean_P + 1.0))
                       if pmeans.size and np.isfinite(log2fc) and np.isfinite(mean_P) else np.nan)}
        for p, (m, d, a) in per_pop_stats.items():
            row[f"act_{p}"] = float(a.mean())
            row[f"mean_{p}"] = float(m.mean())
            row[f"detect_{p}"] = float(100 * d.mean())
            row[f"pct_don_{p}"] = float(100 * (d >= 0.1).mean())
            row[f"log2FC_{p}"] = float(np.log2((mean_P + eps) / (m.mean() + eps))) if pmeans.size else np.nan
        if "B" in per_pop_stats and pmeans.size:  # cell-type specificity (P vs other cell types)
            row["celltype_spec"] = _auc(pmeans, per_pop_stats["B"][0])
        if "R" in per_pop_stats and pmeans.size:  # disease specificity (P vs related diseases)
            row["disease_spec"] = _auc(pmeans, per_pop_stats["R"][0])
        # detection-difference specificity (logFC x delta detection). Surfaces near-binary
        # markers that abundance-weighted DSS misses (e.g. ADAM12: 92% in P vs 3% in ref).
        ref_det = None
        if "H" in per_pop_stats:
            ref_det = 100.0 * per_pop_stats["H"][1].mean()
        elif per_pop_stats:
            ref_det = float(np.mean([100.0 * d.mean() for (_, d, _) in per_pop_stats.values()]))
        if pmeans.size and ref_det is not None and np.isfinite(detect_P) and np.isfinite(log2fc):
            row["delta_detect"] = float(detect_P - ref_det)
            # up-regulated only (a sensor needs the target HIGH in pathogenic); the
            # positive-part product avoids scoring down-genes (neg x neg = pos).
            row["spec_score"] = float(max(log2fc, 0.0) * max(detect_P - ref_det, 0.0) / 100.0)
        rows.append(row)

    df = pd.DataFrame(rows)
    df["FDR"] = _bh_fdr(df["p_value"].to_numpy()) if "p_value" in df else np.nan
    return df.sort_values("RACS", ascending=False, na_position="last").reset_index(drop=True)
