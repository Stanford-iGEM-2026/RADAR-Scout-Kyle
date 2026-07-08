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
from scipy.stats import rankdata

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

    rows = []
    for j, g in enumerate(gene_names):
        x = expr2d[:, j]
        k_op, jstat = _youden_threshold(x[pos_cells], x[neg_cells], band)

        pos_means = np.array([x[idx].mean() for idx in pos_groups])
        neg_means = np.array([x[idx].mean() for idx in neg_groups_all])
        sep = _auc(pos_means, neg_means) if pos_means.size and neg_means.size else np.nan

        if pos_groups:
            pos_act = np.array([hill_activation(x[idx], params, K=k_op).mean() for idx in pos_groups])
            feas = float(pos_act.mean())
            repro = (float(np.clip(1.0 - pos_act.std(ddof=1) / (pos_act.mean() + 1e-9), 0.0, 1.0))
                     if pos_act.size >= 2 else np.nan)
        else:
            feas = repro = np.nan

        per_pop = {pos_label: feas} if pos_groups else {}
        offvals = []
        for p, idxs in neg_groups_by_pop.items():
            act = float(np.mean([hill_activation(x[idx], params, K=k_op).mean() for idx in idxs]))
            per_pop[p] = act
            offvals.append(act)
        offmax = float(max(offvals)) if offvals else np.nan

        row = {"gene": str(g), "RACS": racs(sep, feas, repro, offmax, weights),
               "Sep": sep, "Feas": feas, "Repro": repro, "OffMax": offmax,
               "k_op": k_op, "Youden_J": jstat, "n_donors": n_pos_donors}
        row.update({f"act_{p}": v for p, v in per_pop.items()})
        rows.append(row)

    df = pd.DataFrame(rows)
    return df.sort_values("RACS", ascending=False, na_position="last").reset_index(drop=True)
