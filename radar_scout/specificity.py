"""Specificity indices used for annotation and ablations.

These are *not* the core RACS separability term (which is a donor-aware AUC in
``scoring.py``); they are provided so we can compare RACS against classical
specificity scores in validation V4, and to annotate candidates.
"""

from __future__ import annotations

import numpy as np


def tau_specificity(group_means) -> float:
    """Yanai's tissue-specificity index tau over per-group mean expression.

    tau = sum_i (1 - x_i / x_max) / (N - 1),  in [0, 1].
    tau = 0 means uniformly expressed; tau = 1 means expressed in a single group.

    Parameters
    ----------
    group_means : array-like of non-negative per-group mean expression.
    """
    x = np.asarray(group_means, dtype=float)
    x = np.clip(x, 0.0, None)
    n = x.size
    if n < 2:
        return np.nan
    xmax = x.max()
    if xmax <= 0:
        return np.nan
    return float(np.sum(1.0 - x / xmax) / (n - 1))


def significance_score(target_values, background_values) -> float:
    """Standardized enrichment of a target population vs background (Cohen's d).

    A specificity-style significance score: the standardized mean difference
    between the target population and the pooled background. Positive = enriched
    in the target. This is a transparent stand-in for the significance score of
    Lu et al. (Bioinformatics 2014); swap in their exact formulation before the
    paper freeze (tracked in ROADMAP.md).
    """
    t = np.asarray(target_values, dtype=float)
    b = np.asarray(background_values, dtype=float)
    if t.size < 2 or b.size < 2:
        return np.nan
    nt, nb = t.size, b.size
    sp2 = ((nt - 1) * t.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (nt + nb - 2)
    sp = np.sqrt(sp2)
    if sp == 0:
        return np.nan
    return float((t.mean() - b.mean()) / sp)
