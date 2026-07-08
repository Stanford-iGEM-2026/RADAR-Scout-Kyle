"""Correctness tests for the RACS scoring logic.

These use small toy arrays with a fixed seed purely to verify the math behaves
as designed (an ideal target beats a housekeeping gene and a too-lowly-expressed
specific gene). No biological conclusions are drawn from synthetic data — all
real target prioritization runs on imported single-cell datasets.
"""

import numpy as np
import pytest

from radar_scout import (
    hill_activation,
    DEFAULT_HILL,
    separability,
    racs,
    score_gene,
    tau_specificity,
)
from radar_scout.scoring import score_matrix
from radar_scout.pseudobulk import aggregate_pseudobulk, cpm_normalize

POPS = ["P", "H", "B", "R"]
N_DONORS = 3
N_PER = 80


def _draw(rng, n, frac, level, scale=0.4):
    """n cells: each 0 with prob (1-frac), else lognormal around `level` (CP10k)."""
    expressed = rng.random(n) < frac
    vals = rng.lognormal(mean=np.log(level), sigma=scale, size=n)
    return np.where(expressed, vals, 0.0)


def _make_dataset():
    rng = np.random.default_rng(0)
    donor, pop = [], []
    ideal, hk, lowspec = [], [], []
    specs = {
        "P": dict(ideal=(0.9, 25.0), hk=(0.95, 25.0), lowspec=(0.7, 0.6)),
        "H": dict(ideal=(0.1, 1.0), hk=(0.95, 25.0), lowspec=(0.1, 0.1)),
        "B": dict(ideal=(0.1, 1.0), hk=(0.95, 25.0), lowspec=(0.1, 0.1)),
        "R": dict(ideal=(0.1, 1.0), hk=(0.95, 25.0), lowspec=(0.1, 0.1)),
    }
    for d in range(N_DONORS):
        dfac = 1.0 + 0.1 * (d - 1)  # mild donor effect -> reproducibility < 1
        for p in POPS:
            f_i, l_i = specs[p]["ideal"]
            f_h, l_h = specs[p]["hk"]
            f_l, l_l = specs[p]["lowspec"]
            ideal.append(_draw(rng, N_PER, f_i, l_i * dfac))
            hk.append(_draw(rng, N_PER, f_h, l_h * dfac))
            lowspec.append(_draw(rng, N_PER, f_l, l_l * dfac))
            donor.append(np.full(N_PER, f"donor{d}"))
            pop.append(np.full(N_PER, p))
    donor = np.concatenate(donor)
    pop = np.concatenate(pop)
    expr = np.column_stack([np.concatenate(ideal), np.concatenate(hk), np.concatenate(lowspec)])
    return expr, ["IDEAL", "HK", "LOWSPEC"], donor, pop


# --------------------------------------------------------------------------- #
def test_hill_bounds_and_monotonic():
    x = np.linspace(0, 200, 50)
    a = hill_activation(x, DEFAULT_HILL)
    assert np.all(a >= DEFAULT_HILL.L - 1e-9)
    assert np.all(a < 1.0)
    assert np.all(np.diff(a) >= -1e-12)  # monotone non-decreasing
    assert hill_activation(0.0, DEFAULT_HILL) == pytest.approx(DEFAULT_HILL.L)


def test_racs_nan_propagation():
    assert np.isnan(racs(np.nan, 0.5, 0.5, 0.1))
    # NaN reproducibility is neutral, not fatal
    assert racs(0.9, 0.8, np.nan, 0.05) == pytest.approx(0.9 * 0.8 * (1 - 0.05))


def test_ideal_beats_housekeeping_and_lowabundance():
    expr, genes, donor, pop = _make_dataset()
    s = {g: score_gene(g, expr[:, j], donor, pop, pos_label="P")
         for j, g in enumerate(genes)}

    # Ideal target: separable, clears threshold, low off-target.
    assert s["IDEAL"].sep > 0.9
    assert s["IDEAL"].feas > 0.6
    assert s["IDEAL"].offmax < 0.2

    # Housekeeping: not separable (AUC ~ 0.5) despite high abundance.
    assert abs(s["HK"].sep - 0.5) < 0.15
    # Low-abundance specific: killed by feasibility (can't clear the floor).
    assert s["LOWSPEC"].feas < s["IDEAL"].feas

    # Ranking: the ideal RADAR target wins on the composite score.
    assert s["IDEAL"].racs > s["HK"].racs
    assert s["IDEAL"].racs > s["LOWSPEC"].racs
    for g in genes:
        assert 0.0 <= s[g].racs <= 1.0
        assert 0.0 <= s[g].repro <= 1.0 or np.isnan(s[g].repro)
        assert s[g].n_donors == N_DONORS


def test_score_matrix_ranks_ideal_first():
    expr, genes, donor, pop = _make_dataset()
    df = score_matrix(expr, genes, donor, pop, pos_label="P")
    assert list(df["gene"])[0] == "IDEAL"
    assert df["RACS"].is_monotonic_decreasing


def test_separability_direction():
    expr, genes, donor, pop = _make_dataset()
    sep_ideal = separability(expr[:, 0], donor, pop, "P")
    sep_hk = separability(expr[:, 1], donor, pop, "P")
    assert sep_ideal > sep_hk


def test_tau_extremes():
    assert tau_specificity([1, 0, 0, 0]) == pytest.approx(1.0)
    assert tau_specificity([5, 5, 5, 5]) == pytest.approx(0.0)


def test_pseudobulk_shapes():
    expr, genes, donor, pop = _make_dataset()
    pb = aggregate_pseudobulk(expr, donor, pop, genes)
    assert pb.shape[0] == N_DONORS * len(POPS)  # one row per (donor, pop)
    assert all(g in pb.columns for g in genes)
    norm = cpm_normalize(pb)
    assert norm[genes].to_numpy().min() >= 0.0


def test_nested_donor_design():
    """Real designs nest donors within condition (keloid patients vs healthy
    donors), so no donor has both P and O cells. Donor-level separability must
    still work — the old within-donor AUC would return NaN here and zero out RACS.
    """
    rng = np.random.default_rng(1)
    donor, pop, ideal, hk = [], [], [], []
    conditions = [("P", 4, (0.9, 25.0)), ("H", 3, (0.1, 1.0)), ("R", 3, (0.15, 1.5))]
    hk_spec = (0.95, 25.0)
    did = 0
    for cond, ndon, (frac, level) in conditions:
        for _ in range(ndon):
            ideal.append(_draw(rng, N_PER, frac, level))
            hk.append(_draw(rng, N_PER, *hk_spec))
            donor.append(np.full(N_PER, f"d{did}"))
            pop.append(np.full(N_PER, cond))
            did += 1
    donor = np.concatenate(donor)
    pop = np.concatenate(pop)
    expr = np.column_stack([np.concatenate(ideal), np.concatenate(hk)])

    si = score_gene("IDEAL", expr[:, 0], donor, pop, pos_label="P")
    sh = score_gene("HK", expr[:, 1], donor, pop, pos_label="P")

    assert not np.isnan(si.sep)      # the within-donor bug would make this NaN
    assert si.sep > 0.9              # keloid donors separate from healthy/related donors
    assert sh.sep < si.sep           # housekeeping is a weaker separator (noisy w/ 4 donors)
    assert si.n_donors == 4          # 4 pathogenic donors
    assert si.racs > 3 * sh.racs     # composite strongly favors the ideal target
