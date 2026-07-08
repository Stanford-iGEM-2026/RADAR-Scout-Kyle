"""Correctness tests for donor-aware differential expression.

Toy NESTED-donor raw-count data with a fixed seed: each donor belongs to exactly
one condition (the real design — keloid patients vs healthy donors), so the DE
unit is the donor, not the cell. gene0 is truly up in the "P" condition; gene1 is
null (same distribution in both). No biology is inferred from synthetic data;
this only checks the statistics behave as designed.
"""

import numpy as np
import pytest

from radar_scout.de import pseudobulk_de, mixedlm_de, forest_data
from radar_scout.pseudobulk import aggregate_pseudobulk, cpm_normalize

N_POS_DONORS = 5   # condition "P"
N_NEG_DONORS = 6   # condition "N"
N_PER = 80         # cells per donor
GENES = ["gene0", "gene1", "gene2", "gene3"]


def _make_dataset():
    """Nested-donor raw counts (n_cells, 4 genes).

    gene0: strongly higher counts in P than N (true positive).
    gene1: identical distribution in P and N (true null).
    gene2, gene3: high-abundance stable genes, same in both conditions. They
        dominate the library so that gene0's up-regulation is only a small
        fraction of total depth — otherwise CPM normalization of one dominant DE
        gene compositionally deflates every other gene (a real single-cell
        pitfall we deliberately keep out of the toy null).
    A per-donor size factor scales all genes so donors differ in depth, exercising
    the CPM normalization and giving realistic between-donor variance.
    """
    rng = np.random.default_rng(0)
    donor, condition, rows = [], [], []
    did = 0
    for cond, ndon in [("P", N_POS_DONORS), ("N", N_NEG_DONORS)]:
        for _ in range(ndon):
            sf = rng.uniform(0.7, 1.4)  # donor depth factor (affects all genes)
            # Independent per-donor biological offsets on the signal genes give the
            # mixed model a real, non-singular donor random intercept to estimate
            # (without them the random-effect variance collapses to ~0 and the fit
            # sits on the boundary / fails to converge).
            g0_off = rng.lognormal(mean=0.0, sigma=0.25)
            g1_off = rng.lognormal(mean=0.0, sigma=0.25)
            g0_mean = (40.0 if cond == "P" else 8.0) * g0_off   # gene0 up in P
            g1_mean = 15.0 * g1_off                              # gene1 null
            g2_mean = 800.0                           # abundant stable gene
            g3_mean = 400.0                           # abundant stable gene
            g0 = rng.poisson(g0_mean * sf, size=N_PER)
            g1 = rng.poisson(g1_mean * sf, size=N_PER)
            g2 = rng.poisson(g2_mean * sf, size=N_PER)
            g3 = rng.poisson(g3_mean * sf, size=N_PER)
            rows.append(np.column_stack([g0, g1, g2, g3]))
            donor.append(np.full(N_PER, f"d{did}"))
            condition.append(np.full(N_PER, cond))
            did += 1
    counts = np.vstack(rows).astype(float)
    donor = np.concatenate(donor)
    condition = np.concatenate(condition)
    return counts, donor, condition


# --------------------------------------------------------------------------- #
def test_pseudobulk_de_columns_and_shape():
    counts, donor, condition = _make_dataset()
    df = pseudobulk_de(counts, donor, condition, GENES, pos_label="P")
    assert list(df.columns) == [
        "gene", "log2FC", "mean_pos", "mean_neg", "t_stat", "p_value",
        "p_mannwhitney", "FDR", "n_pos_donors", "n_neg_donors",
    ]
    assert len(df) == len(GENES)
    assert set(df["gene"]) == set(GENES)
    # donor is the experimental unit: 5 P donors, 6 N donors.
    assert (df["n_pos_donors"] == N_POS_DONORS).all()
    assert (df["n_neg_donors"] == N_NEG_DONORS).all()
    # sorted by FDR ascending
    fdr = df["FDR"].to_numpy()
    assert np.all(np.diff(fdr) >= -1e-12)


def test_pseudobulk_de_detects_true_positive_and_null():
    counts, donor, condition = _make_dataset()
    df = pseudobulk_de(counts, donor, condition, GENES, pos_label="P")
    row = df.set_index("gene")

    # gene0 is truly up in P: significant after FDR, positive fold change.
    assert row.loc["gene0", "FDR"] < 0.05
    assert row.loc["gene0", "log2FC"] > 0

    # gene1 is null: should not survive FDR.
    assert row.loc["gene1", "FDR"] > 0.1


def test_pseudobulk_de_default_neg_pools_all_others():
    """Omitting neg_labels pools every non-pos condition as the reference."""
    counts, donor, condition = _make_dataset()
    df_default = pseudobulk_de(counts, donor, condition, GENES, pos_label="P")
    df_explicit = pseudobulk_de(counts, donor, condition, GENES, pos_label="P",
                                neg_labels=["N"])
    a = df_default.set_index("gene").loc["gene0"]
    b = df_explicit.set_index("gene").loc["gene0"]
    assert a["n_neg_donors"] == b["n_neg_donors"] == N_NEG_DONORS
    assert a["p_value"] == pytest.approx(b["p_value"])


def _gene_expr_cells(counts, gene_col):
    """log1p-CP10k per-cell expression for one gene (a reasonable MixedLM scale)."""
    lib = counts.sum(axis=1)
    lib[lib == 0] = 1.0
    cp10k = counts[:, gene_col] / lib * 1e4
    return np.log1p(cp10k)


def test_mixedlm_de_true_positive():
    counts, donor, condition = _make_dataset()
    expr = _gene_expr_cells(counts, 0)  # gene0
    res = mixedlm_de(expr, donor, condition, pos_label="P")
    assert res["converged"] is True
    assert res["p_value"] < 0.05
    assert res["effect"] > 0


def test_mixedlm_de_null_graceful():
    counts, donor, condition = _make_dataset()
    expr = _gene_expr_cells(counts, 1)  # gene1 (null)
    res = mixedlm_de(expr, donor, condition, pos_label="P")
    # Either it converges and finds nothing, or it fails to converge — both are
    # handled gracefully (no exception, well-formed dict).
    assert set(res) == {"effect", "se", "p_value", "converged"}
    assert (res["converged"] is True and res["p_value"] > 0.1) or (res["converged"] is False)


def test_mixedlm_de_returns_nan_dict_on_empty():
    counts, donor, condition = _make_dataset()
    expr = _gene_expr_cells(counts, 0)
    # no cells match this pos_label -> graceful NaN dict, converged False
    res = mixedlm_de(expr, donor, condition, pos_label="DOES_NOT_EXIST",
                     neg_labels=["N"])
    assert res["converged"] is False
    assert np.isnan(res["effect"])


def test_forest_data_gene0():
    counts, donor, condition = _make_dataset()
    fd = forest_data(counts, donor, condition, "gene0", GENES, pos_label="P")
    assert len(fd["pos_values"]) == N_POS_DONORS
    assert len(fd["neg_values"]) == N_NEG_DONORS
    assert len(fd["pos_donors"]) == N_POS_DONORS
    assert len(fd["neg_donors"]) == N_NEG_DONORS
    assert fd["diff"] > 0                       # gene0 up in P on log-CPM too
    assert fd["ci_low"] <= fd["diff"] <= fd["ci_high"]


def test_forest_data_deterministic_and_by_index():
    """default_rng(0) bootstrap is reproducible; index and name select the same gene."""
    counts, donor, condition = _make_dataset()
    fd_name = forest_data(counts, donor, condition, "gene0", GENES, pos_label="P")
    fd_idx = forest_data(counts, donor, condition, 0, GENES, pos_label="P")
    assert fd_name["ci_low"] == pytest.approx(fd_idx["ci_low"])
    assert fd_name["ci_high"] == pytest.approx(fd_idx["ci_high"])
    # null gene CI should straddle zero
    fd_null = forest_data(counts, donor, condition, "gene1", GENES, pos_label="P")
    assert fd_null["ci_low"] < 0 < fd_null["ci_high"]
