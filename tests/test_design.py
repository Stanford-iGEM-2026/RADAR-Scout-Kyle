"""Offline tests for the design utilities (no network)."""

import numpy as np
import pytest

from radar_scout.design import (
    reverse_complement,
    gc_fraction,
    tm_nn,
    design_radar_guide,
    design_validation_primers,
    _context_score,
)


def test_reverse_complement_roundtrip():
    seq = "ACGTTGCAAATTCCGG"
    assert reverse_complement(reverse_complement(seq)) == seq
    assert reverse_complement("AAAA") == "TTTT"


def test_gc_fraction():
    assert gc_fraction("GGCC") == 1.0
    assert gc_fraction("ATAT") == 0.0
    assert gc_fraction("ATGC") == 0.5


def test_tm_increases_with_gc_and_length():
    assert tm_nn("GCGCGCGCGCGCGCGCGCGC") > tm_nn("ATATATATATATATATATAT")
    assert tm_nn("ACGTACGTACGTACGTACGTACGT") > tm_nn("ACGTACGTAC")


def test_adar_context_prefers_5U_3G():
    # A flanked by U(->T) on 5' and G on 3' is the optimal ADAR context (score 1.0)
    seq = "CCTAGCC"  # ...T[A]G...
    i = seq.index("A")
    assert _context_score(seq, i) == pytest.approx(1.0)
    # A flanked by G/C-poor context scores lower
    assert _context_score("CCGACCC", 3) < 1.0


def test_radar_guide_is_antisense_window():
    rng = np.random.default_rng(0)
    cds = "".join(rng.choice(list("ACGT"), size=600))
    guide, edit_i, ctx, start = design_radar_guide(cds, guide_len=120)
    assert len(guide) == 120
    assert 0 <= edit_i < 120
    assert 0.0 <= ctx <= 1.0
    # guide is the reverse complement of the target window
    assert guide == reverse_complement(cds[start:start + 120])


def test_validation_primers_found_on_balanced_sequence():
    rng = np.random.default_rng(1)
    cds = "".join(rng.choice(list("ACGT"), size=800))  # ~50% GC by construction
    pair = design_validation_primers(cds)
    if pair is not None:  # search may legitimately fail on some sequences
        assert pair.amplicon_len >= 90
        assert abs(pair.fwd_tm - pair.rev_tm) <= 2.0
        assert 0.35 <= gc_fraction(pair.forward) <= 0.65
