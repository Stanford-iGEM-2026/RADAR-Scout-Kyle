"""Live tests for ontology harmonization (skip gracefully when offline).

These hit the EBI OLS API; they are skipped rather than failed when there is no
network, so the suite stays green in offline/CI environments.
"""

import pytest

from radar_scout.ontology import resolve_disease, resolve_cell_type


def _online(hit):
    if hit is None:
        pytest.skip("OLS unreachable (offline) — skipping live ontology test")
    return hit


def test_keloid_resolves_to_mondo():
    hit = _online(resolve_disease("keloid"))
    assert hit.id.startswith("MONDO:")
    assert hit.label.lower() == "keloid"


def test_disease_synonym_collapses():
    # 'mammary cancer' should collapse to the same canonical term as 'breast cancer'
    a = _online(resolve_disease("breast cancer"))
    b = _online(resolve_disease("mammary cancer"))
    assert a.id == b.id


def test_fibroblast_resolves_to_cl():
    hit = _online(resolve_cell_type("fibroblast"))
    assert hit.id.startswith("CL:")
    assert "fibroblast" in hit.label.lower()
