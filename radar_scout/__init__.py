"""RADAR-Scout: single-cell RNA target prioritization for ADAR-based RNA sensors.

The public surface is intentionally small and array-based so that the scoring
logic is testable without heavy single-cell dependencies. The Modal pipeline
(``modal_app/``) extracts per-cell arrays from AnnData and calls into here.

See ``docs/RACS_framework.md`` for the theory these functions implement.
"""

from .hill import hill_activation, HillParams, DEFAULT_HILL, reachable_band
from .scoring import (
    separability,
    feasibility,
    off_target_max,
    reproducibility,
    racs,
    score_gene,
    GeneScore,
)
from .specificity import significance_score, tau_specificity
from .ontology import resolve, resolve_disease, resolve_cell_type, OntologyHit
from .genesets import is_technical, filter_technical
from .de import pseudobulk_de, mixedlm_de, forest_data
from .design import design_target, RadarDesign

__version__ = "0.1.0"

__all__ = [
    "hill_activation",
    "HillParams",
    "DEFAULT_HILL",
    "reachable_band",
    "separability",
    "feasibility",
    "off_target_max",
    "reproducibility",
    "racs",
    "score_gene",
    "GeneScore",
    "significance_score",
    "tau_specificity",
    "resolve",
    "resolve_disease",
    "resolve_cell_type",
    "OntologyHit",
    "is_technical",
    "filter_technical",
    "pseudobulk_de",
    "mixedlm_de",
    "forest_data",
    "design_target",
    "RadarDesign",
    "__version__",
]
