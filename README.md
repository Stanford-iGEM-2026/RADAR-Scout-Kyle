# RADAR-Scout

**A single-cell RNA target-prioritization platform for RADAR-based, disease-agnostic therapeutics.**

RADAR-Scout ranks RNA transcripts by how reliably they will
activate a **RADAR** RNA sensor (ADAR-based) inside pathogenic cells while staying
silent in healthy and off-target cells. Inputs: a **disease** and a **cell type**.
Output: a ranked list of candidate targets with a **Gene of interest → Plasmid + Primer**.

Everything is built on **single-cell / single-nucleus RNA-seq only**.

---

## The core idea: RACS

The scientific heart of the project is the **RADAR Activation Compatibility Score
(RACS)**, a signal-detection framework that formalizes the trade-off between
transcript **abundance** and **specificity**. RADAR needs a target to exceed an
activation threshold, but the most abundant transcripts are rarely specific and
the most specific are often too scarce to fire the sensor. 

> **Full derivation:** [`docs/RACS_framework.md`](docs/RACS_framework.md).
> Result in one line: *RADAR target selection = maximizing the donor-aware AUC
> between pathogenic and off-target single-cell abundance at a physically
> reachable activation threshold.*

$$
\mathrm{RACS}(g)=\mathrm{Sep}(g)\cdot\mathrm{Feas}(g)\cdot\mathrm{Repro}(g)\cdot\big(1-\mathrm{OffMax}(g)\big)
$$

| Term | Meaning | Spec requirement |
|---|---|---|
| **Sep** | donor-level AUC, pathogenic vs off-target donors | cell-type + disease specificity, detection |
| **Feas** | on-target activation at the reachable threshold | abundance ≥ activation threshold |
| **Repro** | donor/cohort consistency of the window | reproducibility (no pseudoreplication) |
| **OffMax** | worst-case off-target activation | healthy penalty + off-target + cross-disease |

**Donors, not cells, are the experimental unit.** Every quantity is computed
within a donor and averaged across donors — the structural fix for
pseudoreplication bias.

---

## Repository layout

```
radar_scout/           # core Python package (array-based, unit-tested)
  hill.py              # RADAR/ADAR activation model + threshold band
  scoring.py           # RACS + donor-aware component scores (donor-level AUC)
  specificity.py       # tau / significance-score indices (annotation, ablations)
  pseudobulk.py        # donor-aware aggregation
  de.py                # donor-aware differential expression (pseudobulk + MixedLM)
  ontology.py          # free-text -> MONDO/DOID/CL harmonization (EBI OLS)
  genesets.py          # technical-gene filter (sex/ribosomal/mito/ncRNA/IEG)
  design.py            # Gene -> Plasmid + Primer hand-off (Ensembl + primer design)
docs/RACS_framework.md # the mathematical framework (primary deliverable)
modal_app/census_pull.py # Modal: Census scan/probe + build_and_score (heavy compute)
scripts/figures.py     # publication-quality figures from a RACS table
dashboard/             # React + recharts interactive platform (built; brand palette)
figures/               # generated example figures (keloid vertical)
tests/                 # correctness tests (39 assertions across scoring/de/ontology)
ROADMAP.md             # 2-week plan + the batch-confounding finding
```

## Quickstart

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pytest -q                                   # run the scoring tests

# heavy compute on Modal (single-cell pulls stay off the laptop):
modal run modal_app/census_pull.py          # cheap: scan what's in the Census
```

## Data sources

Single-cell (ranking): CELLxGENE Census, Human Cell Atlas, Broad Single Cell
Portal, GEO, EBI Single Cell Expression Atlas, Tabula Sapiens (healthy reference).
Annotation only (never quantification): Open Targets, Human Protein Atlas,
Reactome, Cell Ontology (CL) + MONDO/DOID for name harmonization.

## First results — keloid (provisional)

First end-to-end vertical on **real** CELLxGENE Census data: keloid skin
fibroblasts (P, 4 donors) vs normal skin fibroblasts (H, 214 donors) and related
fibrotic/wound conditions (R: localized scleroderma, injury), 3,487 candidate genes.

- After filtering technical genes, credible keloid biology surfaces —
  **HAS2, COL5A2, COL6A3, ZEB2, ITGB1** — alongside residual broadly-expressed
  genes, reflecting a real limitation: **keloid is a single cohort, so disease is
  confounded with batch.** See the "KEY FINDING" in [`ROADMAP.md`](ROADMAP.md).
- Scores are **provisional** pending Hill-parameter calibration from the RADAR
  dose-response; the *specificity ranking* is calibration-independent, the
  *activation magnitudes* are not.
- The dashboard loads these results out of the box; `scripts/figures.py` regenerates
  the figures in `figures/`.

Run the dashboard:
```bash
cd dashboard && npm install && npm run dev
```

## Design system

Brand palette: `#8e1918` (crimson), `#1c7170` (teal), white, black.

---
*Status: v0.1 — framework, donor-aware scorer, DE, ontology, technical filter, Modal
Census pipeline, design hand-off, figures, and a working React dashboard all landed
and tested on the keloid vertical. See [`ROADMAP.md`](ROADMAP.md).*
