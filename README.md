# RADAR-Scout

**A single-cell RNA target-prioritization platform for RADAR-based, disease-agnostic therapeutics.**

Stanford iGEM 2026. RADAR-Scout ranks RNA transcripts by how reliably they will
activate a **RADAR** RNA sensor (ADAR-based) inside pathogenic cells while staying
silent in healthy and off-target cells. Inputs: a **disease** and a **cell type**.
Output: a ranked list of candidate targets with a **Gene of interest → Plasmid + Primer**
design hand-off.

Everything is built on **single-cell / single-nucleus RNA-seq only** — bulk RNA-seq
is never used for ranking (it averages over cell composition and cannot tell us
whether a transcript clears the sensor's activation threshold *inside a cell*).

---

## The core idea: RACS

The scientific heart of the project is the **RADAR Activation Compatibility Score
(RACS)** — a signal-detection framework that formalizes the trade-off between
transcript **abundance** and **specificity**. RADAR needs a target to exceed an
activation threshold, but the most abundant transcripts are rarely specific and
the most specific are often too scarce to fire the sensor. RACS derives the
optimal balance from the ADAR activation kinetics rather than assuming it.

> **Full derivation:** [`docs/RACS_framework.md`](docs/RACS_framework.md).
> Result in one line: *RADAR target selection = maximizing the donor-aware AUC
> between pathogenic and off-target single-cell abundance at a physically
> reachable activation threshold.*

$$
\mathrm{RACS}(g)=\mathrm{Sep}(g)\cdot\mathrm{Feas}(g)\cdot\mathrm{Repro}(g)\cdot\big(1-\mathrm{OffMax}(g)\big)
$$

| Term | Meaning | Spec requirement |
|---|---|---|
| **Sep** | donor-mean AUC, pathogenic vs off-target | cell-type + disease specificity, detection |
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
  scoring.py           # RACS + donor-aware component scores
  specificity.py       # tau / significance-score indices (annotation, ablations)
  pseudobulk.py        # donor-aware aggregation for mixed-effects DE
docs/
  RACS_framework.md    # the mathematical framework (primary deliverable)
modal_app/
  census_pull.py       # Modal: CELLxGENE Census pull + scoring (heavy compute)
tests/                 # correctness tests (toy arrays only)
dashboard/             # React/JS interactive platform (WIP)
ROADMAP.md             # 2-week plan mapped to the spec
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

## Design system

Brand palette: `#8e1918` (crimson), `#1c7170` (teal), white, black.

---
*Status: v0.1 — framework + scorer + Modal pipeline landed and tested. See
[`ROADMAP.md`](ROADMAP.md).*
