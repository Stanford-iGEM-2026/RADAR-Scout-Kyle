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

## Results — disease-agnostic, validated across diseases

The platform runs on any (disease, cell type) pair. Demonstrated on **real**
CELLxGENE Census + GEO data:

**Melanoma** (P = malignant cell, 64 donors; B = tumor microenvironment) — top
RADAR targets are bona-fide melanoma markers/therapeutic antigens: **PRAME, GPNMB,
S100B, SERPINE2, PLP1, GPM6B** (log2FC 3–4.6). This is the disease-agnostic engine
working end-to-end with a strong cell-type-specificity axis.

**Keloid** — the hard case (small cohorts). Three things unlock it:
1. **Pathogenic subpopulation identification** (spec Task 4): score the keloid-enriched
   *mesenchymal* fibroblast state, not all fibroblasts — POSTN detection 30% → 61%.
2. **Detection-difference specificity** (`spec_score`) — the sensor metric that surfaces
   near-binary markers (**ADAM12**: 92% of keloid MFB vs 3% of normal scar).
3. **Pooled cross-cohort ranking** across **3 cohorts** (CELLxGENE + GEO **GSE163973**
   + Deng-annotated MFB): **POSTN #1, ASPN #2, ADAM12 #6**, each ≈100th percentile in
   *every* cohort. This **reproduces the team's independent wet-lab sensor nomination
   (ADAM12 + POSTN)** end-to-end through the platform.

The dashboard has **one tab per disease** with a **Pooled** consensus ranking (default),
plus RACS / DSS / **Sensor** toggles and a per-gene "reproducibility across cohorts" view.

> Activation *magnitudes* remain provisional pending Hill calibration from the
> RADAR dose-response; the specificity/DSS *rankings* are calibration-independent.

Regenerate everything:
```bash
# score any disease on Modal:
modal run modal_app/census_pull.py::build_and_score --disease melanoma \
    --pathogenic-cell-types "malignant cell" --subcluster
# ingest an independent GEO cohort:
modal run modal_app/census_pull.py::ingest_and_score_geo --gse GSE163973
# figures, cross-cohort, dashboard data:
python scripts/figures.py outputs/<name>_racs.parquet --umap outputs/<name>_umap.parquet
python scripts/cross_cohort.py CELLxGENE=... GEO=... --rank-by DSS
python scripts/build_dashboard_data.py
cd dashboard && npm install && npm run dev   # multi-disease dashboard
```

## Design system

Brand palette: `#8e1918` (crimson), `#1c7170` (teal), white, black.

---
*Status: v0.3 — disease-agnostic pipeline (dedicated pathogenic fetch → subpopulation
ID), full metric set + RACS/DSS/spec_score/pooled rankings, donor-aware DE, ontology
harmonization, technical filter, doublet removal (scrublet) + batch correction (Harmony),
Modal Census + GEO ingestion, cross-cohort **pooled** ranking, the complete spec
visualization set (UMAP/volcano/scatter/bar/heatmap/violin/dot/forest/Manhattan/PAGA),
and a **one-tab-per-disease** React dashboard. Validated on melanoma, keloid (3 cohorts),
and pulmonary fibrosis — reproduces the wet-lab ADAM12+POSTN nomination. Plasmid step
intentionally out of scope. Remaining: Hill calibration (activation magnitudes only).
See [`ROADMAP.md`](ROADMAP.md).*
