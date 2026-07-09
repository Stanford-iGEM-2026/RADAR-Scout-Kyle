# RADAR-Scout — 2-week roadmap

Mapped to the project spec. `[x]` done, `[~]` in progress, `[ ]` todo.

## Week 1 — framework + data + scoring

### Foundation (done)
- [x] Repo scaffold, env, Modal wiring
- [x] **RACS mathematical framework** — theory + abundance/specificity optimum
      (`docs/RACS_framework.md`) — *primary research objective*
- [x] Core scorer: Sep / Feas / Repro / OffMax, donor-aware (`radar_scout/scoring.py`)
- [x] RADAR activation (Hill) model + reachable-threshold band (`radar_scout/hill.py`)
- [x] Donor-aware pseudobulk aggregation (`radar_scout/pseudobulk.py`)
- [x] Specificity indices for ablation (tau, significance-score) (`radar_scout/specificity.py`)
- [x] Unit tests (ideal > housekeeping > low-abundance) — all green
- [x] Modal Census pipeline: cheap disease scan + `build_and_score` (`modal_app/census_pull.py`)

### Data (done)
- [x] Census scan — keloid IS in Census (1,270 cells, 4 donors, 1 dataset)
- [x] First vertical end-to-end: **keloid** skin fibroblast (P, 4 donors) vs **normal**
      (H, 214 donors) + **localized scleroderma & injury** (R) — 3,487 candidate genes
- [x] Name harmonization: disease→MONDO/DOID, cell type→CL (`radar_scout/ontology.py`)
- [x] Technical-gene filter (`radar_scout/genesets.py`) — sex/ribosomal/mito/ncRNA/IEG
- [ ] GEO / EBI SCEA loaders for datasets not in Census (single-cell only)
- [ ] QC / doublet removal / **batch correction** (scanpy + scVI on Modal) — see finding below

### Stats (done)
- [x] Donor-aware DE: pseudobulk Welch/MWU + BH-FDR, cell-level MixedLM `(1|donor)`,
      forest-plot data (`radar_scout/de.py`)
- [ ] Wire forest plots into the figure/dashboard layer

## Week 2 — validation + dashboard

### Validation (paper-facing, all on real data)
- [x] V1 recover markers — melanoma → PRAME/GPNMB/S100B/SERPINE2; keloid (via
      **subpopulation ID**) → POSTN/ASPN/CTHRC1 + collagens by disease-specificity
- [x] V2 regime/knee plot (abundance vs specificity) — `scripts/figures.py`, `figures/`
- [x] V3 cross-cohort rank stability — `scripts/cross_cohort.py` (CELLxGENE vs GEO)
- [ ] V4 ablations: RACS vs abundance-only / specificity-only / DE-only
- [x] V5 cross-cohort/disease consensus — keloid consensus reproducible across 2 cohorts

### Dashboard (done — React + recharts, palette #8e1918 / #1c7170)
- [x] Disease·cohort selector (multi-disease) → ranked table
- [x] RACS ↔ DSS ranking toggle (RADAR-target vs disease-specificity views)
- [x] Filters (RACS, Feas, Sep, Repro, OffMax, detect_P, log2FC); compare up to 4
- [x] Viz: RACS breakdown, activation/window, knee, **UMAP, volcano, dot, heatmap**,
      cross-cohort consensus, full metric table
- [x] Export CSV + per-chart SVG/PNG; Gene → Plasmid + Primer hand-off
- [ ] Forest plots + PAGA into the dashboard (data prep exists in `de.py`)

## KEY FINDING — batch confounding (RESOLVED; write this up for the judges)
Naive whole-cell-type scoring on a single small cohort is dominated by technical
confounders and by cell-STATE dilution. Two fixes make it robust and are the
scientific story:
1. **Pathogenic subpopulation identification** (Task 4): Leiden-cluster, find the
   disease-enriched state, score *it*. Keloid POSTN detection 30%→61%; the ranking
   becomes the mesenchymal program (POSTN/ASPN/CTHRC1/collagens) instead of artifacts.
2. **Cross-cohort consensus** (§6): the same targets must be high in independent
   cohorts. Keloid CELLxGENE + GEO GSE163973 agree on COL1A1/POSTN/COL3A1/ASPN/…
   despite low overall rank correlation — that's the point of consensus.
Remaining hardening: pan-tissue off-target reference; batch integration (scVI/harmony)
before the UMAP; a 3rd keloid cohort (GSE181318/GSE220300).

## Flagged risks / must-verify before freeze
- [ ] **Calibrate the Hill parameters** (`K, n, L, K_lo, K_hi`) from the published
      RADAR/RADARS dose-response. Current values in `hill.py` are documented
      placeholders; all activation-dependent numbers (Feas/OffMax/RACS magnitudes,
      not the specificity ranking) are provisional until then. **Top blocker.**
- [ ] **Pan-tissue off-target** expansion to suppress ubiquitous-gene false positives.
- [ ] **Verify references** in `docs/RACS_framework.md` (RADAR kinetics, significance score).
- [ ] Replace `significance_score` stand-in with the exact Lu et al. (2014) formula.

## Cost discipline (Modal)
Smoke-test everything; the Census scan is near-free (reads a precomputed summary
table). Confirm before any large `build_and_score` / scVI run.
