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
- [~] V1 recover keloid markers — **partial**: after filtering, real fibrosis genes
      surface (HAS2, COL5A2/6A3, ZEB2, ITGB1) mixed with batch-correlated ubiquitous
      genes; see finding below
- [x] V2 regime/knee plot (abundance vs specificity) — `scripts/figures.py`, `figures/`
- [ ] V3 donor-holdout ranking stability
- [ ] V4 ablations: RACS vs abundance-only / specificity-only / DE-only
- [ ] V5 cross-disease specificity (R is wired; formalize the comparison)

### Dashboard (done — React + recharts, palette #8e1918 / #1c7170)
- [x] Disease + cell-type input → ranked RACS table (loads real results)
- [x] Filters (RACS, Feas, Sep, Repro, OffMax)
- [x] Multi-transcript compare (up to 4)
- [x] Viz: RACS bar breakdown, activation/therapeutic-window, abundance-vs-specificity knee, off-target
- [x] Export CSV + per-chart SVG/PNG
- [x] Gene → Plasmid + Primer hand-off (`radar_scout/design.py`)
- [ ] Additional viz: UMAP, volcano, dot/heatmap, forest, PAGA (needs cell-level payloads to dashboard)

## KEY FINDING — batch confounding (write this up for the judges)
The unfiltered ranking was dominated by **technical confounders** (XIST → sex;
RPL*/RPS* → seq depth; JUN/NR4A1 → dissociation stress; NEAT1/MALAT1 → nuclear
lncRNA) because **keloid comes from a single cohort**, so disease is perfectly
confounded with batch — any technical signature gives donor-level AUC ≈ 1.0. The
gene filter removes these (and they are poor RADAR targets anyway). After filtering,
credible keloid biology appears (**HAS2, COL5A2, COL6A3, ZEB2, ITGB1**) alongside
residual ubiquitous/RNA-binding genes (PLCG2, FUS, splicing factors). Mitigations:
(a) require low off-target activation across the **related-disease** cohorts (already
in OffMax); (b) add a **pan-tissue** off-target reference to penalize broadly-expressed
genes; (c) obtain a second keloid cohort. This is an honest, defensible limitation of
single-cohort target discovery — and a strength to surface, not hide.

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
