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

### Data (next)
- [ ] **Run the Census scan** — confirm whether keloid is in CELLxGENE Census
      (else fall back to GEO loader or a related fibrosis for the first vertical)
- [ ] First vertical end-to-end: keloid fibroblast vs normal/hypertrophic scar + healthy skin
- [ ] Name harmonization: disease→MONDO/DOID, cell type→Cell Ontology (CL)
- [ ] GEO / EBI SCEA loaders for datasets not in Census (single-cell only)
- [ ] QC / doublet removal / batch correction (scanpy + scVI on Modal) — donor metadata preserved

### Stats
- [ ] Mixed-effects DE with `(1|donor)` (pseudobulk) as orthogonal DE readout
- [ ] Forest plots of donor-level effect sizes across cohorts

## Week 2 — validation + dashboard

### Validation (paper-facing, all on real data)
- [ ] V1 recover known keloid markers; housekeeping scores low; ultra-low-specific scores low
- [ ] V2 regime/knee plot (abundance vs specificity, K_lo overlay)
- [ ] V3 donor-holdout ranking stability
- [ ] V4 ablations: RACS vs abundance-only / specificity-only / DE-only
- [ ] V5 cross-disease specificity

### Dashboard (React / JS / HTML / CSS — palette #8e1918 / #1c7170)
- [ ] Disease + cell-type input → ranked RACS table
- [ ] Filters (abundance, threshold, detection, specificity, off-target)
- [ ] Multi-transcript compare
- [ ] Viz: UMAP, volcano, RACS bar breakdown, abundance-vs-specificity, dot/heatmap, forest, PAGA
- [ ] Export publication-quality figures + tables
- [ ] Gene → Plasmid + Primer design hand-off

## Flagged risks / must-verify before freeze
- [ ] **Calibrate the Hill parameters** (`K, n, L, K_lo, K_hi`) from the published
      RADAR/RADARS dose-response. Current values in `hill.py` are documented
      placeholders; every activation-dependent number depends on them.
- [ ] **Verify references** in `docs/RACS_framework.md` (RADAR kinetics paper,
      significance-score paper) — exact venue/authors before wiki/paper.
- [ ] **Keloid coverage in CELLxGENE Census** — confirm empirically via the scan.
- [ ] Replace `significance_score` stand-in with the exact Lu et al. (2014) formula.

## Cost discipline (Modal)
Smoke-test everything; the Census scan is near-free (reads a precomputed summary
table). Confirm before any large `build_and_score` / scVI run.
