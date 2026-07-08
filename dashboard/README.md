# RADAR-Scout Dashboard

Interactive front-end for **RADAR-Scout**: a clean, editorial view over the
**RADAR Activation Compatibility Score (RACS)** rankings. Pick a disease and cell
type, browse candidate target genes ranked by RACS, and inspect each gene's
signal-detection profile (separability, feasibility, reproducibility, off-target
leakage) and its per-population activation window.

> **Scores follow `docs/RACS_framework.md`:**
> `RACS = Sep · Feas · Repro · (1 − OffMax)`.
> The bundled data in `public/` is **placeholder** sample data (~25 genes,
> keloid skin fibroblasts) so the app renders standalone. Wire it to real Modal
> output by replacing `public/racs_results.json` and `public/meta.json`.

## Run

```bash
cd dashboard
npm install          # install deps (react, react-dom, recharts, vite)
npm run dev          # dev server → http://localhost:5173
```

Build and preview the production bundle:

```bash
npm run build        # outputs to dashboard/dist/
npm run preview      # serve the built bundle → http://localhost:4173
```

Requires Node 18+ (developed on Node 22, npm 10).

## Data contract

The app fetches two files from the site root (Vite serves `public/` at `/`):

**`/racs_results.json`** — array of gene objects, one per candidate transcript:

| field | type | meaning |
|---|---|---|
| `gene` | str | gene symbol |
| `RACS` | float [0,1] | composite score |
| `Sep` | float [0,1] | donor-level AUC, pathogenic vs off-target (specificity) |
| `Feas` | float [0,1] | on-target activation at the operating threshold (abundance) |
| `Repro` | float [0,1] | donor/cohort reproducibility |
| `OffMax` | float [0,1] | worst-case off-target activation |
| `k_op` | float | ROC-optimal operating threshold (CP10k) |
| `Youden_J` | float [0,1] | Youden J at `k_op` |
| `n_donors` | int | pathogenic donors backing the score |
| `act_P` / `act_H` / `act_R` | float [0,1] | activation in pathogenic / healthy / related populations |

These match `radar_scout.scoring.GeneScore.as_row()` exactly, so real Modal output
(`modal_app/census_pull.py`) drops in unchanged.

**`/meta.json`** — context for the header strip:

```json
{
  "disease": "keloid",
  "tissue": "skin",
  "pathogenic_cell_types": ["skin fibroblast", "myofibroblast"],
  "cell_counts": { "P": 18432, "H": 24107, "R": 9563 },
  "n_donors": { "P": 9, "H": 9, "R": 6 },
  "census_version": "2025-11-08"
}
```

## Features

- **Header** — disease + cell-type selectors, gene search, and a stats strip
  (candidate count, donors per population, tissue, Census version).
- **Ranked list** (left) — italic gene name, a slim crimson RACS bar, and the
  score. Sortable by RACS or any component (click a sort chip to toggle
  direction; OffMax defaults to ascending since lower is better).
- **Filters** (collapsible) — min sliders for RACS / Feas / Sep / Repro and a
  max slider for OffMax.
- **Detail panel** (right) — four charts for the selected gene:
  1. **RACS breakdown** — Sep, Feas, Repro, (1 − OffMax) in teal, RACS in crimson.
  2. **Per-population activation** — on-target (crimson) vs off-target (teal); the
     therapeutic-window view, labeled with `k_op`.
  3. **Abundance vs specificity** — scatter of all genes (x = Feas, y = Sep, size
     ∝ RACS) with the high-Sep/high-Feas **RADAR sweet spot** annotated; selected
     gene highlighted in crimson.
  4. **Off-target check** — on-target `act_P` vs `max(act_H, act_R)`.
- **Compare mode** — select up to 4 genes → a comparison table + a grouped bar of
  their component scores.
- **Export** — “Download CSV” of the filtered table; every chart has **SVG** and
  **PNG** download buttons (recharts renders SVG; PNG is rasterized via canvas).

## Design

Nature/Cell editorial aesthetic. Brand palette only — crimson `#8e1918`
(primary), teal `#1c7170` (secondary), white background, near-black `#141414`
text, `#e6e6e6` borders, `#6b6b6b` muted. No shadows, gradients, 3D, or
chartjunk; thin 1px borders and generous whitespace throughout.

## Component structure

```
dashboard/
  index.html                  # Vite entry
  vite.config.js              # @vitejs/plugin-react
  package.json
  public/
    racs_results.json         # PLACEHOLDER sample: ~25 genes
    meta.json                 # PLACEHOLDER sample: disease/tissue/donors
  src/
    main.jsx                  # React root
    styles.css                # single stylesheet (brand tokens + layout)
    App.jsx                   # state: data load, filter/sort/search, layout
    lib/
      utils.js                # COLORS, formatting, CSV + chart SVG/PNG export
    components/
      Header.jsx              # selectors, search, stats strip
      GeneList.jsx            # ranked list + sort chips
      Filters.jsx             # collapsible threshold sliders
      DetailPanel.jsx         # 4 charts for the selected gene
      ComparePanel.jsx        # compare table + grouped bar
      ChartCard.jsx           # titled chart wrapper + SVG/PNG download
```

## Wiring to real data

`modal_app/census_pull.py::build_and_score` writes a ranked table from the
donor-aware RACS scorer. Export it as JSON with the fields above (it already uses
`GeneScore.as_row()`), drop the file at `public/racs_results.json`, write a
matching `public/meta.json`, and rebuild. No code changes required.
