# The RADAR Activation Compatibility Score (RACS)

**A signal-detection framework for prioritizing single-cell RNA targets for ADAR-based RNA sensors.**

Status: v0.1 (theory + reference implementation). This document is the primary
research deliverable of RADAR-Scout and is written to be self-contained,
reproducible, and citable. Symbols in `code font` map 1:1 to `radar_scout/`.

---

## 1. Problem statement

RADAR-class therapeutics (reprogrammable ADAR sensors — RADARS, CellREADR, and
relatives [1,2]) translate a payload **only when a target transcript is present
in the same cell**. A guide RNA hybridizes to the target mRNA; the resulting
double-stranded region recruits endogenous ADAR, which edits an in-frame stop
codon (UAG → UGG) and de-represses payload translation. Payload output is
therefore a **saturating, thresholded function of target abundance in the
individual cell** — not of the tissue average.

This makes target selection fundamentally different from conventional drug-target
discovery:

- The experimental unit that matters is the **single cell**, so target
  quantification must come from single-cell / single-nucleus RNA-seq. Bulk
  RNA-seq averages over cell composition and cannot tell us whether a transcript
  clears the activation threshold *inside a pathogenic cell*. Bulk is used here
  only as orthogonal validation, never for ranking.
- The design goal is **discrimination**, not maximization: activate in
  pathogenic cells, stay silent in healthy and off-target cells.
- The right statistical unit for *reproducibility* is the **biological donor**,
  not the cell. Thousands of cells from one patient are pseudoreplicates, not
  independent replicates [3,4]. All donor-level quantities below are computed
  per donor and then aggregated across donors.

The central scientific question — *what is the optimal relationship between
transcript abundance and specificity for RNA sensing?* — is answered in §4 as a
consequence of the activation model, not assumed.

---

## 2. The RADAR activation model

Let $x$ be the abundance of the target transcript in a given cell (normalized
counts, e.g. CP10k; see §6 for units and calibration). ADAR-sensor editing yield
as a function of target abundance is well described by a Hill function:

$$
a(x;\,K,n,L) \;=\; L \;+\; (1-L)\,\frac{x^{n}}{K^{n}+x^{n}}
$$

- $K$ — half-activation abundance (the **activation threshold**). The sensor
  designer can tune $K$ within a reachable band $[K_{\text{lo}}, K_{\text{hi}}]$
  by changing guide length, editing-site context, and ADAR recruitment.
- $n$ — effective steepness / cooperativity of the sensor response.
- $L$ — basal leak (payload with no target).

`a(x)` is **monotone increasing in $x$**. A cell "activates" with probability
$\approx a(x)$. This single quantity already folds together the two properties
the spec lists separately — *abundance* (through the magnitude of $x$) and
*detection frequency* (cells with $x=0$ contribute $a=L\approx 0$) — so we never
have to combine them with an ad-hoc weight.

Calibration of $(K,n,L)$ from the published RADAR/RADARS dose–response is
described in §6 and implemented in `radar_scout/hill.py`.

---

## 3. Donor-aware activation, and the ROC equivalence

### 3.1 Populations

For a chosen **disease** $d$ and **pathogenic cell type** $c$ (a Cell Ontology
term), define:

- $P$ — **on-target**: cells of type $c$ in disease $d$.
- Off-target set, partitioned so each maps to a spec requirement:
  - $H$ — healthy cells of type $c$ (healthy-tissue penalty),
  - $B$ — bystander cell types in the diseased tissue (cell-type specificity),
  - $R$ — cells of type $c$ in **biologically related diseases** (disease
    specificity; e.g. for keloid: normal scar, hypertrophic scar, other fibrosis).

### 3.2 Per-donor activation fractions

Let $S_j$ be the cells of population $S$ from donor $j$. The donor-level
activation fraction is

$$
\bar A_S^{(j)}(g;K,n) \;=\; \frac{1}{|S_j|}\sum_{i \in S_j} a\!\left(x_{gi};K,n\right).
$$

Each donor contributes **one** number per population; donors are then summarized
by their mean and between-donor standard error. Cells are never pooled across
donors as if independent — this is the anti-pseudoreplication core of the method.

### 3.3 Therapeutic window

The per-donor **therapeutic window** (worst-case over off-target populations) is

$$
W^{(j)}(g;K,n) \;=\; \bar A_P^{(j)} \;-\; \max_{o\in\{H,B,R\}} \bar A_o^{(j)} .
$$

### 3.4 Key result: threshold sweep = ROC

Because $a(\cdot;K,n)$ is monotone in $x$, and as $K$ sweeps $[0,\infty)$ it acts
as a soft threshold at $K$, the family of achievable operating points
$\big(\bar A_O(K),\,\bar A_P(K)\big)$ traces the **ROC curve** of the abundance
score $x_g$ for the label $P$ vs $O$. In the hard-threshold limit $n\to\infty$,
$a(x;K,\infty)=\mathbb{1}[x>K]$ and the correspondence is exact:

$$
\bar A_P(K) = \mathrm{TPR}(K), \qquad \bar A_O(K) = \mathrm{FPR}(K).
$$

Two consequences follow directly:

1. **Best achievable window = Youden's $J$** of the abundance classifier,
   attained at the ROC-optimal threshold $K^\star$ — *provided $K^\star$ is
   physically reachable*, $K^\star\in[K_\text{lo},K_\text{hi}]$:
   $$
   \max_K W(g;K) = \max_K\big[\mathrm{TPR}(K)-\mathrm{FPR}(K)\big] = J(g).
   $$
2. **Tuning-independent separability = AUC**:
   $$
   \mathrm{AUC}(g) = \Pr\big(x_{g}\!\mid\!P \;>\; x_{g}\!\mid\!O\big),
   $$
   the Mann–Whitney statistic — computed **per donor**, then averaged.

This is the theoretical bridge that justifies the whole platform: *RADAR target
selection is signal-detection over single-cell abundance distributions, gated by
ADAR Hill kinetics.* It is why neither the most-abundant nor the most-specific
transcript is optimal — the objective is **maximum separability at a physically
reachable operating point.**

---

## 4. The abundance–specificity optimum (primary research objective)

We now derive the requested optimal relationship. Model within-population
log-abundance as Gaussian with dropout; let $\mu_S$ be the mean log-abundance in
population $S$ and $\sigma$ a shared scale. Define the two competing properties:

- **Abundance** $\equiv \mu_P$ (how highly the target is expressed in pathogenic cells);
- **Specificity** $\equiv \Delta = \mu_P - \mu_O$, or the standardized separation
  $d = \Delta/\sigma$.

For a hard threshold at $K$ ($n\to\infty$), ignoring dropout,
$\mathrm{TPR}(K)=1-\Phi\!\big(\tfrac{\ln K-\mu_P}{\sigma}\big)$ and
$\mathrm{FPR}(K)=1-\Phi\!\big(\tfrac{\ln K-\mu_O}{\sigma}\big)$. The
Youden-optimal threshold (equal variance) is $\ln K^\star=(\mu_P+\mu_O)/2$,
giving

$$
J^\star \;=\; 2\,\Phi\!\left(\frac{d}{2}\right) - 1 .
$$

**In the specificity-limited regime the window depends only on the standardized
separation $d$ — absolute abundance $\mu_P$ drops out.** But RADAR reachability
constrains $\ln K \in [\ln K_\text{lo}, \ln K_\text{hi}]$. If the target is too
lowly expressed, the optimal threshold falls below the sensor floor
($\ln K^\star < \ln K_\text{lo}$) and we are forced to operate at $K_\text{lo}$:

$$
W(K_\text{lo}) = \Big[1-\Phi\big(\tfrac{\ln K_\text{lo}-\mu_P}{\sigma}\big)\Big]
              - \Big[1-\Phi\big(\tfrac{\ln K_\text{lo}-\mu_O}{\sigma}\big)\Big],
$$

which collapses as $\mu_P$ falls (pathogenic cells can no longer clear the
floor), and dropout makes it worse still. Hence **two regimes and a knee**:

| Regime | Condition | Window behaves like | What to optimize |
|---|---|---|---|
| Specificity-limited | $\mu_P$ high enough that $K^\star$ reachable | $2\Phi(d/2)-1$ | **separation $d$** |
| Threshold/abundance-limited | $\mu_P$ low, $K^\star<K_\text{lo}$ | fraction of $P$ above $K_\text{lo}$ | **abundance $\mu_P$** |

The RADAR-optimal target maximizes separation **subject to** the pathogenic
abundance clearing the reachable floor with margin. Writing the feasibility as a
sigmoid gate on abundance, the optimum is

$$
\boxed{\;g^\star=\arg\max_g\; \underbrace{\big[2\Phi(d_g/2)-1\big]}_{\text{specificity}}\;\cdot\;\underbrace{\Phi\!\left(\frac{\mu_{P,g}-\ln K_\text{lo}}{\sigma_g}-z_{\text{margin}}\right)}_{\text{abundance feasibility}}\;}
$$

This is a concrete, falsifiable prediction: on a plot of genes in
(abundance, specificity) space, the RADAR optimum is **not the top-right corner**
but the **knee** where the pathogenic lower tail just clears $K_\text{lo}$ while
separation is maximal. Overlaying the calibrated $K_\text{lo}$ on real data
(§7, validation V2) tests it. The boxed expression is exactly the
$\text{Sep}\cdot\text{Feas}$ core of RACS below — the score is the optimum, not a
heuristic bolted on afterward.

---

## 5. The RACS score

$$
\mathrm{RACS}(g)\;=\;\mathrm{Sep}(g)^{\alpha}\,\cdot\,\mathrm{Feas}(g)^{\beta}\,\cdot\,\mathrm{Repro}(g)^{\gamma}\,\cdot\,\big(1-\mathrm{OffMax}(g)\big)^{\delta}\;\in[0,1]
$$

with defaults $\alpha=\beta=\gamma=\delta=1$ (weights are configurable and can be
*fit* if labeled RADAR outcomes become available). Every factor is estimated
per donor and aggregated, and each maps to one spec requirement:

| Factor | Definition (donor-aware) | Spec requirement covered |
|---|---|---|
| $\mathrm{Sep}(g)$ | **donor-level** $\mathrm{AUC}$: pseudobulk each donor, then $P$-donors vs $O$-donors | cell-type + disease specificity, detection freq (threshold-free) |
| $\mathrm{Feas}(g)$ | donor-mean $\bar A_P$ at best reachable threshold $K^\star\wedge K_\text{lo}$ | abundance ≥ activation threshold, detection frequency |
| $\mathrm{Repro}(g)$ | $1-\mathrm{CV}$ of on-target activation across pathogenic donors | donor / cohort reproducibility |
| $\mathrm{OffMax}(g)$ | $\max_{o\in\{H,B,R\}}$ donor-mean $\bar A_o$ at $K^\star$ | healthy penalty + off-target + cross-disease |

`Sep` also reports the Youden $J$ and the ROC-optimal reachable threshold; a
plug-in alternative (the significance score of Lu et al. [5]) is provided in
`radar_scout/specificity.py` for ablations.

---

## 6. Estimation & calibration (implementation contract)

- **Normalization**: per-cell library-size normalization to CP10k, `log1p` for
  the Gaussian-regime analysis of §4; the Hill activation of §2 operates on
  **linear** CP10k. Genes measured on the same platform where possible; platform
  is carried as a covariate.
- **Pseudobulk**: sum raw counts per (donor, population), then normalize →
  one profile per donor per population. Used for mixed-effects differential
  expression (`statsmodels` MixedLM with a `(1|donor)` random effect, or a
  pseudobulk limma/DESeq2-style test) as an orthogonal DE readout.
- **Separability (AUC)** is estimated at the **donor level**: each donor is
  summarized by its mean expression (pseudobulk), then $\mathrm{AUC}$ compares
  pathogenic donors against off-target donors. This is the correct unit for
  cross-condition designs, where pathogenic and healthy cells come from
  *different individuals* — a within-donor $P$-vs-$O$ AUC is undefined. (When a
  donor does contribute cells to both sides, e.g. pathogenic vs bystander cells
  in the same patient, it contributes a summary to each set.)
- **Activation / feasibility / off-target** are per-donor means over that donor's
  cells, averaged across the donors of the relevant condition; SEs from
  between-donor variance. **Reproducibility** is the coefficient of variation of
  on-target activation across pathogenic donors.
- **Threshold calibration** $(K_\text{lo}, K_\text{hi}, n, L)$: digitized from the
  published RADAR/RADARS dose–response (sensor output vs. target abundance) and
  converted into CP10k units. Current values in `hill.py` are **documented
  placeholders** pending digitization — every downstream number that depends on
  them is flagged. This is the single most important empirical input to finalize.

---

## 7. Validation plan (reproducible, judge-facing)

All on **real** data (keloid is the first vertical — it is the running example in
the spec and connects to the team's antifibrotic circuit work):

- **V1 — Recover known biology.** High-RACS genes should include curated keloid
  fibroblast markers (cross-checked against Open Targets / literature);
  housekeeping genes should score low *despite* high abundance (killed by `Sep`);
  ultra-specific but very low genes should score low (killed by `Feas`).
- **V2 — Regime plot.** (abundance vs specificity) scatter colored by RACS with
  $K_\text{lo}$ overlaid → show the **knee** predicted in §4.
- **V3 — Donor holdout.** RACS ranking stable under donor subsampling and across
  independent cohorts (this is what `Repro` operationalizes).
- **V4 — Ablations.** RACS vs abundance-only vs specificity-only vs DE-only
  rankings; show RACS Pareto-dominates on (on-target activation, off-target leak).
- **V5 — Cross-disease.** keloid vs hypertrophic/normal scar → `OffMax` and the
  $R$ population drive disease specificity.

---

## 8. Assumptions & limitations

- Hill/Gaussian forms are modeling choices; §4's regime result is robust to them
  (it needs only monotone $a$ and stochastically-ordered $P>O$), but the exact
  optimum coordinates depend on calibration (§6).
- Normalized scRNA abundance is a proxy for the absolute intracellular target
  concentration the sensor actually sees; the calibration step is what ties the
  proxy to the activation axis and is the main source of systematic uncertainty.
- ADAR editing efficiency varies by sequence context and cell-intrinsic ADAR
  level; v0.1 treats these as absorbed into $(K,n,L)$. A per-cell ADAR covariate
  is a planned extension.

---

## References

1. Kaseniit KE, et al. *Modular, programmable RNA sensing using ADAR editing in
   living cells.* Nature Biotechnology (2023). (RADARS)
2. Qian Y, et al. *Programmable RNA sensing for cell monitoring and
   manipulation.* Nature (2022). (CellREADR)
3. Squair JW, et al. *Confronting false discoveries in single-cell differential
   expression.* Nature Communications (2021).
4. Zimmerman KD, et al. *A practical solution to pseudoreplication bias in
   single-cell studies.* Nature Communications (2021).
5. Lu Y, Yi Y, Liu P, et al. *A novel significance score for gene selection and
   ranking.* Bioinformatics (2014) 30(6):801. (specificity plug-in)
6. Crowell HL, et al. *muscat detects subpopulation-specific state transitions
   from multi-sample multi-condition single-cell transcriptomics data.* Nature
   Communications (2020). (pseudobulk DE)
7. Youden WJ. *Index for rating diagnostic tests.* Cancer (1950).

> Reference numbers 1, 2, 5 are the load-bearing external claims (RADAR kinetics,
> significance score). Verify exact venues/authors before the wiki/paper freeze;
> flagged in `ROADMAP.md`.
