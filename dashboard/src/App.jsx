import { useEffect, useMemo, useState } from 'react'
import Header from './components/Header'
import GeneList from './components/GeneList'
import Filters from './components/Filters'
import DetailPanel from './components/DetailPanel'
import ComparePanel from './components/ComparePanel'
import UmapPanel from './components/UmapPanel'
import VolcanoPanel from './components/VolcanoPanel'
import CrossCohortPanel from './components/CrossCohortPanel'
import { RANK_METRICS, exportCsv, hasCol, extent } from './lib/utils'

const BASE_FILTERS = { RACS: 0, Feas: 0, Sep: 0, Repro: 0, OffMax: 1, detect_P: 0, log2FC: 0 }

// Non-gene keys inside a umap record; everything else is an expression column.
const UMAP_META_KEYS = new Set(['UMAP1', 'UMAP2', 'cell_type', 'pop', 'donor'])

export default function App() {
  // manifest + per-disease payloads
  const [diseases, setDiseases] = useState(null)
  const [activeKey, setActiveKey] = useState(null)
  const [genes, setGenes] = useState(null)
  const [umap, setUmap] = useState(null)
  const [crossCohort, setCrossCohort] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)

  // view state — pooled consensus ranking is the default lens
  const [rankMetric, setRankMetric] = useState('pooled_score')
  const [query, setQuery] = useState('')
  const [sortKey, setSortKey] = useState('pooled_score')
  const [sortDir, setSortDir] = useState('desc')
  const [filters, setFilters] = useState(BASE_FILTERS)
  const [selected, setSelected] = useState(null)
  const [compareOpen, setCompareOpen] = useState(false)
  const [compareIds, setCompareIds] = useState([])

  const base = import.meta.env.BASE_URL || '/'

  // ---- load the disease manifest + shared cross-cohort once ---------------- //
  useEffect(() => {
    Promise.all([
      fetch(`${base}data/diseases.json`).then((r) => {
        if (!r.ok) throw new Error(`diseases.json ${r.status}`)
        return r.json()
      }),
      fetch(`${base}data/cross_cohort.json`)
        .then((r) => (r.ok ? r.json() : null))
        .catch(() => null),
    ])
      .then(([manifest, xc]) => {
        const list = manifest?.diseases || []
        setDiseases(list)
        setCrossCohort(xc)
        if (list.length) setActiveKey(list[0].key)
      })
      .catch((e) => setError(e.message))
  }, [])

  const info = useMemo(
    () => (diseases && activeKey ? diseases.find((d) => d.key === activeKey) : null),
    [diseases, activeKey],
  )

  // ---- load the selected disease's gene records (+ umap if present) --------- //
  useEffect(() => {
    if (!info) return
    let cancelled = false
    setLoading(true)
    setGenes(null)
    setUmap(null)

    const jobs = [
      fetch(`${base}data/${info.key}_racs.json`).then((r) => {
        if (!r.ok) throw new Error(`${info.key}_racs.json ${r.status}`)
        return r.json()
      }),
    ]
    // The pooled key has no umap of its own — fetch the disease's umap_key file.
    const umapKey = info.has_umap && info.umap_key ? info.umap_key : null
    jobs.push(
      umapKey
        ? fetch(`${base}data/${umapKey}_umap.json`)
            .then((r) => (r.ok ? r.json() : null))
            .catch(() => null)
        : Promise.resolve(null),
    )

    Promise.all(jobs)
      .then(([racs, um]) => {
        if (cancelled) return
        setGenes(racs)
        setUmap(um)
        setSelected(racs.length ? racs[0].gene : null)
        setCompareIds([])
        setLoading(false)
      })
      .catch((e) => {
        if (!cancelled) {
          setError(e.message)
          setLoading(false)
        }
      })
    return () => {
      cancelled = true
    }
  }, [info])

  const hasFDR = useMemo(() => hasCol(genes, 'FDR'), [genes])
  const isKeloid = info?.disease === 'keloid'

  // if the active rank metric isn't available for this disease, fall back to
  // the pooled score (always present in the pooled files).
  useEffect(() => {
    if (genes && rankMetric !== 'pooled_score' && !hasCol(genes, rankMetric)) {
      setRankMetric('pooled_score')
    }
  }, [genes, rankMetric])

  // keep the display sort valid when switching diseases (a metric column may
  // vanish, e.g. DSS/Repro absent for some cohorts).
  useEffect(() => {
    if (genes && sortKey && !hasCol(genes, sortKey)) setSortKey(rankMetric)
  }, [genes, sortKey, rankMetric])

  // when the rank metric changes, default the display sort to match it
  useEffect(() => {
    setSortKey(rankMetric)
    setSortDir('desc')
  }, [rankMetric])

  // filter bounds for the dynamic sliders (detect_P, log2FC)
  const bounds = useMemo(() => {
    if (!genes) return null
    return {
      detect_P: extent(genes, 'detect_P'),
      log2FC: extent(genes, 'log2FC'),
    }
  }, [genes])

  // rank by the active metric (stable, independent of the display sort)
  const ranked = useMemo(() => {
    if (!genes) return []
    const mk = rankMetric
    const byMetric = [...genes].sort((a, b) => (b[mk] ?? -Infinity) - (a[mk] ?? -Infinity))
    const rankMap = new Map(byMetric.map((g, i) => [g.gene, i + 1]))
    return genes.map((g) => ({ ...g, _rank: rankMap.get(g.gene) }))
  }, [genes, rankMetric])

  // ---- filter + search + sort --------------------------------------------- //
  const filtered = useMemo(() => {
    const q = query.trim().toLowerCase()
    // Component metrics may be absent for some cohorts (e.g. Repro in PF); treat
    // a missing value as "passes" so the default (0) filters never hide genes.
    const geq = (v, t) => v === undefined || v === null || v >= t
    const leq = (v, t) => v === undefined || v === null || v <= t
    let out = ranked.filter(
      (g) =>
        geq(g.RACS, filters.RACS) &&
        geq(g.Feas, filters.Feas) &&
        geq(g.Sep, filters.Sep) &&
        geq(g.Repro, filters.Repro) &&
        leq(g.OffMax, filters.OffMax) &&
        geq(g.detect_P, filters.detect_P) &&
        geq(g.log2FC, filters.log2FC),
    )
    if (q) out = out.filter((g) => g.gene.toLowerCase().includes(q))
    const dir = sortDir === 'desc' ? -1 : 1
    out = [...out].sort((a, b) => {
      const av = a[sortKey]
      const bv = b[sortKey]
      if (av === bv) return a._rank - b._rank
      if (av === undefined || av === null) return 1
      if (bv === undefined || bv === null) return -1
      return (av < bv ? -1 : 1) * dir
    })
    return out
  }, [ranked, filters, query, sortKey, sortDir])

  const maxMetric = useMemo(() => {
    const mk = rankMetric
    return ranked.length ? Math.max(...ranked.map((g) => g[mk] ?? 0)) : 0
  }, [ranked, rankMetric])

  const activeFilterCount = useMemo(() => {
    let n = 0
    if (filters.RACS > 0) n++
    if (filters.Feas > 0) n++
    if (filters.Sep > 0) n++
    if (filters.Repro > 0) n++
    if (filters.OffMax < 1) n++
    if (filters.detect_P > 0) n++
    if (filters.log2FC > 0) n++
    return n
  }, [filters])

  const selectedGene = useMemo(
    () => ranked.find((g) => g.gene === selected) || null,
    [ranked, selected],
  )
  const compareGenes = useMemo(
    () => compareIds.map((id) => ranked.find((g) => g.gene === id)).filter(Boolean),
    [ranked, compareIds],
  )

  // umap expression columns (genes embedded in the umap json)
  const umapGeneCols = useMemo(() => {
    if (!umap || !umap.length) return []
    return Object.keys(umap[0]).filter((k) => !UMAP_META_KEYS.has(k))
  }, [umap])

  // the cohort whose single-cell data backs the umap (for the panel label)
  const umapCohort = useMemo(
    () => info?.cohorts?.find((c) => c.key === info.umap_key) || null,
    [info],
  )

  // keep selection valid when filters hide the current gene
  useEffect(() => {
    if (selected && !filtered.some((g) => g.gene === selected) && filtered.length) {
      setSelected(filtered[0].gene)
    }
  }, [filtered, selected])

  // ---- handlers ------------------------------------------------------------ //
  const onSort = (key) => {
    if (key === sortKey) {
      setSortDir((d) => (d === 'desc' ? 'asc' : 'desc'))
    } else {
      setSortKey(key)
      // OffMax reads best ascending (lower is better); others descending.
      setSortDir(key === 'OffMax' ? 'asc' : 'desc')
    }
  }

  const toggleCompare = (id) =>
    setCompareIds((ids) => {
      if (ids.includes(id)) return ids.filter((x) => x !== id)
      if (ids.length >= 4) return ids // cap at 4
      return [...ids, id]
    })

  const selectGene = (g) => {
    setSelected(g)
    setCompareOpen(false)
  }

  if (error) {
    return (
      <div className="app">
        <div className="error">
          Failed to load data: {error}
          <div className="muted" style={{ marginTop: 8, fontSize: 13 }}>
            Ensure <code>public/data/diseases.json</code> and the per-disease files exist, then
            reload.
          </div>
        </div>
      </div>
    )
  }

  if (!diseases) {
    return (
      <div className="app">
        <div className="loading">Loading disease manifest…</div>
      </div>
    )
  }

  return (
    <div className="app">
      <Header
        diseases={diseases}
        activeKey={activeKey}
        setActiveKey={setActiveKey}
        info={info}
        query={query}
        setQuery={setQuery}
        nCandidates={filtered.length}
        nTotal={ranked.length}
        loading={loading}
      />

      {/* rank-metric toggle */}
      <div className="rank-toggle-bar">
        <div className="rank-toggle">
          <span className="rt-label">Rank by</span>
          <div className="seg">
            {Object.values(RANK_METRICS).map((m) => (
              <button
                key={m.key}
                className={`seg-btn${rankMetric === m.key ? ' active' : ''}`}
                onClick={() => setRankMetric(m.key)}
                disabled={!hasCol(genes, m.key)}
                title={!hasCol(genes, m.key) ? `${m.label} not available for this cohort` : undefined}
              >
                {m.label}
              </button>
            ))}
          </div>
        </div>
        <span className="rt-hint">
          <strong style={{ color: RANK_METRICS[rankMetric].color }}>
            {RANK_METRICS[rankMetric].long}
          </strong>{' '}
          — surfacing {RANK_METRICS[rankMetric].blurb}.
        </span>
      </div>

      {loading ? (
        <div className="loading">
          Loading {info?.disease}
          {info?.n_cohorts ? ` · ${info.n_cohorts} cohort${info.n_cohorts > 1 ? 's' : ''}` : ''}…
        </div>
      ) : (
        <>
          <div className="main-grid">
            {/* left column: list + filters */}
            <div>
              <GeneList
                genes={filtered}
                maxMetric={maxMetric}
                rankMetric={rankMetric}
                selected={selected}
                onSelect={setSelected}
                sortKey={sortKey}
                sortDir={sortDir}
                onSort={onSort}
                compareIds={compareIds}
                onToggleCompare={toggleCompare}
                compareEnabled={compareOpen}
              />
              <Filters
                filters={filters}
                setFilters={setFilters}
                activeCount={activeFilterCount}
                defaults={BASE_FILTERS}
                bounds={bounds}
              />
              <div className="panel">
                <div className="panel-body" style={{ display: 'flex', gap: 10, flexWrap: 'wrap' }}>
                  <button className="chip" onClick={() => exportCsv(filtered, `${activeKey}_filtered.csv`)}>
                    Download CSV ({filtered.length})
                  </button>
                  <button
                    className={`chip${compareOpen ? ' active' : ''}`}
                    onClick={() => setCompareOpen((o) => !o)}
                  >
                    {compareOpen ? 'Exit compare mode' : 'Compare mode'}
                  </button>
                </div>
              </div>
            </div>

            {/* right column: detail or compare */}
            <div>
              {compareOpen ? (
                <div className="panel">
                  <div className="panel-head">
                    <h2>Compare targets</h2>
                    <span className="hint">up to 4 genes</span>
                  </div>
                  <ComparePanel
                    genes={compareGenes}
                    onRemove={(id) => toggleCompare(id)}
                    onClear={() => setCompareIds([])}
                  />
                </div>
              ) : (
                <>
                  {selectedGene ? (
                    <div style={{ marginBottom: 14, display: 'flex', gap: 10, flexWrap: 'wrap' }}>
                      <button
                        className="chip"
                        onClick={() => {
                          toggleCompare(selectedGene.gene)
                          setCompareOpen(true)
                        }}
                      >
                        + Add “{selectedGene.gene}” to compare
                      </button>
                    </div>
                  ) : null}
                  <DetailPanel
                    gene={selectedGene}
                    allGenes={ranked}
                    rank={selectedGene?._rank}
                    cohorts={info?.cohorts || []}
                  />
                </>
              )}
            </div>
          </div>

          {/* full-width analysis panels */}
          <div className="panels-stack">
            {hasFDR || hasCol(genes, 'p_value') ? (
              <div className="panel viz-panel">
                <div className="panel-head">
                  <h2>Differential expression</h2>
                  <span className="hint">volcano · click to select</span>
                </div>
                <div className="panel-body">
                  <VolcanoPanel
                    genes={genes}
                    selectedGene={selected}
                    onSelect={selectGene}
                    hasFDR={hasFDR}
                  />
                </div>
              </div>
            ) : null}

            {info?.has_umap && info?.umap_key && umap && umap.length ? (
              <div className="panel viz-panel">
                <div className="panel-head">
                  <h2>Cell embedding</h2>
                  <span className="hint">{umapCohort?.cell_type ?? info.disease} · {umapCohort?.cohort ?? 'single-cell atlas'}</span>
                </div>
                <div className="panel-body">
                  <UmapPanel points={umap} geneCols={umapGeneCols} selectedGene={selected} />
                </div>
              </div>
            ) : null}

            {isKeloid && crossCohort ? (
              <CrossCohortPanel rows={crossCohort} onSelect={selectGene} />
            ) : null}
          </div>
        </>
      )}

      <footer style={{ marginTop: 40, paddingTop: 16, borderTop: '1px solid var(--border)' }}>
        <span className="muted" style={{ fontSize: 12 }}>
          RADAR-Scout · {info?.disease ?? '—'}
          {info?.n_cohorts ? ` · pooled across ${info.n_cohorts} cohort${info.n_cohorts > 1 ? 's' : ''}` : ''} —
          scores follow docs/RACS_framework.md. Pooled = cross-cohort consensus; RACS = RADAR
          compatibility; DSS = disease specificity; Sensor = detection specificity.
        </span>
      </footer>
    </div>
  )
}
