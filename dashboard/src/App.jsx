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

  // view state
  const [rankMetric, setRankMetric] = useState('RACS') // 'RACS' | 'DSS'
  const [query, setQuery] = useState('')
  const [sortKey, setSortKey] = useState('RACS')
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
    jobs.push(
      info.has_umap
        ? fetch(`${base}data/${info.key}_umap.json`)
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

  // is DSS available for this disease?
  const hasDSS = useMemo(() => hasCol(genes, 'DSS'), [genes])
  const hasFDR = useMemo(() => hasCol(genes, 'FDR'), [genes])
  const isKeloid = info?.disease === 'keloid'

  // if a disease has no DSS, force the RACS view
  useEffect(() => {
    if (rankMetric !== 'RACS' && !hasCol(genes, rankMetric)) setRankMetric('RACS')
  }, [genes, hasDSS, rankMetric])

  // keep the display sort valid when switching diseases (e.g. DSS column gone)
  useEffect(() => {
    if (sortKey === 'DSS' && !hasDSS) setSortKey('RACS')
  }, [hasDSS, sortKey])

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
    let out = ranked.filter(
      (g) =>
        g.RACS >= filters.RACS &&
        g.Feas >= filters.Feas &&
        g.Sep >= filters.Sep &&
        g.Repro >= filters.Repro &&
        g.OffMax <= filters.OffMax &&
        (g.detect_P ?? Infinity) >= filters.detect_P &&
        (g.log2FC ?? Infinity) >= filters.log2FC,
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
          — surfacing {RANK_METRICS[rankMetric].blurb}
          {rankMetric === 'DSS'
            ? ' (high fold-change × abundance).'
            : ' (specific, abundant, reproducible, low off-target).'}
        </span>
      </div>

      {loading ? (
        <div className="loading">Loading {info?.disease} · {info?.cohort}…</div>
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
                hasDSS={hasDSS}
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
                  <DetailPanel gene={selectedGene} allGenes={ranked} rank={selectedGene?._rank} />
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

            {info?.has_umap && umap && umap.length ? (
              <div className="panel viz-panel">
                <div className="panel-head">
                  <h2>Cell embedding</h2>
                  <span className="hint">{info.cell_type} · {info.cohort}</span>
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
          RADAR-Scout · {info?.disease ?? '—'} · {info?.cell_type ?? ''} · {info?.cohort ?? ''} —
          scores follow docs/RACS_framework.md. RACS = RADAR compatibility; DSS = disease
          specificity.
        </span>
      </footer>
    </div>
  )
}
