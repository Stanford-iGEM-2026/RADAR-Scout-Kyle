import { useEffect, useMemo, useState } from 'react'
import Header from './components/Header'
import GeneList from './components/GeneList'
import Filters from './components/Filters'
import DetailPanel from './components/DetailPanel'
import ComparePanel from './components/ComparePanel'
import { exportCsv } from './lib/utils'

const DEFAULT_FILTERS = { RACS: 0, Feas: 0, Sep: 0, Repro: 0, OffMax: 1 }

export default function App() {
  const [genes, setGenes] = useState(null)
  const [meta, setMeta] = useState(null)
  const [error, setError] = useState(null)

  const [disease, setDisease] = useState('keloid')
  const [cellType, setCellType] = useState('skin fibroblast')
  const [query, setQuery] = useState('')
  const [sortKey, setSortKey] = useState('RACS')
  const [sortDir, setSortDir] = useState('desc')
  const [filters, setFilters] = useState(DEFAULT_FILTERS)
  const [selected, setSelected] = useState(null)
  const [compareOpen, setCompareOpen] = useState(false)
  const [compareIds, setCompareIds] = useState([])

  // ---- load sample data ---------------------------------------------------- //
  useEffect(() => {
    const base = import.meta.env.BASE_URL || '/'
    Promise.all([
      fetch(`${base}racs_results.json`).then((r) => {
        if (!r.ok) throw new Error(`racs_results.json ${r.status}`)
        return r.json()
      }),
      fetch(`${base}meta.json`)
        .then((r) => (r.ok ? r.json() : null))
        .catch(() => null),
    ])
      .then(([data, m]) => {
        setGenes(data)
        setMeta(m)
        if (data.length) setSelected(data[0].gene)
        if (m?.pathogenic_cell_types?.length) setCellType(m.pathogenic_cell_types[0])
      })
      .catch((e) => setError(e.message))
  }, [])

  // rank by RACS (stable, independent of the display sort)
  const ranked = useMemo(() => {
    if (!genes) return []
    const byRacs = [...genes].sort((a, b) => b.RACS - a.RACS)
    const rankMap = new Map(byRacs.map((g, i) => [g.gene, i + 1]))
    return genes.map((g) => ({ ...g, _rank: rankMap.get(g.gene) }))
  }, [genes])

  // ---- filter + search + sort --------------------------------------------- //
  const filtered = useMemo(() => {
    const q = query.trim().toLowerCase()
    let out = ranked.filter(
      (g) =>
        g.RACS >= filters.RACS &&
        g.Feas >= filters.Feas &&
        g.Sep >= filters.Sep &&
        g.Repro >= filters.Repro &&
        g.OffMax <= filters.OffMax,
    )
    if (q) out = out.filter((g) => g.gene.toLowerCase().includes(q))
    const dir = sortDir === 'desc' ? -1 : 1
    out = [...out].sort((a, b) => {
      const av = a[sortKey]
      const bv = b[sortKey]
      if (av === bv) return a._rank - b._rank
      return (av < bv ? -1 : 1) * dir
    })
    return out
  }, [ranked, filters, query, sortKey, sortDir])

  const maxRacs = useMemo(
    () => (ranked.length ? Math.max(...ranked.map((g) => g.RACS)) : 0),
    [ranked],
  )

  const activeFilterCount = useMemo(() => {
    let n = 0
    if (filters.RACS > 0) n++
    if (filters.Feas > 0) n++
    if (filters.Sep > 0) n++
    if (filters.Repro > 0) n++
    if (filters.OffMax < 1) n++
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

  if (error) {
    return (
      <div className="app">
        <div className="error">
          Failed to load data: {error}
          <div className="muted" style={{ marginTop: 8, fontSize: 13 }}>
            Ensure <code>public/racs_results.json</code> exists, then reload.
          </div>
        </div>
      </div>
    )
  }

  if (!genes) {
    return (
      <div className="app">
        <div className="loading">Loading RACS results…</div>
      </div>
    )
  }

  return (
    <div className="app">
      <Header
        meta={meta}
        disease={disease}
        setDisease={setDisease}
        cellType={cellType}
        setCellType={setCellType}
        query={query}
        setQuery={setQuery}
        nCandidates={filtered.length}
        nTotal={ranked.length}
      />

      <div className="main-grid">
        {/* left column: list + filters */}
        <div>
          <GeneList
            genes={filtered}
            maxRacs={maxRacs}
            selected={selected}
            onSelect={setSelected}
            sortKey={sortKey}
            sortDir={sortDir}
            onSort={onSort}
            compareIds={compareIds}
            onToggleCompare={toggleCompare}
            compareEnabled={compareOpen}
          />
          <Filters filters={filters} setFilters={setFilters} activeCount={activeFilterCount} />
          <div className="panel">
            <div className="panel-body" style={{ display: 'flex', gap: 10, flexWrap: 'wrap' }}>
              <button className="chip" onClick={() => exportCsv(filtered, 'racs_filtered.csv')}>
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
              />
            </>
          )}
        </div>
      </div>

      <footer style={{ marginTop: 40, paddingTop: 16, borderTop: '1px solid var(--border)' }}>
        <span className="muted" style={{ fontSize: 12 }}>
          RADAR-Scout · {meta?.disease ?? disease} · {cellType} — sample data (placeholder RACS
          scores; not calibrated). Scores follow docs/RACS_framework.md.
        </span>
      </footer>
    </div>
  )
}
