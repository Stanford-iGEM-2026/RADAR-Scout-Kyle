import { COMPONENTS, RANK_METRICS, fmt } from '../lib/utils'

// Candidate sort columns, in display order. Only those present in the loaded
// records are shown as chips (diseases carry different metric subsets).
const SORT_CANDIDATES = [
  { key: 'pooled_score', label: 'Pooled' },
  { key: 'consensus_pct', label: 'Consensus' },
  { key: 'DSS', label: 'DSS' },
  { key: 'spec_score', label: 'Sensor' },
  ...COMPONENTS.filter((c) => c.key !== 'RACS'),
  { key: 'RACS', label: 'RACS' },
]

// Left ranked list. Rows show rank, italic gene name, a slim bar for the active
// rank metric, and its numeric value. The bar/value follow whichever metric is
// active (crimson for pooled/RACS, teal for DSS). Sortable by any present
// column; clicking selects.
export default function GeneList({
  genes,
  maxMetric,
  rankMetric,
  selected,
  onSelect,
  sortKey,
  sortDir,
  onSort,
  compareIds,
  onToggleCompare,
  compareEnabled,
}) {
  const rm = RANK_METRICS[rankMetric]
  const probe = genes[0]
  const sortChips = probe
    ? SORT_CANDIDATES.filter((c) => probe[c.key] !== undefined && probe[c.key] !== null)
    : SORT_CANDIDATES

  return (
    <div className="panel">
      <div className="panel-head">
        <h2>Ranked targets</h2>
        <span className="hint">by {rm.label} · {genes.length} shown</span>
      </div>

      <div className="list-toolbar">
        <span className="sort-label">Sort</span>
        {sortChips.map((c) => (
          <button
            key={c.key}
            className={`chip${sortKey === c.key ? ' active' : ''}`}
            onClick={() => onSort(c.key)}
          >
            {c.label}
            {sortKey === c.key ? <span className="dir">{sortDir === 'desc' ? '↓' : '↑'}</span> : null}
          </button>
        ))}
      </div>

      {genes.length === 0 ? (
        <div className="empty-note">No genes match the current filters or search.</div>
      ) : (
        <ul className="gene-list">
          {genes.map((g) => {
            const isSel = selected === g.gene
            const inCompare = compareIds.includes(g.gene)
            const val = g[rm.key]
            const width = maxMetric > 0 && Number.isFinite(val) ? Math.max(2, (val / maxMetric) * 100) : 0
            return (
              <li
                key={g.gene}
                className={`gene-row${isSel ? ' selected' : ''}`}
                onClick={() => onSelect(g.gene)}
              >
                <span className="rank">{g._rank}</span>
                <span className="gene-name">
                  {g.gene}
                  <span className="gene-meta">
                    n={g.n_donors}
                    {compareEnabled ? (
                      <>
                        {'  '}
                        <button
                          className="link-btn"
                          style={{ fontSize: 11 }}
                          onClick={(e) => {
                            e.stopPropagation()
                            onToggleCompare(g.gene)
                          }}
                        >
                          {inCompare ? '− compare' : '+ compare'}
                        </button>
                      </>
                    ) : null}
                  </span>
                </span>
                <span className="racs-cell">
                  <span className="racs-bar-track">
                    <span
                      className="racs-bar-fill"
                      style={{ width: `${width}%`, background: rm.color }}
                    />
                  </span>
                  <span className="racs-value">{fmt(val, rm.d)}</span>
                </span>
              </li>
            )
          })}
        </ul>
      )}
    </div>
  )
}
