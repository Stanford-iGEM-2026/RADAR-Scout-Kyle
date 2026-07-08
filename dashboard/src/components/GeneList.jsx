import { COMPONENTS, fmt } from '../lib/utils'

// Left ranked list. Rows show rank, italic gene name, a slim crimson RACS bar,
// and the numeric score. Sortable by RACS or any component. Clicking selects.
export default function GeneList({
  genes,
  maxRacs,
  selected,
  onSelect,
  sortKey,
  sortDir,
  onSort,
  compareIds,
  onToggleCompare,
  compareEnabled,
}) {
  return (
    <div className="panel">
      <div className="panel-head">
        <h2>Ranked targets</h2>
        <span className="hint">{genes.length} shown</span>
      </div>

      <div className="list-toolbar">
        <span className="sort-label">Sort</span>
        {COMPONENTS.map((c) => (
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
            const width = maxRacs > 0 ? Math.max(2, (g.RACS / maxRacs) * 100) : 0
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
                    <span className="racs-bar-fill" style={{ width: `${width}%` }} />
                  </span>
                  <span className="racs-value">{fmt(g.RACS, 3)}</span>
                </span>
              </li>
            )
          })}
        </ul>
      )}
    </div>
  )
}
