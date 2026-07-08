import { useState } from 'react'

// Collapsible threshold filters. Four MIN sliders (RACS, Feas, Sep, Repro) and
// one MAX slider (OffMax). Reports how many filters are non-default.
const MIN_FIELDS = [
  { key: 'RACS', label: 'RACS' },
  { key: 'Feas', label: 'Feas' },
  { key: 'Sep', label: 'Sep' },
  { key: 'Repro', label: 'Repro' },
]

export default function Filters({ filters, setFilters, activeCount }) {
  const [open, setOpen] = useState(false)

  const update = (key, value) => setFilters((f) => ({ ...f, [key]: value }))
  const reset = () =>
    setFilters({ RACS: 0, Feas: 0, Sep: 0, Repro: 0, OffMax: 1 })

  return (
    <div className="panel">
      <button className={`filters-toggle${open ? ' open' : ''}`} onClick={() => setOpen((o) => !o)}>
        <span className="caret">▶</span>
        Filters
        {activeCount > 0 ? <span className="count">{activeCount} active</span> : null}
      </button>

      {open ? (
        <div className="filters-body">
          {MIN_FIELDS.map((f) => (
            <div className="slider-row" key={f.key}>
              <label htmlFor={`flt-${f.key}`}>{f.label} ≥</label>
              <input
                id={`flt-${f.key}`}
                type="range"
                min="0"
                max="1"
                step="0.01"
                value={filters[f.key]}
                onChange={(e) => update(f.key, Number(e.target.value))}
              />
              <span className="thresh">{filters[f.key].toFixed(2)}</span>
            </div>
          ))}
          <div className="slider-row">
            <label htmlFor="flt-OffMax">OffMax ≤</label>
            <input
              id="flt-OffMax"
              type="range"
              min="0"
              max="1"
              step="0.01"
              value={filters.OffMax}
              onChange={(e) => update('OffMax', Number(e.target.value))}
            />
            <span className="thresh">{filters.OffMax.toFixed(2)}</span>
          </div>
          <div className="filters-actions">
            <span className="muted" style={{ fontSize: 12 }}>
              Off-target ceiling keeps housekeeping genes out.
            </span>
            <button className="link-btn" onClick={reset}>
              Reset
            </button>
          </div>
        </div>
      ) : null}
    </div>
  )
}
