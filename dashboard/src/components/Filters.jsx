import { useState } from 'react'

// Collapsible threshold filters.
//  - Four MIN sliders on the 0–1 RACS components (RACS, Feas, Sep, Repro)
//  - one MAX slider (OffMax)
//  - detection-frequency MIN (detect_P, %) and log2FC MIN, scaled to the data.
// Reports how many filters are non-default.
const MIN_FIELDS = [
  { key: 'RACS', label: 'RACS' },
  { key: 'Feas', label: 'Feas' },
  { key: 'Sep', label: 'Sep' },
  { key: 'Repro', label: 'Repro' },
]

export default function Filters({ filters, setFilters, activeCount, defaults, bounds }) {
  const [open, setOpen] = useState(false)

  const update = (key, value) => setFilters((f) => ({ ...f, [key]: value }))
  const reset = () => setFilters({ ...defaults })

  const dp = bounds?.detect_P || [0, 100]
  const fc = bounds?.log2FC || [0, 6]

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

          <div className="filters-divider" />

          <div className="slider-row">
            <label htmlFor="flt-detect">detect_P ≥</label>
            <input
              id="flt-detect"
              type="range"
              min={Math.floor(dp[0])}
              max={Math.ceil(dp[1])}
              step="1"
              value={filters.detect_P}
              onChange={(e) => update('detect_P', Number(e.target.value))}
            />
            <span className="thresh">{filters.detect_P.toFixed(0)}%</span>
          </div>
          <div className="slider-row">
            <label htmlFor="flt-log2fc">log2FC ≥</label>
            <input
              id="flt-log2fc"
              type="range"
              min={Math.floor(fc[0])}
              max={Math.ceil(fc[1])}
              step="0.1"
              value={filters.log2FC}
              onChange={(e) => update('log2FC', Number(e.target.value))}
            />
            <span className="thresh">{filters.log2FC.toFixed(1)}</span>
          </div>

          <div className="filters-actions">
            <span className="muted" style={{ fontSize: 12 }}>
              Off-target ceiling + detection floor keep housekeeping / sparse genes out.
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
