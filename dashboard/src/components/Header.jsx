// Slim top header: title, disease + cell-type selectors, gene search, stats strip.
export default function Header({
  meta,
  disease,
  setDisease,
  cellType,
  setCellType,
  query,
  setQuery,
  nCandidates,
  nTotal,
}) {
  const diseaseOptions = ['keloid']
  const cellOptions =
    meta?.pathogenic_cell_types?.length ? meta.pathogenic_cell_types : ['skin fibroblast']
  const donors = meta?.n_donors || {}

  return (
    <header className="header">
      <div className="header-top">
        <div className="brand">
          <h1>
            RADAR<span className="accent">·</span>Scout
          </h1>
          <span className="tagline">RNA target prioritization by RACS</span>
        </div>
        <div className="controls">
          <div className="control">
            <label htmlFor="sel-disease">Disease</label>
            <select id="sel-disease" value={disease} onChange={(e) => setDisease(e.target.value)}>
              {diseaseOptions.map((d) => (
                <option key={d} value={d}>
                  {d}
                </option>
              ))}
            </select>
          </div>
          <div className="control">
            <label htmlFor="sel-cell">Cell type</label>
            <select id="sel-cell" value={cellType} onChange={(e) => setCellType(e.target.value)}>
              {cellOptions.map((c) => (
                <option key={c} value={c}>
                  {c}
                </option>
              ))}
            </select>
          </div>
          <div className="control">
            <label htmlFor="search">Search gene</label>
            <input
              id="search"
              className="search"
              type="search"
              placeholder="e.g. POSTN"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              autoComplete="off"
            />
          </div>
        </div>
      </div>

      <div className="stats-strip">
        <div className="stat">
          <span className="value">
            {nCandidates}
            {nCandidates !== nTotal ? <span className="sub"> / {nTotal}</span> : null}
          </span>
          <span className="label">Candidate genes</span>
        </div>
        <div className="stat">
          <span className="value">{donors.P ?? '—'}</span>
          <span className="label">Donors · pathogenic</span>
        </div>
        <div className="stat">
          <span className="value">{donors.H ?? '—'}</span>
          <span className="label">Donors · healthy</span>
        </div>
        <div className="stat">
          <span className="value">{donors.R ?? '—'}</span>
          <span className="label">Donors · related</span>
        </div>
        <div className="stat">
          <span className="value" style={{ textTransform: 'capitalize' }}>
            {meta?.tissue ?? '—'}
          </span>
          <span className="label">Tissue</span>
        </div>
        <div className="stat">
          <span className="value" style={{ fontSize: 13, fontWeight: 500 }}>
            {meta?.census_version ?? '—'}
          </span>
          <span className="label">Census version</span>
        </div>
      </div>
    </header>
  )
}
