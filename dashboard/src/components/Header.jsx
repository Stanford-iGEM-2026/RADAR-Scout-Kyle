// Slim top header: title, disease/cohort selector (from diseases.json),
// gene search, and a stats strip driven by the selected disease's manifest entry.
export default function Header({
  diseases,
  activeKey,
  setActiveKey,
  info,
  query,
  setQuery,
  nCandidates,
  nTotal,
  loading,
}) {
  const donors = info?.n_donors || {}
  const subpop = info?.subpop && Object.keys(info.subpop).length ? info.subpop : null

  // Human label for each manifest entry in the switcher.
  const optLabel = (d) => `${d.disease} · ${d.cell_type} · ${d.cohort}`

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
            <label htmlFor="sel-disease">Disease · cohort</label>
            <select
              id="sel-disease"
              value={activeKey}
              onChange={(e) => setActiveKey(e.target.value)}
              style={{ minWidth: 280 }}
            >
              {diseases.map((d) => (
                <option key={d.key} value={d.key}>
                  {optLabel(d)}
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
            {loading ? '…' : nCandidates}
            {!loading && nCandidates !== nTotal ? <span className="sub"> / {nTotal}</span> : null}
          </span>
          <span className="label">Genes scored</span>
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
          <span className="value">{donors.B ?? '—'}</span>
          <span className="label">Donors · background</span>
        </div>
        <div className="stat">
          <span className="value" style={{ textTransform: 'capitalize' }}>
            {info?.cell_type ?? '—'}
          </span>
          <span className="label">Cell type</span>
        </div>
        <div className="stat">
          <span className="value" style={{ fontSize: 13, fontWeight: 500 }}>
            {info?.cohort ?? '—'}
          </span>
          <span className="label">Cohort</span>
        </div>
        {subpop ? (
          <div className="stat">
            <span className="value" style={{ fontSize: 13, fontWeight: 600, color: 'var(--crimson)' }}>
              cluster {subpop.path_cluster}
            </span>
            <span className="label">
              pathogenic subpopulation
              {subpop.n_P_after != null ? ` · ${subpop.n_P_after} cells` : ''}
            </span>
          </div>
        ) : null}
      </div>
    </header>
  )
}
