// Slim top header: title, disease selector (one option per disease, pooled),
// gene search, a stats strip, and a compact list of the pooled cohorts so users
// can see what went into the consensus ranking.
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
  const cohorts = info?.cohorts || []

  // One option per disease: the disease name, plus a cohort count when >1.
  const optLabel = (d) =>
    d.n_cohorts > 1 ? `${d.disease} · ${d.n_cohorts} cohorts` : d.disease

  return (
    <header className="header">
      <div className="header-top">
        <div className="brand">
          <h1>
            RADAR<span className="accent">·</span>Scout
          </h1>
          <span className="tagline">pooled RNA target consensus</span>
        </div>
        <div className="controls">
          <div className="control">
            <label htmlFor="sel-disease">Disease</label>
            <select
              id="sel-disease"
              value={activeKey}
              onChange={(e) => setActiveKey(e.target.value)}
              style={{ minWidth: 220, textTransform: 'capitalize' }}
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
          <span className="value" style={{ textTransform: 'capitalize' }}>
            {info?.disease ?? '—'}
          </span>
          <span className="label">Disease</span>
        </div>
        <div className="stat">
          <span className="value">{info?.n_cohorts ?? '—'}</span>
          <span className="label">Cohorts pooled</span>
        </div>
        <div className="stat">
          <span className="value">
            {loading ? '…' : nCandidates}
            {!loading && nCandidates !== nTotal ? <span className="sub"> / {nTotal}</span> : null}
          </span>
          <span className="label">Genes ranked</span>
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
      </div>

      {cohorts.length ? (
        <div className="cohort-strip">
          <span className="cohort-lead">Pooled from</span>
          {cohorts.map((c) => (
            <span className="cohort-pill" key={c.key} title={c.cell_type}>
              {c.cohort}
              <span className="cohort-ct">{c.cell_type}</span>
            </span>
          ))}
        </div>
      ) : null}
    </header>
  )
}
