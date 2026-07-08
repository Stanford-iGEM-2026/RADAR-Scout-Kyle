import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  ResponsiveContainer,
  Tooltip,
  Legend,
} from 'recharts'
import ChartCard from './ChartCard'
import { COLORS, fmt } from '../lib/utils'

// Palette for up to 4 compared genes: crimson + teal + softened variants,
// all within the brand family (no foreign accent colors).
const SERIES_COLORS = [COLORS.crimson, COLORS.teal, COLORS.crimsonSoft, COLORS.tealSoft]

const ROWS = [
  { key: 'RACS', label: 'RACS', racs: true },
  { key: 'Sep', label: 'Sep' },
  { key: 'Feas', label: 'Feas' },
  { key: 'Repro', label: 'Repro' },
  { key: 'OffMax', label: 'OffMax' },
  { key: 'k_op', label: 'k_op', d: 1 },
  { key: 'Youden_J', label: 'Youden J', d: 2 },
  { key: 'n_donors', label: 'Donors', d: 0 },
]

function CmpTooltip({ active, payload, label }) {
  if (!active || !payload?.length) return null
  return (
    <div className="rc-tooltip">
      <div className="tt-title" style={{ fontStyle: 'normal' }}>
        {label}
      </div>
      {payload.map((s) => (
        <div className="tt-row" key={s.dataKey}>
          <span className="tt-k" style={{ color: s.color, fontStyle: 'italic' }}>
            {s.name}
          </span>
          <span>{fmt(s.value, 3)}</span>
        </div>
      ))}
    </div>
  )
}

export default function ComparePanel({ genes, onRemove, onClear }) {
  if (genes.length === 0) {
    return (
      <div className="panel-body">
        <p className="muted" style={{ margin: 0, fontSize: 13 }}>
          Add up to 4 genes with the <strong>+ compare</strong> links in the ranked list, or by
          selecting a gene and pressing “Add to compare”.
        </p>
      </div>
    )
  }

  // Grouped bar: one category per component, one bar per gene.
  const comps = ['Sep', 'Feas', 'Repro', 'OffMax', 'RACS']
  const chartData = comps.map((c) => {
    const row = { component: c }
    genes.forEach((g) => {
      row[g.gene] = g[c]
    })
    return row
  })

  return (
    <div className="panel-body">
      <div className="compare-bar" style={{ marginBottom: 18 }}>
        <span className="lead">Comparing {genes.length} of 4</span>
        <div className="compare-pills">
          {genes.map((g) => (
            <span className="pill" key={g.gene}>
              {g.gene}
              <button onClick={() => onRemove(g.gene)} title="Remove">
                ×
              </button>
            </span>
          ))}
        </div>
        <span className="spacer" />
        <button className="link-btn" onClick={onClear}>
          Clear all
        </button>
      </div>

      <div className="compare-table-wrap" style={{ marginBottom: 22 }}>
        <table className="compare-table">
          <thead>
            <tr>
              <th>Metric</th>
              {genes.map((g) => (
                <th key={g.gene} style={{ fontStyle: 'italic' }}>
                  {g.gene}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {ROWS.map((r) => (
              <tr key={r.key}>
                <td>{r.label}</td>
                {genes.map((g) => (
                  <td key={g.gene} className={r.racs ? 'racs' : undefined}>
                    {fmt(g[r.key], r.d ?? 3)}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <ChartCard
        title="Component comparison"
        sub="grouped by RACS factor"
        filename="compare_components"
        note="Higher is better for Sep / Feas / Repro / RACS; lower is better for OffMax."
      >
        <ResponsiveContainer width="100%" height={280}>
          <BarChart data={chartData} margin={{ top: 12, right: 16, bottom: 4, left: -8 }} barGap={2}>
            <CartesianGrid vertical={false} stroke={COLORS.border} />
            <XAxis dataKey="component" tick={{ fontSize: 12 }} stroke={COLORS.border} />
            <YAxis domain={[0, 1]} tick={{ fontSize: 11 }} stroke={COLORS.border} />
            <Tooltip content={<CmpTooltip />} cursor={{ fill: 'rgba(0,0,0,0.03)' }} />
            <Legend
              wrapperStyle={{ fontSize: 12, fontStyle: 'italic' }}
              iconType="square"
              iconSize={10}
            />
            {genes.map((g, i) => (
              <Bar
                key={g.gene}
                dataKey={g.gene}
                fill={SERIES_COLORS[i % SERIES_COLORS.length]}
                radius={[2, 2, 0, 0]}
                isAnimationActive={false}
                maxBarSize={26}
              />
            ))}
          </BarChart>
        </ResponsiveContainer>
      </ChartCard>
    </div>
  )
}
