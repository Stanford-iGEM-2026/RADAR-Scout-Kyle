import { useMemo } from 'react'
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ReferenceLine,
  LabelList,
} from 'recharts'
import ChartCard from './ChartCard'
import { COLORS, fmt } from '../lib/utils'

// Keloid cross-cohort robustness: consensus expression percentiles across the
// CELLxGENE atlas and the GEO GSE163973 cohort. Points near the diagonal are
// concordant; top consensus targets (POSTN, ASPN, collagens) sit top-right.
// A compact table lists the leading consensus genes.

function XCTooltip({ active, payload }) {
  if (!active || !payload?.length) return null
  const p = payload[0].payload
  return (
    <div className="rc-tooltip">
      <div className="tt-title">{p.gene}</div>
      <div className="tt-row">
        <span className="tt-k">CELLxGENE %ile</span>
        <span>{fmt(p.CELLxGENE, 1)}</span>
      </div>
      <div className="tt-row">
        <span className="tt-k">GEO %ile</span>
        <span>{fmt(p.GEO, 1)}</span>
      </div>
      <div className="tt-row">
        <span className="tt-k">consensus</span>
        <span>{fmt(p.consensus, 1)}</span>
      </div>
    </div>
  )
}

// Canonical keloid ECM markers we prefer to call out (they cluster top-right).
const PREFERRED = ['POSTN', 'ASPN', 'CTHRC1', 'COL1A1', 'FN1']

export default function CrossCohortPanel({ rows, onSelect }) {
  if (!rows || rows.length === 0) return null

  const sorted = useMemo(
    () => [...rows].sort((a, b) => b.consensus - a.consensus),
    [rows],
  )
  const top = sorted.slice(0, 12)

  // percentile axes are tight near the top; zoom the visible window a bit
  const lo = useMemo(() => {
    const mn = Math.min(...rows.map((r) => Math.min(r.CELLxGENE, r.GEO)))
    return Math.max(0, Math.floor(mn - 0.3))
  }, [rows])

  // Because the consensus targets pile into the top-right corner, give each
  // labelled gene a staggered vertical offset (and a right anchor) so the
  // callouts stack cleanly instead of overprinting.
  const labelOffsets = useMemo(() => {
    const present = PREFERRED.filter((g) => rows.some((r) => r.gene === g))
    const m = new Map()
    present.forEach((g, i) => m.set(g, -6 - i * 15))
    return m
  }, [rows])

  const renderLabel = (props) => {
    const { x, y, value } = props
    const dy = labelOffsets.get(value)
    if (dy === undefined) return null
    // These genes sit at the far right (≈100th pct), so extend the callout
    // leftward into the plot and draw a thin leader back to the point.
    return (
      <g>
        <line
          x1={x}
          y1={y}
          x2={x - 8}
          y2={y + dy + 3}
          stroke={COLORS.crimsonSoft}
          strokeWidth={0.75}
        />
        <text
          x={x - 10}
          y={y + dy}
          textAnchor="end"
          style={{ fontSize: 11, fontStyle: 'italic', fill: COLORS.crimson, fontWeight: 600 }}
        >
          {value}
        </text>
      </g>
    )
  }

  return (
    <div className="panel" style={{ marginTop: 24 }}>
      <div className="panel-head">
        <h2>Cross-cohort consensus</h2>
        <span className="hint">keloid · CELLxGENE vs GEO GSE163973</span>
      </div>
      <div className="panel-body xc-grid">
        <ChartCard
          title="Percentile concordance"
          sub="each point a gene · diagonal = perfect agreement"
          filename="cross_cohort_concordance"
          note="Genes hugging the top-right diagonal are robustly high in both cohorts — the consensus targets. Labelled genes are canonical keloid ECM markers."
        >
          <ResponsiveContainer width="100%" height={330}>
            <ScatterChart margin={{ top: 16, right: 20, bottom: 32, left: 8 }}>
              <CartesianGrid stroke={COLORS.border} />
              <ReferenceLine
                segment={[
                  { x: lo, y: lo },
                  { x: 100, y: 100 },
                ]}
                stroke={COLORS.borderStrong}
                strokeDasharray="4 4"
                ifOverflow="extendDomain"
              />
              <XAxis
                type="number"
                dataKey="CELLxGENE"
                domain={[lo, 100]}
                tick={{ fontSize: 11 }}
                stroke={COLORS.border}
                label={{
                  value: 'CELLxGENE percentile →',
                  position: 'insideBottom',
                  offset: -14,
                  fontSize: 11,
                  fill: COLORS.muted,
                }}
              />
              <YAxis
                type="number"
                dataKey="GEO"
                domain={[lo, 100]}
                tick={{ fontSize: 11 }}
                stroke={COLORS.border}
                label={{
                  value: 'GEO percentile →',
                  angle: -90,
                  position: 'insideLeft',
                  offset: 16,
                  fontSize: 11,
                  fill: COLORS.muted,
                }}
              />
              <Tooltip content={<XCTooltip />} cursor={{ strokeDasharray: '3 3' }} />
              <Scatter
                data={rows}
                fill={COLORS.teal}
                fillOpacity={0.7}
                isAnimationActive={false}
                onClick={(d) => d && onSelect?.(d.gene)}
                style={{ cursor: 'pointer' }}
              >
                <LabelList dataKey="gene" content={renderLabel} />
              </Scatter>
            </ScatterChart>
          </ResponsiveContainer>
        </ChartCard>

        <div className="xc-table-card">
          <div className="chart-head">
            <h3>
              Top consensus targets
              <span className="sub">  ranked by min(CELLxGENE, GEO)</span>
            </h3>
          </div>
          <div className="compare-table-wrap" style={{ padding: '4px 14px 12px' }}>
            <table className="compare-table">
              <thead>
                <tr>
                  <th>Gene</th>
                  <th>CxG</th>
                  <th>GEO</th>
                  <th>min</th>
                </tr>
              </thead>
              <tbody>
                {top.map((r) => (
                  <tr
                    key={r.gene}
                    onClick={() => onSelect?.(r.gene)}
                    style={{ cursor: onSelect ? 'pointer' : 'default' }}
                  >
                    <td>{r.gene}</td>
                    <td>{fmt(r.CELLxGENE, 1)}</td>
                    <td>{fmt(r.GEO, 1)}</td>
                    <td className="racs">{fmt(r.min_pct, 1)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>
  )
}
