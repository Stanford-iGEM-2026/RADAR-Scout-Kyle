import { useMemo } from 'react'
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  ZAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ReferenceLine,
  LabelList,
} from 'recharts'
import ChartCard from './ChartCard'
import { COLORS, fmt } from '../lib/utils'

// Differential-expression volcano: x = log2FC (pathogenic vs reference),
// y = -log10(FDR) (falls back to p_value if FDR absent). Significant genes
// (FDR<0.05 & |log2FC|>1) are crimson; the rest teal-grey. Top genes labelled.
// Clicking any point selects that gene in the rest of the dashboard.

const FDR_CUT = 0.05
const FC_CUT = 1

function VolcanoTooltip({ active, payload, yField }) {
  if (!active || !payload?.length) return null
  const p = payload[0].payload
  return (
    <div className="rc-tooltip">
      <div className="tt-title">{p.gene}</div>
      <div className="tt-row">
        <span className="tt-k">log2FC</span>
        <span>{fmt(p.log2FC, 2)}</span>
      </div>
      <div className="tt-row">
        <span className="tt-k">{yField}</span>
        <span>{p.rawP < 1e-4 ? p.rawP.toExponential(1) : fmt(p.rawP, 4)}</span>
      </div>
      <div className="tt-row">
        <span className="tt-k">-log10</span>
        <span>{fmt(p.y, 2)}</span>
      </div>
    </div>
  )
}

export default function VolcanoPanel({ genes, selectedGene, onSelect, hasFDR }) {
  const yField = hasFDR ? 'FDR' : 'p_value'

  const { sig, nonsig, xDomain, labelSet } = useMemo(() => {
    const pts = genes
      .map((g) => {
        const raw = g[yField]
        if (raw === undefined || raw === null || g.log2FC === undefined || g.log2FC === null) return null
        // clamp p=0 to a floor so -log10 is finite
        const p = Math.max(raw, 1e-10)
        const y = -Math.log10(p)
        const significant = raw < FDR_CUT && Math.abs(g.log2FC) > FC_CUT
        return {
          gene: g.gene,
          log2FC: g.log2FC,
          y,
          rawP: raw,
          significant,
          z: g.gene === selectedGene ? 200 : 60,
        }
      })
      .filter(Boolean)

    const xs = pts.map((p) => p.log2FC)
    const xmax = Math.max(1.5, ...xs.map(Math.abs))
    // Label the strongest markers by fold-change (many share the FDR ceiling, so
    // ranking by y alone piles labels on top of each other). Take the top few by
    // |log2FC| among significant, then thin any that would collide horizontally.
    const strong = [...pts.filter((p) => p.significant)].sort(
      (a, b) => Math.abs(b.log2FC) - Math.abs(a.log2FC),
    )
    const labels = new Set()
    const usedX = []
    const minDx = (2 * (xmax * 1.05)) / 12 // ~12 label slots across the axis
    for (const p of strong) {
      if (labels.size >= 6) break
      if (usedX.some((x) => Math.abs(x - p.log2FC) < minDx)) continue
      labels.add(p.gene)
      usedX.push(p.log2FC)
    }
    if (selectedGene) labels.add(selectedGene)

    const xLim = Math.ceil(xmax)
    return {
      sig: pts.filter((p) => p.significant),
      nonsig: pts.filter((p) => !p.significant),
      xDomain: [-xLim, xLim],
      labelSet: labels,
    }
  }, [genes, yField, selectedGene])

  const renderLabel = (props) => {
    const { x, y, value } = props
    if (!labelSet.has(value)) return null
    return (
      <text
        x={x}
        y={y - 7}
        textAnchor="middle"
        style={{ fontSize: 11, fontStyle: 'italic', fill: COLORS.crimson, fontWeight: 600 }}
      >
        {value}
      </text>
    )
  }

  return (
    <ChartCard
      title="Differential expression"
      sub={`volcano · pathogenic vs reference · y = −log10(${yField})`}
      filename="volcano"
      legend={
        <>
          <span className="lk">
            <span className="swatch" style={{ background: COLORS.crimson }} /> significant (FDR&lt;0.05, |log2FC|&gt;1)
          </span>
          <span className="lk">
            <span className="swatch" style={{ background: COLORS.tealSoft }} /> not significant
          </span>
        </>
      }
      note="Click any point to select that gene. Up-regulated disease markers sit top-right; these are what DSS ranks."
      full
    >
      <ResponsiveContainer width="100%" height={360}>
        <ScatterChart margin={{ top: 16, right: 24, bottom: 32, left: 8 }}>
          <CartesianGrid stroke={COLORS.border} />
          <ReferenceLine x={FC_CUT} stroke={COLORS.borderStrong} strokeDasharray="3 3" />
          <ReferenceLine x={-FC_CUT} stroke={COLORS.borderStrong} strokeDasharray="3 3" />
          <ReferenceLine
            y={-Math.log10(FDR_CUT)}
            stroke={COLORS.borderStrong}
            strokeDasharray="3 3"
          />
          <XAxis
            type="number"
            dataKey="log2FC"
            domain={xDomain}
            allowDecimals={false}
            tick={{ fontSize: 11 }}
            stroke={COLORS.border}
            label={{
              value: 'log2 fold-change  (pathogenic / reference →)',
              position: 'insideBottom',
              offset: -14,
              fontSize: 11,
              fill: COLORS.muted,
            }}
          />
          <YAxis
            type="number"
            dataKey="y"
            tick={{ fontSize: 11 }}
            stroke={COLORS.border}
            label={{
              value: `−log10(${yField})`,
              angle: -90,
              position: 'insideLeft',
              offset: 16,
              fontSize: 11,
              fill: COLORS.muted,
            }}
          />
          <ZAxis type="number" dataKey="z" range={[30, 200]} />
          <Tooltip content={<VolcanoTooltip yField={yField} />} cursor={{ strokeDasharray: '3 3' }} />
          <Scatter
            data={nonsig}
            fill={COLORS.tealSoft}
            fillOpacity={0.5}
            isAnimationActive={false}
            onClick={(d) => d && onSelect(d.gene)}
            style={{ cursor: 'pointer' }}
          />
          <Scatter
            data={sig}
            fill={COLORS.crimson}
            fillOpacity={0.85}
            isAnimationActive={false}
            onClick={(d) => d && onSelect(d.gene)}
            style={{ cursor: 'pointer' }}
          >
            <LabelList dataKey="gene" content={renderLabel} />
          </Scatter>
        </ScatterChart>
      </ResponsiveContainer>
    </ChartCard>
  )
}
