import {
  BarChart,
  Bar,
  Cell,
  XAxis,
  YAxis,
  CartesianGrid,
  ResponsiveContainer,
  Tooltip,
  ScatterChart,
  Scatter,
  ZAxis,
  ReferenceArea,
  LabelList,
} from 'recharts'
import ChartCard from './ChartCard'
import { COLORS, fmt } from '../lib/utils'

// ---- custom tooltips (editorial, no recharts default chrome) -------------- //
function BreakdownTooltip({ active, payload }) {
  if (!active || !payload?.length) return null
  const p = payload[0].payload
  return (
    <div className="rc-tooltip">
      <div className="tt-title">{p.name}</div>
      <div className="tt-row">
        <span className="tt-k">value</span>
        <span>{fmt(p.value, 3)}</span>
      </div>
    </div>
  )
}

function ActTooltip({ active, payload, label }) {
  if (!active || !payload?.length) return null
  return (
    <div className="rc-tooltip">
      <div className="tt-title">{label}</div>
      {payload.map((s) => (
        <div className="tt-row" key={s.dataKey}>
          <span className="tt-k" style={{ color: s.color }}>
            {s.name}
          </span>
          <span>{fmt(s.value, 3)}</span>
        </div>
      ))}
    </div>
  )
}

function ScatterTooltip({ active, payload }) {
  if (!active || !payload?.length) return null
  const p = payload[0].payload
  return (
    <div className="rc-tooltip">
      <div className="tt-title">{p.gene}</div>
      <div className="tt-row">
        <span className="tt-k">Feas (abundance)</span>
        <span>{fmt(p.Feas, 3)}</span>
      </div>
      <div className="tt-row">
        <span className="tt-k">Sep (specificity)</span>
        <span>{fmt(p.Sep, 3)}</span>
      </div>
      <div className="tt-row">
        <span className="tt-k">RACS</span>
        <span>{fmt(p.RACS, 3)}</span>
      </div>
    </div>
  )
}

const Swatch = ({ color }) => <span className="swatch" style={{ background: color }} />

export default function DetailPanel({ gene, allGenes, rank }) {
  if (!gene) {
    return (
      <div className="detail-empty">
        <p style={{ margin: 0, fontSize: 15 }}>Select a gene to inspect its RACS profile.</p>
        <p className="muted" style={{ marginTop: 8, fontSize: 13 }}>
          Breakdown, therapeutic-window activation, and the abundance–specificity sweet spot.
        </p>
      </div>
    )
  }

  // (a) RACS breakdown: the four multiplicative factors, RACS highlighted.
  const breakdown = [
    { name: 'Sep', value: gene.Sep, kind: 'comp' },
    { name: 'Feas', value: gene.Feas, kind: 'comp' },
    { name: 'Repro', value: gene.Repro, kind: 'comp' },
    { name: '1 − OffMax', value: 1 - gene.OffMax, kind: 'comp' },
    { name: 'RACS', value: gene.RACS, kind: 'racs' },
  ]

  // (b) per-population activation (therapeutic window).
  const activation = [
    { pop: 'Pathogenic', key: 'P', P: gene.act_P, off: null },
    { pop: 'Healthy', key: 'H', P: null, off: gene.act_H },
    { pop: 'Related', key: 'R', P: null, off: gene.act_R },
  ]

  // (c) abundance-vs-specificity scatter for ALL genes.
  const scatterData = allGenes.map((g) => ({
    gene: g.gene,
    Feas: g.Feas,
    Sep: g.Sep,
    RACS: g.RACS,
    z: Math.max(30, g.RACS * 340),
    selected: g.gene === gene.gene,
  }))
  const others = scatterData.filter((d) => !d.selected)
  const selPoint = scatterData.filter((d) => d.selected)

  // (d) off-target view: on-target vs worst-case off-target.
  const offMaxPop = gene.act_H >= gene.act_R ? 'H' : 'R'
  const offTarget = [
    {
      name: gene.gene,
      onTarget: gene.act_P,
      offTarget: Math.max(gene.act_H, gene.act_R),
    },
  ]

  return (
    <div>
      <div className="detail-header">
        <div className="detail-title">
          <h2>{gene.gene}</h2>
          <div className="detail-racs">
            <span className="num">{fmt(gene.RACS, 3)}</span>
            <span className="cap">RACS</span>
          </div>
          <span className="rank-badge">rank #{rank}</span>
        </div>
        <div className="kv-strip">
          <div className="kv">
            <span className="k">k_op</span>
            <span className="v">{fmt(gene.k_op, 1)}</span>
          </div>
          <div className="kv">
            <span className="k">Youden J</span>
            <span className="v">{fmt(gene.Youden_J, 2)}</span>
          </div>
          <div className="kv">
            <span className="k">Donors</span>
            <span className="v">{gene.n_donors}</span>
          </div>
        </div>
      </div>

      <div className="charts-grid">
        {/* (a) RACS breakdown */}
        <ChartCard
          title="RACS breakdown"
          sub="multiplicative factors"
          filename={`${gene.gene}_racs_breakdown`}
          note="RACS = Sep · Feas · Repro · (1 − OffMax). Teal = components; crimson = composite."
        >
          <ResponsiveContainer width="100%" height={210}>
            <BarChart
              data={breakdown}
              layout="vertical"
              margin={{ top: 4, right: 44, bottom: 4, left: 8 }}
              barCategoryGap={8}
            >
              <CartesianGrid horizontal={false} stroke={COLORS.border} />
              <XAxis type="number" domain={[0, 1]} tick={{ fontSize: 11 }} stroke={COLORS.border} />
              <YAxis
                type="category"
                dataKey="name"
                width={78}
                tick={{ fontSize: 12 }}
                stroke={COLORS.border}
              />
              <Tooltip content={<BreakdownTooltip />} cursor={{ fill: 'rgba(0,0,0,0.03)' }} />
              <Bar dataKey="value" radius={[0, 2, 2, 0]} isAnimationActive={false}>
                {breakdown.map((d) => (
                  <Cell key={d.name} fill={d.kind === 'racs' ? COLORS.crimson : COLORS.teal} />
                ))}
                <LabelList
                  dataKey="value"
                  position="right"
                  formatter={(v) => fmt(v, 2)}
                  style={{ fontSize: 11, fill: COLORS.muted }}
                />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </ChartCard>

        {/* (b) per-population activation */}
        <ChartCard
          title="Per-population activation"
          sub={`therapeutic window · k_op = ${fmt(gene.k_op, 1)}`}
          filename={`${gene.gene}_activation_window`}
          legend={
            <>
              <span className="lk">
                <Swatch color={COLORS.crimson} /> Pathogenic (on-target)
              </span>
              <span className="lk">
                <Swatch color={COLORS.teal} /> Off-target (healthy / related)
              </span>
            </>
          }
          note="A wide gap (crimson high, teal low) is the therapeutic window RADAR exploits."
        >
          <ResponsiveContainer width="100%" height={210}>
            <BarChart data={activation} margin={{ top: 8, right: 12, bottom: 4, left: -12 }}>
              <CartesianGrid vertical={false} stroke={COLORS.border} />
              <XAxis dataKey="pop" tick={{ fontSize: 12 }} stroke={COLORS.border} />
              <YAxis domain={[0, 1]} tick={{ fontSize: 11 }} stroke={COLORS.border} />
              <Tooltip content={<ActTooltip />} cursor={{ fill: 'rgba(0,0,0,0.03)' }} />
              <Bar
                dataKey="P"
                name="On-target"
                fill={COLORS.crimson}
                radius={[2, 2, 0, 0]}
                isAnimationActive={false}
                maxBarSize={64}
              />
              <Bar
                dataKey="off"
                name="Off-target"
                fill={COLORS.teal}
                radius={[2, 2, 0, 0]}
                isAnimationActive={false}
                maxBarSize={64}
              />
            </BarChart>
          </ResponsiveContainer>
        </ChartCard>

        {/* (c) abundance-vs-specificity scatter (full width) */}
        <ChartCard
          title="Abundance vs specificity"
          sub="the RADAR sweet spot · point size ∝ RACS"
          filename={`${gene.gene}_sweet_spot`}
          full
          note="High-Sep / high-Feas is the knee where a target is both specific and abundant enough to fire. The selected gene is crimson."
        >
          <ResponsiveContainer width="100%" height={340}>
            <ScatterChart margin={{ top: 16, right: 24, bottom: 32, left: 8 }}>
              <CartesianGrid stroke={COLORS.border} />
              {/* sweet-spot annotation */}
              <ReferenceArea
                x1={0.6}
                x2={1}
                y1={0.75}
                y2={1}
                fill={COLORS.tealFaint}
                fillOpacity={0.35}
                stroke={COLORS.tealSoft}
                strokeDasharray="3 3"
                label={{
                  value: 'RADAR sweet spot',
                  position: 'insideTopRight',
                  fontSize: 11,
                  fill: COLORS.teal,
                }}
              />
              <XAxis
                type="number"
                dataKey="Feas"
                name="Feas"
                domain={[0, 1]}
                tick={{ fontSize: 11 }}
                stroke={COLORS.border}
                label={{
                  value: 'Feas  (abundance proxy →)',
                  position: 'insideBottom',
                  offset: -14,
                  fontSize: 11,
                  fill: COLORS.muted,
                }}
              />
              <YAxis
                type="number"
                dataKey="Sep"
                name="Sep"
                domain={[0, 1]}
                tick={{ fontSize: 11 }}
                stroke={COLORS.border}
                label={{
                  value: 'Sep  (specificity →)',
                  angle: -90,
                  position: 'insideLeft',
                  offset: 16,
                  fontSize: 11,
                  fill: COLORS.muted,
                }}
              />
              <ZAxis type="number" dataKey="z" range={[30, 340]} />
              <Tooltip content={<ScatterTooltip />} cursor={{ strokeDasharray: '3 3' }} />
              <Scatter data={others} fill={COLORS.tealSoft} fillOpacity={0.55} isAnimationActive={false} />
              <Scatter data={selPoint} fill={COLORS.crimson} isAnimationActive={false}>
                <LabelList
                  dataKey="gene"
                  position="top"
                  style={{ fontSize: 12, fontStyle: 'italic', fill: COLORS.crimson, fontWeight: 600 }}
                />
              </Scatter>
            </ScatterChart>
          </ResponsiveContainer>
        </ChartCard>

        {/* (d) off-target view */}
        <ChartCard
          title="Off-target check"
          sub={`on-target vs max(act_H, act_R) — worst is ${offMaxPop}`}
          filename={`${gene.gene}_off_target`}
          legend={
            <>
              <span className="lk">
                <Swatch color={COLORS.crimson} /> On-target act_P
              </span>
              <span className="lk">
                <Swatch color={COLORS.teal} /> Worst off-target
              </span>
            </>
          }
          note="Compact leakage view: the taller the crimson over teal, the cleaner the target."
        >
          <ResponsiveContainer width="100%" height={210}>
            <BarChart data={offTarget} margin={{ top: 8, right: 12, bottom: 4, left: -12 }}>
              <CartesianGrid vertical={false} stroke={COLORS.border} />
              <XAxis dataKey="name" tick={{ fontSize: 12, fontStyle: 'italic' }} stroke={COLORS.border} />
              <YAxis domain={[0, 1]} tick={{ fontSize: 11 }} stroke={COLORS.border} />
              <Tooltip content={<ActTooltip />} cursor={{ fill: 'rgba(0,0,0,0.03)' }} />
              <Bar
                dataKey="onTarget"
                name="On-target"
                fill={COLORS.crimson}
                radius={[2, 2, 0, 0]}
                isAnimationActive={false}
                maxBarSize={70}
              />
              <Bar
                dataKey="offTarget"
                name="Worst off-target"
                fill={COLORS.teal}
                radius={[2, 2, 0, 0]}
                isAnimationActive={false}
                maxBarSize={70}
              />
            </BarChart>
          </ResponsiveContainer>
        </ChartCard>
      </div>
    </div>
  )
}
