import { useEffect, useMemo, useRef, useState } from 'react'
import { COLORS, PALETTE, reds, fmt } from '../lib/utils'

// UMAP embedding scatter, rendered on <canvas> for 8k-point performance.
// Color-by toggle: population (P crimson / others teal-grey), cell type
// (categorical), or expression of a chosen gene (Reds colormap).
// The gene options are the expression columns embedded in the umap json.

// Population -> color. Pathogenic pops in crimson family, references in teal/grey.
const POP_COLORS = {
  P: COLORS.crimson,
  B: '#9aa0a0',
  H: COLORS.tealSoft,
  R: '#c0a35a',
}
const POP_LABEL = { P: 'Pathogenic', B: 'Background', H: 'Healthy', R: 'Related' }

const CANVAS_W = 640
const CANVAS_H = 460
const PAD = 18

export default function UmapPanel({ points, geneCols, selectedGene }) {
  const canvasRef = useRef(null)
  const [colorBy, setColorBy] = useState('pop') // 'pop' | 'cell_type' | 'expr'
  // default the expression gene to the selected gene if it's in the umap, else first col
  const exprGene = useMemo(() => {
    if (selectedGene && geneCols.includes(selectedGene)) return selectedGene
    return geneCols[0] || null
  }, [selectedGene, geneCols])
  const [exprCol, setExprCol] = useState(exprGene)
  useEffect(() => setExprCol(exprGene), [exprGene])

  const pops = useMemo(() => [...new Set(points.map((p) => p.pop))], [points])
  const cellTypes = useMemo(() => [...new Set(points.map((p) => p.cell_type))], [points])
  const cellTypeColor = useMemo(() => {
    const m = new Map()
    cellTypes.forEach((ct, i) => m.set(ct, PALETTE[i % PALETTE.length]))
    return m
  }, [cellTypes])

  // coordinate transform (data -> canvas), memoized on the point cloud
  const transform = useMemo(() => {
    let x0 = Infinity
    let x1 = -Infinity
    let y0 = Infinity
    let y1 = -Infinity
    for (const p of points) {
      if (p.UMAP1 < x0) x0 = p.UMAP1
      if (p.UMAP1 > x1) x1 = p.UMAP1
      if (p.UMAP2 < y0) y0 = p.UMAP2
      if (p.UMAP2 > y1) y1 = p.UMAP2
    }
    const sx = (CANVAS_W - 2 * PAD) / (x1 - x0 || 1)
    const sy = (CANVAS_H - 2 * PAD) / (y1 - y0 || 1)
    return {
      px: (x) => PAD + (x - x0) * sx,
      // flip y so +UMAP2 is up
      py: (y) => CANVAS_H - PAD - (y - y0) * sy,
    }
  }, [points])

  const exprMax = useMemo(() => {
    if (colorBy !== 'expr' || !exprCol) return 1
    let m = 0
    for (const p of points) {
      const v = p[exprCol]
      if (typeof v === 'number' && v > m) m = v
    }
    return m || 1
  }, [points, exprCol, colorBy])

  // draw
  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return
    const dpr = Math.min(window.devicePixelRatio || 1, 2)
    canvas.width = CANVAS_W * dpr
    canvas.height = CANVAS_H * dpr
    const ctx = canvas.getContext('2d')
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
    ctx.clearRect(0, 0, CANVAS_W, CANVAS_H)

    const colorOf = (p) => {
      if (colorBy === 'pop') return POP_COLORS[p.pop] || COLORS.muted
      if (colorBy === 'cell_type') return cellTypeColor.get(p.cell_type) || COLORS.muted
      return reds((p[exprCol] || 0) / exprMax)
    }

    // For expression view, draw low-expression points first so highlights sit on top.
    const order =
      colorBy === 'expr' && exprCol
        ? [...points].sort((a, b) => (a[exprCol] || 0) - (b[exprCol] || 0))
        : points

    ctx.globalAlpha = 0.72
    for (const p of order) {
      ctx.fillStyle = colorOf(p)
      ctx.beginPath()
      ctx.arc(transform.px(p.UMAP1), transform.py(p.UMAP2), 2.1, 0, 2 * Math.PI)
      ctx.fill()
    }
    ctx.globalAlpha = 1
  }, [points, transform, colorBy, cellTypeColor, exprCol, exprMax])

  const modes = [
    { k: 'pop', label: 'Population' },
    { k: 'cell_type', label: 'Cell type' },
    { k: 'expr', label: 'Expression' },
  ]

  return (
    <div className="chart-card full">
      <div className="chart-head">
        <h3>
          UMAP embedding
          <span className="sub">  {points.length.toLocaleString()} cells · color by {modes.find((m) => m.k === colorBy).label.toLowerCase()}</span>
        </h3>
        <div className="seg">
          {modes.map((m) => (
            <button
              key={m.k}
              className={`seg-btn${colorBy === m.k ? ' active' : ''}`}
              onClick={() => setColorBy(m.k)}
              disabled={m.k === 'cell_type' && cellTypes.length < 2}
              title={m.k === 'cell_type' && cellTypes.length < 2 ? 'single cell type in this cohort' : undefined}
            >
              {m.label}
            </button>
          ))}
        </div>
      </div>

      <div className="chart-body" style={{ display: 'flex', gap: 16, alignItems: 'flex-start', flexWrap: 'wrap' }}>
        <div style={{ flex: '1 1 380px', minWidth: 300 }}>
          <canvas
            ref={canvasRef}
            style={{ width: '100%', maxWidth: CANVAS_W, height: 'auto', aspectRatio: `${CANVAS_W} / ${CANVAS_H}`, display: 'block' }}
          />
        </div>

        {/* legend / controls */}
        <div className="umap-legend">
          {colorBy === 'pop' ? (
            <div className="legend-list">
              {pops.map((p) => (
                <span className="lk" key={p}>
                  <span className="swatch" style={{ background: POP_COLORS[p] || COLORS.muted }} />
                  {POP_LABEL[p] || p}
                </span>
              ))}
            </div>
          ) : null}

          {colorBy === 'cell_type' ? (
            <div className="legend-list">
              {cellTypes.map((ct) => (
                <span className="lk" key={ct} title={ct}>
                  <span className="swatch" style={{ background: cellTypeColor.get(ct) }} />
                  <span className="ct-name">{ct}</span>
                </span>
              ))}
            </div>
          ) : null}

          {colorBy === 'expr' ? (
            <div>
              <label className="umap-gene-label" htmlFor="umap-gene">Gene</label>
              <select
                id="umap-gene"
                value={exprCol || ''}
                onChange={(e) => setExprCol(e.target.value)}
                style={{ minWidth: 130, marginBottom: 12 }}
              >
                {geneCols.map((g) => (
                  <option key={g} value={g}>
                    {g}
                  </option>
                ))}
              </select>
              <div className="expr-ramp">
                <span className="ramp-bar" />
                <div className="ramp-ticks">
                  <span>0</span>
                  <span>{fmt(exprMax, 1)}</span>
                </div>
                <span className="ramp-cap">{exprCol} expr</span>
              </div>
            </div>
          ) : null}
        </div>
      </div>
      <div className="chart-note">
        Subsampled embedding ({points.length.toLocaleString()} of the full cohort). Pathogenic cells sit in crimson under the population view; switch to Expression to see where a candidate gene concentrates.
      </div>
    </div>
  )
}
