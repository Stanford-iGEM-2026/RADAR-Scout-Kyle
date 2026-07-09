// Shared helpers: brand colors, formatting, CSV export, chart-SVG/PNG download.

export const COLORS = {
  crimson: '#8e1918',
  crimsonSoft: '#b45150',
  teal: '#1c7170',
  tealSoft: '#5aa3a2',
  tealFaint: '#cfe3e3',
  ink: '#141414',
  muted: '#6b6b6b',
  border: '#e6e6e6',
  borderStrong: '#d5d5d5',
}

// Component fields shown across list / filters / compare.
export const COMPONENTS = [
  { key: 'RACS', label: 'RACS' },
  { key: 'Sep', label: 'Sep' },
  { key: 'Feas', label: 'Feas' },
  { key: 'Repro', label: 'Repro' },
  { key: 'OffMax', label: 'OffMax' },
]

// The two ways a target can be ranked.
export const RANK_METRICS = {
  RACS: {
    key: 'RACS',
    label: 'RACS',
    long: 'RADAR compatibility',
    blurb: 'good RADAR targets',
    d: 3,
    color: '#8e1918',
  },
  DSS: {
    key: 'DSS',
    label: 'DSS',
    long: 'Disease specificity',
    blurb: 'markers of this disease',
    d: 2,
    color: '#1c7170',
  },
  spec_score: {
    key: 'spec_score',
    label: 'Sensor',
    long: 'Detection specificity',
    blurb: 'near-binary sensor markers (POSTN/ADAM12)',
    d: 2,
    color: '#8e1918',
  },
}

// Categorical palette for cell types etc. — brand-forward, muted, no neon.
export const PALETTE = [
  '#8e1918', '#1c7170', '#b45150', '#5aa3a2', '#c98a2b',
  '#6b6b6b', '#4a6d8c', '#8a6ea3', '#a3a3a3', '#3f7d54',
  '#c25b7a', '#7a8a3f', '#2f5d5c', '#b8823a', '#555555',
]

export const fmt = (v, d = 3) =>
  v === null || v === undefined || Number.isNaN(v) ? '—' : Number(v).toFixed(d)

// True if the field exists (and isn't null) on the first record.
export const hasCol = (rows, key) =>
  Array.isArray(rows) && rows.length > 0 && rows[0][key] !== undefined && rows[0][key] !== null

// Reds-style sequential colormap on t in [0,1] → 'rgb(...)'. White → crimson.
export function reds(t) {
  const x = Math.max(0, Math.min(1, Number.isFinite(t) ? t : 0))
  // interpolate #fff5f0 (near-white) → #8e1918 (crimson), gamma-ish for punch
  const g = Math.pow(x, 0.85)
  const r = Math.round(255 + (142 - 255) * g)
  const gg = Math.round(245 + (25 - 245) * g)
  const b = Math.round(240 + (24 - 240) * g)
  return `rgb(${r},${gg},${b})`
}

// Min/max of a numeric field, robust to missing values.
export function extent(rows, key) {
  let lo = Infinity
  let hi = -Infinity
  for (const r of rows) {
    const v = r[key]
    if (typeof v === 'number' && Number.isFinite(v)) {
      if (v < lo) lo = v
      if (v > hi) hi = v
    }
  }
  if (lo === Infinity) return [0, 1]
  if (lo === hi) return [lo - 0.5, hi + 0.5]
  return [lo, hi]
}

// Trigger a browser download of a Blob.
function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url
  a.download = filename
  document.body.appendChild(a)
  a.click()
  document.body.removeChild(a)
  URL.revokeObjectURL(url)
}

// Export an array of gene objects to CSV (respects a column order).
export function exportCsv(rows, filename = 'racs_table.csv') {
  if (!rows.length) return
  const cols = [
    'gene', 'RACS', 'DSS', 'Sep', 'Feas', 'Repro', 'OffMax',
    'detect_P', 'log2FC', 'FDR', 'celltype_spec', 'disease_spec',
    'k_op', 'Youden_J', 'n_donors', 'act_P', 'act_H', 'act_B', 'act_R',
  ]
  const present = cols.filter((c) => c in rows[0])
  const header = present.join(',')
  const body = rows
    .map((r) =>
      present
        .map((c) => {
          const val = r[c]
          if (typeof val === 'string' && (val.includes(',') || val.includes('"'))) {
            return `"${val.replace(/"/g, '""')}"`
          }
          return val ?? ''
        })
        .join(','),
    )
    .join('\n')
  downloadBlob(new Blob([`${header}\n${body}\n`], { type: 'text/csv;charset=utf-8' }), filename)
}

// Locate the <svg> a recharts container rendered inside `container`, clone it,
// inline a white background + font, and hand back a serialized standalone SVG.
function serializeChartSvg(container) {
  const svg = container?.querySelector('svg')
  if (!svg) return null
  const clone = svg.cloneNode(true)
  const width = svg.clientWidth || svg.getBoundingClientRect().width || 600
  const height = svg.clientHeight || svg.getBoundingClientRect().height || 360
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
  clone.setAttribute('width', width)
  clone.setAttribute('height', height)
  clone.setAttribute('viewBox', `0 0 ${width} ${height}`)
  // white background rect so exported file isn't transparent
  const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect')
  rect.setAttribute('width', '100%')
  rect.setAttribute('height', '100%')
  rect.setAttribute('fill', '#ffffff')
  clone.insertBefore(rect, clone.firstChild)
  const style = document.createElementNS('http://www.w3.org/2000/svg', 'style')
  style.textContent =
    "text{font-family:-apple-system,BlinkMacSystemFont,'Inter',system-ui,sans-serif;}"
  clone.insertBefore(style, clone.firstChild)
  return { markup: new XMLSerializer().serializeToString(clone), width, height }
}

export function downloadChartSvg(container, filename = 'chart.svg') {
  const out = serializeChartSvg(container)
  if (!out) return
  downloadBlob(new Blob([out.markup], { type: 'image/svg+xml;charset=utf-8' }), filename)
}

// Rasterize the chart SVG to PNG via a canvas (2x for crispness).
export function downloadChartPng(container, filename = 'chart.png') {
  const out = serializeChartSvg(container)
  if (!out) return
  const { markup, width, height } = out
  const scale = 2
  const img = new Image()
  const svgBlob = new Blob([markup], { type: 'image/svg+xml;charset=utf-8' })
  const url = URL.createObjectURL(svgBlob)
  img.onload = () => {
    const canvas = document.createElement('canvas')
    canvas.width = width * scale
    canvas.height = height * scale
    const ctx = canvas.getContext('2d')
    ctx.scale(scale, scale)
    ctx.drawImage(img, 0, 0)
    URL.revokeObjectURL(url)
    canvas.toBlob((blob) => {
      if (blob) downloadBlob(blob, filename)
    }, 'image/png')
  }
  img.onerror = () => URL.revokeObjectURL(url)
  img.src = url
}
