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

export const fmt = (v, d = 3) =>
  v === null || v === undefined || Number.isNaN(v) ? '—' : Number(v).toFixed(d)

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
    'gene', 'RACS', 'Sep', 'Feas', 'Repro', 'OffMax',
    'k_op', 'Youden_J', 'n_donors', 'act_P', 'act_H', 'act_R',
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
