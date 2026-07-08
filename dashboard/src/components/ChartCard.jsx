import { useRef } from 'react'
import { downloadChartSvg, downloadChartPng } from '../lib/utils'

// Wraps a recharts chart with a titled header and SVG/PNG download buttons.
// The chart is rendered as children inside a ref'd container so we can reach
// into the DOM to serialize the <svg> recharts produces.
export default function ChartCard({ title, sub, filename, children, note, legend, full = false }) {
  const ref = useRef(null)
  const base = filename || (title || 'chart').toLowerCase().replace(/[^a-z0-9]+/g, '_')

  return (
    <div className={`chart-card${full ? ' full' : ''}`}>
      <div className="chart-head">
        <h3>
          {title}
          {sub ? <span className="sub">  {sub}</span> : null}
        </h3>
        <div className="nowrap">
          <button
            className="dl-btn"
            onClick={() => downloadChartSvg(ref.current, `${base}.svg`)}
            title="Download as SVG"
          >
            SVG
          </button>{' '}
          <button
            className="dl-btn"
            onClick={() => downloadChartPng(ref.current, `${base}.png`)}
            title="Download as PNG"
          >
            PNG
          </button>
        </div>
      </div>
      <div className="chart-body" ref={ref}>
        {children}
      </div>
      {legend ? <div className="legend-inline">{legend}</div> : null}
      {note ? <div className="chart-note">{note}</div> : null}
    </div>
  )
}
