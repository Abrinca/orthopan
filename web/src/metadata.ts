// Genome metadata columns drawn between the tree and the heatmap, one narrow
// column per field, aligned with matrix rows through the shared transform.

import type {Genome, GenomeMetadataInfo} from './api';
import {categoryColor} from './annotations';
import {prepareCanvas, type Transform} from './transform';

export const META_COL_W = 16;
const MISSING = '#f2f4f7';

const NUMERIC_STOPS: [number, number, number][] = [
  [247, 251, 255],
  [198, 219, 239],
  [107, 174, 214],
  [33, 113, 181],
  [8, 48, 107],
];

function colorAt(stops: [number, number, number][], t: number): string {
  const u = Math.min(1, Math.max(0, t)) * (stops.length - 1);
  const a = Math.min(stops.length - 2, Math.floor(u));
  const f = u - a;
  const c = stops[a].map((v, i) => Math.round(v * (1 - f) + stops[a + 1][i] * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}

/** Sequential colour for a numeric field value in [0,1]. */
export function numericColor(t: number): string {
  return colorAt(NUMERIC_STOPS, t);
}

export function metadataValueColor(field: GenomeMetadataInfo, value: string | number | null): string {
  if (value === null || value === undefined || value === '') return MISSING;
  if (field.kind === 'categorical') {
    const idx = (field.categories ?? []).indexOf(String(value));
    return idx < 0 ? MISSING : categoryColor(idx);
  }
  const min = field.min ?? 0;
  const max = field.max ?? 1;
  const span = max - min || 1;
  return numericColor((Number(value) - min) / span);
}

export class Metadata {
  constructor(private canvas: HTMLCanvasElement) {}

  draw(tr: Transform, genomes: Genome[], fields: GenomeMetadataInfo[], pw: number, ph: number): void {
    const ctx = prepareCanvas(this.canvas, pw, ph);
    if (!ctx) return;
    const r0 = Math.max(0, Math.floor(tr.rowAt(0)) - 1);
    const r1 = Math.min(genomes.length - 1, Math.ceil(tr.rowAt(ph)) + 1);
    fields.forEach((field, i) => {
      const x = i * META_COL_W;
      for (let r = r0; r <= r1; r++) {
        const y0 = tr.yOf(r);
        const y1 = tr.yOf(r + 1);
        if (y1 < 0 || y0 > ph) continue;
        ctx.fillStyle = metadataValueColor(field, genomes[r]?.metadata?.[field.name] ?? null);
        ctx.fillRect(x, y0, META_COL_W - 1, Math.max(1, y1 - y0));
      }
    });
    // rotated header labels, drawn downward so narrow columns stay readable
    ctx.font = '10px system-ui, sans-serif';
    ctx.fillStyle = '#475467';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    fields.forEach((field, i) => {
      ctx.save();
      ctx.translate(i * META_COL_W + META_COL_W / 2, 2);
      ctx.rotate(Math.PI / 2);
      ctx.fillText(field.name, 0, 0);
      ctx.restore();
    });
  }

  /** Field and genome row under a pointer (canvas-local CSS px). */
  hit(px: number, py: number, tr: Transform, nGenomes: number): {row: number; field: number} | null {
    const field = Math.floor(px / META_COL_W);
    const row = Math.floor(tr.rowAt(py));
    if (field < 0 || row < 0 || row >= nGenomes) return null;
    return {row, field};
  }
}
