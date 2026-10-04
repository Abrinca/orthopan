// Top strip area: numeric track bars and categorical annotation colour strips,
// aligned with the matrix column axis through the shared transform.

import type {AnnotationTrackResult, TrackInfo, TracksResult} from './api';
import {annotationStripColor} from './annotations';
import {prepareCanvas, type Transform} from './transform';

const PALETTE = ['#2563eb', '#059669', '#d97706'];
const LABEL_H = 13;

export function fmt3(v: number): string {
  if (!Number.isFinite(v)) return '—';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1000 || a < 0.01) return v.toExponential(2);
  return Number(v.toPrecision(3)).toString();
}

/** 1st..99th percentile domain of the finite values, falling back to meta min/max. */
export function robustDomain(values: ArrayLike<number>, info?: TrackInfo): [number, number] {
  const finite: number[] = [];
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (Number.isFinite(v)) finite.push(v);
  }
  const isFlat = (lo: number, hi: number) => !(hi - lo > 1e-9);
  if (finite.length >= 8) {
    finite.sort((a, b) => a - b);
    const lo = finite[Math.floor(0.01 * (finite.length - 1))];
    const hi = finite[Math.ceil(0.99 * (finite.length - 1))];
    if (!isFlat(lo, hi)) return [lo, hi];
  }
  if (info && !isFlat(info.min, info.max)) return [info.min, info.max];
  const v = finite.length ? finite[Math.floor(finite.length / 2)] : 0;
  return [v - 0.5, v + 0.5];
}

export interface NumericStrip {
  kind: 'numeric';
  name: string;
  info?: TrackInfo;
  result: TracksResult;
  world: {og0: number; og1: number};
}

export interface CategoryStrip {
  kind: 'category';
  name: string;
  categories: string[];
  result: AnnotationTrackResult;
  world: {og0: number; og1: number};
}

export type Strip = NumericStrip | CategoryStrip;

export interface StripHit {
  strip: Strip;
  index: number;
  value: number;
}

export function stripValue(strip: Strip, i: number): number {
  return strip.kind === 'numeric'
    ? strip.result.data[strip.result.names.indexOf(strip.name) * strip.result.bins + i]
    : strip.result.data[i];
}

export class Tracks {
  constructor(private canvas: HTMLCanvasElement) {}

  draw(tr: Transform, strips: Strip[], pw: number, ph: number): void {
    const ctx = prepareCanvas(this.canvas, pw, ph);
    if (!ctx) return;
    const n = strips.length;
    if (!n) return;
    const band = ph / n;
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, pw, ph);
    ctx.clip();
    strips.forEach((strip, t) => {
      const top = t * band;
      ctx.fillStyle = t % 2 ? '#fbfcfe' : '#f5f8fc';
      ctx.fillRect(0, top, pw, band);
      ctx.strokeStyle = '#e4e7ec';
      ctx.beginPath();
      ctx.moveTo(0, top + 0.5);
      ctx.lineTo(pw, top + 0.5);
      ctx.stroke();
      if (strip.kind === 'category') this.drawCategory(ctx, tr, strip, top, band, pw);
      else this.drawNumeric(ctx, tr, strip, top, band, pw, t);
    });
    ctx.restore();
  }

  private drawNumeric(
    ctx: CanvasRenderingContext2D,
    tr: Transform,
    strip: NumericStrip,
    top: number,
    band: number,
    pw: number,
    t: number,
  ): void {
    const L = strip.world.og1 - strip.world.og0;
    const src = strip.result.names.indexOf(strip.name);
    if (src < 0) return;
    const values = strip.result.data.subarray(src * strip.result.bins, (src + 1) * strip.result.bins);
    const [lo, hi] = robustDomain(values, strip.info);
    const span = hi - lo || 1;
    ctx.fillStyle = PALETTE[t % PALETTE.length];
    const base = top + band - 3;
    const maxBar = band - LABEL_H - 2;
    for (let i = 0; i < strip.result.bins; i++) {
      const v = values[i];
      if (!Number.isFinite(v)) continue;
      const x0 = tr.xOf(strip.world.og0 + (i * L) / strip.result.bins);
      const x1 = tr.xOf(strip.world.og0 + ((i + 1) * L) / strip.result.bins);
      if (x1 < 0 || x0 > pw) continue;
      const r = Math.min(1, Math.max(0, (v - lo) / span));
      ctx.fillRect(x0, base - Math.max(1, r * maxBar), Math.max(1, x1 - x0), Math.max(1, r * maxBar));
    }
    ctx.fillStyle = '#344054';
    ctx.font = '10px system-ui, sans-serif';
    ctx.textBaseline = 'top';
    ctx.fillText(`${strip.name}  ${fmt3(lo)}–${fmt3(hi)}`, 4, top + 1);
  }

  private drawCategory(
    ctx: CanvasRenderingContext2D,
    tr: Transform,
    strip: CategoryStrip,
    top: number,
    band: number,
    pw: number,
  ): void {
    const L = strip.world.og1 - strip.world.og0;
    const y = top + LABEL_H;
    const h = Math.max(1, band - LABEL_H - 1);
    for (let i = 0; i < strip.result.bins; i++) {
      const v = strip.result.data[i];
      if (v <= 0) continue;
      const x0 = tr.xOf(strip.world.og0 + (i * L) / strip.result.bins);
      const x1 = tr.xOf(strip.world.og0 + ((i + 1) * L) / strip.result.bins);
      if (x1 < 0 || x0 > pw) continue;
      ctx.fillStyle = annotationStripColor(v);
      ctx.fillRect(x0, y, Math.max(1, x1 - x0), h);
    }
    ctx.fillStyle = '#344054';
    ctx.font = '10px system-ui, sans-serif';
    ctx.textBaseline = 'top';
    ctx.fillText(strip.name, 4, top + 1);
  }

  /** Strip/bin under a pointer (canvas-local CSS px), or null if outside. */
  hit(tr: Transform, strips: Strip[], px: number, py: number, ph: number): StripHit | null {
    if (!strips.length || py < 0 || py >= ph) return null;
    const band = ph / strips.length;
    const t = Math.min(strips.length - 1, Math.max(0, Math.floor(py / band)));
    const strip = strips[t];
    const L = strip.world.og1 - strip.world.og0;
    const op = tr.ogAt(px);
    // outside the fetched orthogroup span (e.g. a stale strip) -> no hit
    if (op < strip.world.og0 || op >= strip.world.og1) return null;
    const i = Math.min(strip.result.bins - 1, Math.max(0, Math.floor(((op - strip.world.og0) * strip.result.bins) / L)));
    return {strip, index: i, value: stripValue(strip, i)};
  }
}
