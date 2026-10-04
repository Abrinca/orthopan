// Top track strips (per-orthogroup values) aligned with the matrix column axis.

import type {TrackInfo, TracksResult} from './api';
import type {Transform} from './transform';

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

export class Tracks {
  private readonly ctx: CanvasRenderingContext2D;

  constructor(private canvas: HTMLCanvasElement) {
    const ctx = canvas.getContext('2d');
    if (!ctx) throw new Error('tracks: no 2d context');
    this.ctx = ctx;
  }

  draw(
    tr: Transform,
    result: TracksResult,
    world: {og0: number; og1: number},
    infos: TrackInfo[],
    requested: string[],
    pw: number,
    ph: number,
  ): void {
    const dpr = window.devicePixelRatio || 1;
    if (this.canvas.width !== Math.round(pw * dpr) || this.canvas.height !== Math.round(ph * dpr)) {
      this.canvas.width = Math.round(pw * dpr);
      this.canvas.height = Math.round(ph * dpr);
    }
    const ctx = this.ctx;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, pw, ph);

    const names = requested.filter((n) => result.names.includes(n));
    const nT = names.length;
    if (!nT) return;
    const L = world.og1 - world.og0;
    const band = ph / nT;
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, pw, ph);
    ctx.clip();
    names.forEach((name, t) => {
      const src = result.names.indexOf(name);
      const values = result.data.subarray(src * result.bins, (src + 1) * result.bins);
      const info = infos.find((x) => x.name === name);
      const [lo, hi] = robustDomain(values, info);
      const span = hi - lo || 1;
      const top = t * band;
      ctx.fillStyle = t % 2 ? '#fbfcfe' : '#f5f8fc';
      ctx.fillRect(0, top, pw, band);
      ctx.strokeStyle = '#e4e7ec';
      ctx.beginPath();
      ctx.moveTo(0, top + 0.5);
      ctx.lineTo(pw, top + 0.5);
      ctx.stroke();

      const color = PALETTE[t % PALETTE.length];
      ctx.fillStyle = color;
      const base = top + band - 3;
      const maxBar = band - LABEL_H - 2;
      for (let i = 0; i < result.bins; i++) {
        const v = values[i];
        if (!Number.isFinite(v)) continue;
        const x0 = tr.xOf(world.og0 + (i * L) / result.bins);
        const x1 = tr.xOf(world.og0 + ((i + 1) * L) / result.bins);
        if (x1 < 0 || x0 > pw) continue;
        const r = Math.min(1, Math.max(0, (v - lo) / span));
        const barH = Math.max(1, r * maxBar);
        ctx.fillRect(x0, base - barH, Math.max(1, x1 - x0), barH);
      }
      ctx.fillStyle = '#344054';
      ctx.font = '10px system-ui, sans-serif';
      ctx.textBaseline = 'top';
      ctx.fillText(`${name}  ${fmt3(lo)}–${fmt3(hi)}`, 4, top + 1);
    });
    ctx.restore();
  }
}
