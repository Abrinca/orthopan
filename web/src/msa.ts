// MSA viewer: pure layout/colour/conservation math plus a dependency-free
// Canvas2D renderer that draws only the visible window.

import type {MsaResult, MsaRow} from './api';
import {prepareCanvas} from './transform';

// --- pure math --------------------------------------------------------------

export const MIN_PX_PER_COL = 0.25;
export const MAX_PX_PER_COL = 40;
/** below this column width residues are drawn as blocks, not letters */
export const LETTER_PX = 8;

export const clampPxPerCol = (px: number) => Math.min(MAX_PX_PER_COL, Math.max(MIN_PX_PER_COL, px));
export const showLetters = (pxPerCol: number) => pxPerCol >= LETTER_PX;

/** Half-open visible column range [c0, c1) for a horizontal scroll/zoom state. */
export function visibleCols(scrollCol: number, width: number, pxPerCol: number, length: number): [number, number] {
  const c0 = Math.max(0, Math.floor(scrollCol));
  const c1 = Math.min(length, Math.ceil(scrollCol + width / pxPerCol));
  return [c0, Math.max(c0, c1)];
}

/** Keep the horizontal scroll inside [0, length - viewportCols]. */
export function clampScrollCol(scrollCol: number, width: number, pxPerCol: number, length: number): number {
  const viewportCols = width / pxPerCol;
  return Math.max(0, Math.min(scrollCol, Math.max(0, length - viewportCols)));
}

/** Keep the vertical scroll inside [0, contentH - viewH]. */
export function clampScrollY(scrollY: number, contentH: number, viewH: number): number {
  return Math.max(0, Math.min(scrollY, Math.max(0, contentH - viewH)));
}

/**
 * Zoom keeping the anchor column under its pixel offset: screen x of the anchor
 * stays `anchorPx`, i.e. newScroll = anchorCol - anchorPx / newPxPerCol.
 */
export function zoomedScroll(anchorCol: number, anchorPx: number, newPxPerCol: number, width: number, length: number): number {
  return clampScrollCol(anchorCol - anchorPx / newPxPerCol, width, newPxPerCol, length);
}

/** Default zoom: whole alignment fits, capped at 12 px/col. */
export function fitPxPerCol(width: number, length: number): number {
  return clampPxPerCol(Math.min(12, width / Math.max(1, length)));
}

/** Horizontal scroll from clicking the overview strip at fractional position. */
export function scrollFromBar(frac: number, width: number, pxPerCol: number, length: number): number {
  const viewportCols = width / pxPerCol;
  return clampScrollCol(Math.min(1, Math.max(0, frac)) * length - viewportCols / 2, width, pxPerCol, length);
}

// --- residue colours (Clustal-like) -----------------------------------------

const HYDROPHOBIC = new Set('AVFPMILW'.split(''));
const ACIDIC = new Set('DE'.split(''));
const BASIC = new Set('RK'.split(''));
const POLAR = new Set('STYHCNGQ'.split(''));

export const GAP_COLOR = '#eceff3';

/** Clustal-like residue background colour; gaps are light. */
export function residueColor(residue: string): string {
  if (residue === '-' || residue === '.' || residue === ' ') return GAP_COLOR;
  const r = residue.toUpperCase();
  if (HYDROPHOBIC.has(r)) return '#7fa8d9';
  if (ACIDIC.has(r)) return '#e08bb0';
  if (BASIC.has(r)) return '#e06666';
  if (POLAR.has(r)) return '#8fce8f';
  return '#f0c869';
}

/** Conservation (most-common share of non-gap) and gap fraction per column. */
export function computeConservation(rows: MsaRow[], length: number): {conservation: Float32Array; gap: Float32Array} {
  const conservation = new Float32Array(length);
  const gap = new Float32Array(length);
  const n = rows.length || 1;
  const counts = new Int32Array(length * 26);
  const nonGap = new Int32Array(length);
  for (const row of rows) {
    const seq = row.seq;
    for (let c = 0; c < length; c++) {
      const ch = seq.charCodeAt(c);
      if (ch >= 65 && ch <= 90) {
        counts[c * 26 + (ch - 65)]++;
        nonGap[c]++;
      } else if (ch >= 97 && ch <= 122) {
        counts[c * 26 + (ch - 97)]++;
        nonGap[c]++;
      }
    }
  }
  for (let c = 0; c < length; c++) {
    gap[c] = 1 - nonGap[c] / n;
    if (!nonGap[c]) {
      conservation[c] = 0;
      continue;
    }
    let best = 0;
    for (let k = 0; k < 26; k++) if (counts[c * 26 + k] > best) best = counts[c * 26 + k];
    conservation[c] = best / nonGap[c];
  }
  return {conservation, gap};
}

// --- row layout -------------------------------------------------------------

export const ROW_H = 15;
export const GROUP_H = 22;

export interface MsaGroup {
  key: string;
  genome: string;
  row: number;
  top: number;
  height: number;
  rowStart: number;
  rowEnd: number;
}

export interface MsaLayout {
  rows: MsaRow[];
  ys: number[];
  groups: MsaGroup[];
  totalHeight: number;
}

/** Stack rows in genome groups (rows must already be sorted by genome row). */
export function layoutMsaRows(rows: MsaRow[], rowH = ROW_H, groupH = GROUP_H): MsaLayout {
  const ys: number[] = [];
  const groups: MsaGroup[] = [];
  let y = 0;
  let i = 0;
  while (i < rows.length) {
    const key = `${rows[i].row}|${rows[i].genome}`;
    const start = i;
    y += groupH;
    while (i < rows.length && `${rows[i].row}|${rows[i].genome}` === key) {
      ys[i] = y;
      y += rowH;
      i++;
    }
    groups.push({key, genome: rows[start].genome, row: rows[start].row, top: y - (i - start) * rowH - groupH, height: groupH + (i - start) * rowH, rowStart: start, rowEnd: i});
  }
  return {rows, ys, groups, totalHeight: y};
}

// --- viewer -----------------------------------------------------------------

const LABEL_W = 210;
const BAR_H = 34;
const CONS_H = 12; // conservation bar: y 2..14
const GAP_Y = 16; // gap strip: y 16..26 (>=8 px)
const GAP_H = 10;
const SB_Y = 28; // overview scrollbar: y 28..34
const SB_H = 6;
export const MSA_LABEL_W = LABEL_W;
export const MSA_BAR_H = BAR_H;

export interface MsaHover {
  index: number;
  col: number;
  residue: string;
  gene: string;
  genome: string;
  ungapped: number;
  conservation: number;
  gap: number;
  bar?: boolean;
}

export class MsaView {
  private data: MsaResult | null = null;
  private layout: MsaLayout | null = null;
  private conservation: Float32Array | null = null;
  private gap: Float32Array | null = null;
  private pxPerCol = 8;
  private scrollCol = 0;
  private scrollY = 0;
  private userZoomed = false;
  private dragging = false;
  private observer: ResizeObserver;

  constructor(
    private canvas: HTMLCanvasElement,
    private onHover: (h: MsaHover | null, clientX: number, clientY: number) => void,
  ) {
    canvas.addEventListener('wheel', (e) => this.wheel(e), {passive: false});
    canvas.addEventListener('mousemove', (e) => this.mousemove(e));
    canvas.addEventListener('mouseleave', () => this.onHover(null, 0, 0));
    canvas.addEventListener('pointerdown', (e) => this.pointerdown(e));
    window.addEventListener('pointermove', (e) => this.pointermove(e));
    window.addEventListener('pointerup', () => (this.dragging = false));
    this.observer = new ResizeObserver(() => this.onResize());
    this.observer.observe(canvas.parentElement ?? canvas);
  }

  setData(data: MsaResult | null): void {
    this.data = data;
    this.layout = data ? layoutMsaRows(data.rows) : null;
    if (data) {
      const {conservation, gap} = computeConservation(data.rows, data.length);
      this.conservation = conservation;
      this.gap = gap;
      this.scrollCol = 0;
      this.scrollY = 0;
      this.userZoomed = false;
      this.pxPerCol = fitPxPerCol(this.gridW(), data.length);
    } else {
      this.conservation = this.gap = null;
    }
    this.draw();
  }

  private gridW(): number {
    return Math.max(1, this.canvas.clientWidth - LABEL_W);
  }

  private viewH(): number {
    return Math.max(1, this.canvas.clientHeight - BAR_H);
  }

  private onResize(): void {
    if (!this.data) {
      this.draw();
      return;
    }
    if (!this.userZoomed) this.pxPerCol = fitPxPerCol(this.gridW(), this.data.length);
    this.scrollCol = clampScrollCol(this.scrollCol, this.gridW(), this.pxPerCol, this.data.length);
    this.scrollY = clampScrollY(this.scrollY, this.layout!.totalHeight, this.viewH());
    this.draw();
  }

  get px(): number {
    return this.pxPerCol;
  }

  zoomBy(factor: number, anchorCol?: number, anchorPx?: number): void {
    if (!this.data) return;
    const width = this.gridW();
    const px = anchorPx ?? width / 2;
    const anchor = anchorCol ?? this.scrollCol + px / this.pxPerCol;
    const next = clampPxPerCol(this.pxPerCol * factor);
    if (next === this.pxPerCol) return;
    this.pxPerCol = next;
    this.scrollCol = zoomedScroll(anchor, px, next, width, this.data.length);
    if (factor !== 1) this.userZoomed = true;
    this.draw();
  }

  private wheel(e: WheelEvent): void {
    if (!this.data) return;
    e.preventDefault();
    this.onHover(null, 0, 0);
    if (e.ctrlKey || e.metaKey) {
      const rect = this.canvas.getBoundingClientRect();
      const anchorPx = Math.max(0, e.clientX - rect.left - LABEL_W);
      const anchorCol = this.scrollCol + anchorPx / this.pxPerCol;
      this.zoomBy(2 ** (-e.deltaY * 0.0015), anchorCol, anchorPx);
      return;
    }
    const horizontal = e.shiftKey || Math.abs(e.deltaX) > Math.abs(e.deltaY);
    if (horizontal) {
      const delta = (e.shiftKey ? e.deltaY : e.deltaX) / this.pxPerCol;
      this.scrollCol = clampScrollCol(this.scrollCol + delta, this.gridW(), this.pxPerCol, this.data.length);
    } else {
      this.scrollY = clampScrollY(this.scrollY + e.deltaY, this.layout!.totalHeight, this.viewH());
    }
    this.draw();
  }

  private pointerdown(e: PointerEvent): void {
    if (!this.data || !this.layout) return;
    const rect = this.canvas.getBoundingClientRect();
    const y = e.clientY - rect.top;
    const x = e.clientX - rect.left;
    if (y >= SB_Y && y < SB_Y + SB_H && x >= LABEL_W) {
      this.dragging = true;
      this.scrollCol = scrollFromBar((x - LABEL_W) / this.gridW(), this.gridW(), this.pxPerCol, this.data.length);
      this.draw();
    }
  }

  private pointermove(e: PointerEvent): void {
    if (!this.dragging || !this.data) return;
    const rect = this.canvas.getBoundingClientRect();
    this.scrollCol = scrollFromBar((e.clientX - rect.left - LABEL_W) / this.gridW(), this.gridW(), this.pxPerCol, this.data.length);
    this.draw();
  }

  private mousemove(e: MouseEvent): void {
    if (!this.data || !this.layout) return;
    const rect = this.canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    if (x < LABEL_W || (y >= SB_Y && y < BAR_H)) {
      this.onHover(null, e.clientX, e.clientY);
      return;
    }
    const col = Math.min(this.data.length - 1, Math.max(0, Math.floor(this.scrollCol + (x - LABEL_W) / this.pxPerCol)));
    if (y < GAP_Y + GAP_H && y >= GAP_Y) {
      this.onHover({index: -1, col, residue: '', gene: '', genome: '', ungapped: -1, conservation: this.conservation![col], gap: this.gap![col], bar: true}, e.clientX, e.clientY);
      return;
    }
    if (y < GAP_Y) {
      this.onHover({index: -1, col, residue: '', gene: '', genome: '', ungapped: -1, conservation: this.conservation![col], gap: this.gap![col], bar: true}, e.clientX, e.clientY);
      return;
    }
    const contentY = y - BAR_H + this.scrollY;
    const index = this.rowAt(contentY);
    if (index < 0) {
      this.onHover(null, e.clientX, e.clientY);
      return;
    }
    const row = this.data.rows[index];
    const residue = row.seq[col] ?? '-';
    let ungapped = 0;
    for (let c = 0; c < col; c++) if (row.seq[c] !== '-') ungapped++;
    this.onHover(
      {index, col, residue, gene: row.gene, genome: row.genome, ungapped: residue === '-' ? -1 : ungapped + 1, conservation: this.conservation![col], gap: this.gap![col]},
      e.clientX,
      e.clientY,
    );
  }

  private rowAt(contentY: number): number {
    const ys = this.layout!.ys;
    for (let i = 0; i < ys.length; i++) {
      if (contentY >= ys[i] && contentY < ys[i] + ROW_H) return i;
    }
    return -1;
  }

  draw(): void {
    const cssW = this.canvas.clientWidth;
    const cssH = this.canvas.clientHeight;
    const ctx = prepareCanvas(this.canvas, cssW, cssH);
    if (!ctx) return;
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, cssW, cssH);
    if (!this.data || !this.layout) {
      this.drawLabels(ctx, cssH);
      this.placeholder(ctx, 'no alignment');
      return;
    }
    const data = this.data;
    const layout = this.layout;
    const gridW = this.gridW();
    const rowAreaH = this.viewH();
    const [c0, c1] = visibleCols(this.scrollCol, gridW, this.pxPerCol, data.length);

    this.drawBar(ctx, c0, c1, gridW);

    // grid (own clip)
    ctx.save();
    ctx.beginPath();
    ctx.rect(LABEL_W, BAR_H, gridW, rowAreaH);
    ctx.clip();
    const letters = showLetters(this.pxPerCol);
    ctx.font = '11px ui-monospace, monospace';
    ctx.textBaseline = 'middle';
    ctx.textAlign = letters ? 'center' : 'left';
    for (let i = 0; i < data.rows.length; i++) {
      const y = BAR_H + layout.ys[i] - this.scrollY;
      if (y + ROW_H < BAR_H || y > cssH) continue;
      const row = data.rows[i];
      for (let c = c0; c < c1; c++) {
        const x = LABEL_W + (c - this.scrollCol) * this.pxPerCol;
        const w = Math.max(1, this.pxPerCol);
        const res = row.seq[c] ?? '-';
        ctx.fillStyle = residueColor(res);
        ctx.fillRect(x, y, w, ROW_H - 1);
        if (letters && res !== '-') {
          ctx.fillStyle = '#1f2937';
          ctx.fillText(res, x + this.pxPerCol / 2, y + ROW_H / 2);
        }
      }
    }
    ctx.textAlign = 'left';
    ctx.restore();
    // labels have their own clip (previously drawn inside the grid clip, so invisible)
    this.drawLabels(ctx, cssH);
  }

  private drawBar(ctx: CanvasRenderingContext2D, c0: number, c1: number, gridW: number): void {
    ctx.fillStyle = '#f8fafc';
    ctx.fillRect(LABEL_W, 0, gridW, BAR_H);
    ctx.strokeStyle = '#e4e7ec';
    ctx.beginPath();
    ctx.moveTo(LABEL_W, BAR_H - 0.5);
    ctx.lineTo(LABEL_W + gridW, BAR_H - 0.5);
    ctx.stroke();
    const cons = this.conservation!;
    const gap = this.gap!;
    for (let c = c0; c < c1; c++) {
      const x = LABEL_W + (c - this.scrollCol) * this.pxPerCol;
      const w = Math.max(1, this.pxPerCol);
      const g = Math.round(60 + 170 * cons[c]);
      ctx.fillStyle = `rgb(${Math.round(g * 0.7)},${g},${Math.round(g * 0.7)})`;
      ctx.fillRect(x, 2, w, CONS_H);
      if (gap[c] > 0.005) {
        ctx.fillStyle = `rgba(146,64,14,${0.25 + 0.65 * Math.min(1, gap[c])})`;
        ctx.fillRect(x, GAP_Y, w, GAP_H);
      }
    }
    // overview scrollbar
    ctx.fillStyle = '#eef2f6';
    ctx.fillRect(LABEL_W, SB_Y, gridW, SB_H);
    const viewportFrac = gridW / this.pxPerCol / Math.max(1, this.data!.length);
    const handleW = Math.max(8, gridW * Math.min(1, viewportFrac));
    const frac = this.data!.length > gridW / this.pxPerCol ? this.scrollCol / (this.data!.length - gridW / this.pxPerCol) : 0;
    ctx.fillStyle = '#94a3b8';
    ctx.fillRect(LABEL_W + Math.min(gridW - handleW, frac * (gridW - handleW)), SB_Y, handleW, SB_H);
    ctx.fillStyle = '#475467';
    ctx.font = '10px system-ui, sans-serif';
    ctx.textBaseline = 'top';
    ctx.fillText('conservation / gap', 6, 4);
  }

  private drawLabels(ctx: CanvasRenderingContext2D, cssH: number): void {
    const layout = this.layout;
    ctx.fillStyle = '#f8fafc';
    ctx.fillRect(0, 0, LABEL_W, cssH);
    ctx.strokeStyle = '#e4e7ec';
    ctx.beginPath();
    ctx.moveTo(LABEL_W - 0.5, 0);
    ctx.lineTo(LABEL_W - 0.5, cssH);
    ctx.stroke();
    if (!layout) return;
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, BAR_H, LABEL_W, Math.max(1, cssH - BAR_H));
    ctx.clip();
    ctx.textBaseline = 'middle';
    for (const g of layout.groups) {
      const top = BAR_H + g.top - this.scrollY;
      if (top > cssH || top + g.height < BAR_H) continue;
      ctx.fillStyle = '#111827';
      ctx.font = '11px system-ui, sans-serif';
      ctx.fillText(g.genome || 'unresolved', 6, top + 11);
      for (let i = g.rowStart; i < g.rowEnd; i++) {
        const y = BAR_H + layout.ys[i] - this.scrollY + ROW_H / 2;
        ctx.fillStyle = '#475467';
        ctx.font = '10px system-ui, sans-serif';
        ctx.fillText(this.data!.rows[i].gene, 14, y);
      }
      ctx.strokeStyle = '#eef2f6';
      ctx.beginPath();
      ctx.moveTo(0, top + g.height - 0.5);
      ctx.lineTo(LABEL_W, top + g.height - 0.5);
      ctx.stroke();
    }
    ctx.restore();
  }

  private placeholder(ctx: CanvasRenderingContext2D, text: string): void {
    ctx.fillStyle = '#98a2b3';
    ctx.font = '12px system-ui, sans-serif';
    ctx.textBaseline = 'top';
    ctx.fillText(text, 12, 12);
  }
}
