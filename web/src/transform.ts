// Single source of truth for the world <-> pixel mapping shared by the heatmap,
// tree, labels, tracks and pointer interaction. World X is an orthogroup
// position, world Y is a genome row; a genome row r occupies world y [r, r+1]
// and its centre is r + 0.5. The heatmap bitmap is placed by deck.gl at its
// fetched window rectangle, which is rendered through the same transform.

export interface View {
  cx: number;
  cy: number;
  zx: number;
  zy: number;
}

/** A half-open fetched window plus the requested output-cell counts. */
export interface Window {
  og0: number;
  og1: number;
  g0: number;
  g1: number;
  w: number;
  h: number;
}

const MAX_PIXELS_PER_UNIT = 40;
const MIN_PIXELS_PER_UNIT = 0.01;
// Backend per-request limits (src/orthopan/server.py MAX_CELLS / MAX_CELLS_TOTAL).
export const MAX_CELLS = 4096;
export const MAX_CELLS_TOTAL = 4_000_000;

export const clamp = (v: number, lo: number, hi: number) => Math.min(hi, Math.max(lo, v));

export function clampView(v: View, nOG: number, nG: number, pw: number, ph: number): View {
  const fitX = Math.log2(pw / Math.max(1, nOG));
  const fitY = Math.log2(ph / Math.max(1, nG));
  const zx = clamp(v.zx, Math.max(fitX, Math.log2(MIN_PIXELS_PER_UNIT)), Math.log2(MAX_PIXELS_PER_UNIT));
  const zy = clamp(v.zy, Math.max(fitY, Math.log2(MIN_PIXELS_PER_UNIT)), Math.log2(MAX_PIXELS_PER_UNIT));
  const halfW = pw / 2 / 2 ** zx;
  const halfH = ph / 2 / 2 ** zy;
  const cx = 2 * halfW >= nOG ? nOG / 2 : clamp(v.cx, halfW, nOG - halfW);
  const cy = 2 * halfH >= nG ? nG / 2 : clamp(v.cy, halfH, nG - halfH);
  return {cx, cy, zx, zy};
}

export function fitView(nOG: number, nG: number, pw: number, ph: number): View {
  return {
    cx: nOG / 2,
    cy: nG / 2,
    zx: Math.log2(pw / Math.max(1, nOG)),
    zy: Math.log2(ph / Math.max(1, nG)),
  };
}

/** Reduce w/h, preserving aspect ratio, so that w*h <= MAX_CELLS_TOTAL. */
export function capCellCounts(w: number, h: number): {w: number; h: number} {
  w = Math.max(1, Math.floor(w));
  h = Math.max(1, Math.floor(h));
  if (w * h <= MAX_CELLS_TOTAL) return {w, h};
  const k = Math.sqrt(MAX_CELLS_TOTAL / (w * h));
  let nw = Math.max(1, Math.floor(w * k));
  let nh = Math.max(1, Math.floor(h * k));
  while (nw * nh > MAX_CELLS_TOTAL) {
    if (nw >= nh) nw--;
    else nh--;
  }
  return {w: nw, h: nh};
}

/** Visible half-open window plus the requested output-cell counts. */
export function computeWindow(v: View, nOG: number, nG: number, pw: number, ph: number): Window {
  const halfW = pw / 2 / 2 ** v.zx;
  const halfH = ph / 2 / 2 ** v.zy;
  const og0 = clamp(Math.floor(v.cx - halfW), 0, Math.max(0, nOG - 1));
  const og1 = clamp(Math.ceil(v.cx + halfW), og0 + 1, nOG);
  const g0 = clamp(Math.floor(v.cy - halfH), 0, Math.max(0, nG - 1));
  const g1 = clamp(Math.ceil(v.cy + halfH), g0 + 1, nG);
  const w = clamp(Math.round((og1 - og0) * 2 ** v.zx), 1, Math.min(og1 - og0, pw, MAX_CELLS));
  const h = clamp(Math.round((g1 - g0) * 2 ** v.zy), 1, Math.min(g1 - g0, ph, MAX_CELLS));
  const capped = capCellCounts(w, h);
  return {og0, og1, g0, g1, ...capped};
}

/** Backing-store size for a canvas whose CSS box is cssPw x cssPh at the given dpr. */
export function backingSize(cssPw: number, cssPh: number, dpr: number): {width: number; height: number} {
  return {width: Math.max(1, Math.round(cssPw * dpr)), height: Math.max(1, Math.round(cssPh * dpr))};
}

/**
 * Size a canvas for DPR-correct drawing: backing store = CSS size x dpr, one
 * ctx.setTransform(dpr,...) so all drawing is in CSS px. The canvas's CSS box
 * must already be cssPw x cssPh (callers rely on CSS width/height:100%).
 */
export function prepareCanvas(
  canvas: HTMLCanvasElement,
  cssPw: number,
  cssPh: number,
): CanvasRenderingContext2D | null {
  const dpr = window.devicePixelRatio || 1;
  const {width, height} = backingSize(cssPw, cssPh, dpr);
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }
  const ctx = canvas.getContext('2d');
  if (!ctx) return null;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssPw, cssPh);
  return ctx;
}

export interface Transform {
  sx: number;
  sy: number;
  xOf(x: number): number;
  yOf(r: number): number;
  ogAt(px: number): number;
  rowAt(py: number): number;
}

export function makeTransform(v: View, pw: number, ph: number): Transform {
  const sx = 2 ** v.zx;
  const sy = 2 ** v.zy;
  return {
    sx,
    sy,
    xOf: (x) => (x - v.cx) * sx + pw / 2,
    yOf: (r) => (r - v.cy) * sy + ph / 2,
    ogAt: (px) => (px - pw / 2) / sx + v.cx,
    rowAt: (py) => (py - ph / 2) / sy + v.cy,
  };
}

/** Buffer column j covers OG positions [og0 + floor(j*L/w), og0 + floor((j+1)*L/w)). */
export function colRange(win: Window, j: number): [number, number] {
  const L = win.og1 - win.og0;
  return [win.og0 + Math.floor((j * L) / win.w), win.og0 + Math.floor(((j + 1) * L) / win.w)];
}

/** Buffer row i covers genome rows [g0 + floor(i*M/h), g0 + floor((i+1)*M/h)). */
export function rowRange(win: Window, i: number): [number, number] {
  const M = win.g1 - win.g0;
  return [win.g0 + Math.floor((i * M) / win.h), win.g0 + Math.floor(((i + 1) * M) / win.h)];
}

/** Floor-rounded range of strip bin `i` over the orthogroup span [start, end). */
export function binRange(start: number, end: number, bins: number, i: number): [number, number] {
  return [start + Math.floor((i * (end - start)) / bins), start + Math.floor(((i + 1) * (end - start)) / bins)];
}

/**
 * Pixel rectangle of one buffer cell (used for the hover outline). Uses the
 * bitmap's continuous, equal-width texel boundaries (not the floor-rounded
 * ranges used for tooltip text).
 */
export function cellRect(win: Window, tr: Transform, col: number, row: number): [number, number, number, number] {
  const L = win.og1 - win.og0;
  const M = win.g1 - win.g0;
  const x0 = tr.xOf(win.og0 + (col * L) / win.w);
  const x1 = tr.xOf(win.og0 + ((col + 1) * L) / win.w);
  const y0 = tr.yOf(win.g0 + (row * M) / win.h);
  const y1 = tr.yOf(win.g0 + ((row + 1) * M) / win.h);
  return [x0, y0, Math.max(1, x1 - x0), Math.max(1, y1 - y0)];
}

export interface CellHit {
  row: number;
  col: number;
  genome: number;
  og: number;
  /** true when the buffer cell aggregates more than one genome row / orthogroup */
  binnedX: boolean;
  binnedY: boolean;
}

/**
 * Pointer (matrix-local px/py) -> buffer cell and the single genome/OG under it.
 * Returns null when the pointer maps outside the fetched half-open window, so
 * stale data (e.g. during a blocked fetch or right after an order switch) is
 * not interactive.
 */
export function bufferAt(
  win: Window,
  tr: Transform,
  px: number,
  py: number,
  nG: number,
  nOG: number,
): CellHit | null {
  const gr = tr.rowAt(py);
  const op = tr.ogAt(px);
  if (op < win.og0 || op >= win.og1 || gr < win.g0 || gr >= win.g1) return null;
  const L = win.og1 - win.og0;
  const M = win.g1 - win.g0;
  const row = Math.min(win.h - 1, Math.floor(((gr - win.g0) * win.h) / M));
  const col = Math.min(win.w - 1, Math.floor(((op - win.og0) * win.w) / L));
  return {
    row,
    col,
    genome: clamp(Math.floor(gr), 0, nG - 1),
    og: clamp(Math.floor(op), 0, nOG - 1),
    binnedX: win.w < L,
    binnedY: win.h < M,
  };
}
