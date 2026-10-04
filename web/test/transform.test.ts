import {describe, expect, test} from 'bun:test';
import {
  bufferAt,
  capCellCounts,
  cellRect,
  clampView,
  colRange,
  computeWindow,
  fitView,
  makeTransform,
  rowRange,
  type View,
} from '../src/transform';

const close = (a: number, b: number, digits = 6) => expect(a).toBeCloseTo(b, digits);

describe('transform round-trip', () => {
  const v: View = {cx: 12.3, cy: 4.7, zx: 1.3, zy: -0.4};
  const pw = 1010;
  const ph = 792;
  const tr = makeTransform(v, pw, ph);

  test('xOf/ogAt invert', () => {
    for (const px of [0, 123.4, 505, 1009.9]) close(tr.xOf(tr.ogAt(px)), px);
  });

  test('yOf/rowAt invert', () => {
    for (const py of [0, 300.7, 790]) close(tr.yOf(tr.rowAt(py)), py);
  });

  test('scale equals 2^zoom', () => {
    close(tr.sx, 2 ** v.zx);
    close(tr.sy, 2 ** v.zy);
  });
});

describe('bitmap <-> label alignment', () => {
  test('row i centre matches between bitmap cell and label at full detail', () => {
    const nOG = 5657;
    const nG = 24;
    const pw = 1010;
    const ph = 792;
    const v = fitView(nOG, nG, pw, ph);
    const win = computeWindow(v, nOG, nG, pw, ph);
    expect(win.g0).toBe(0);
    expect(win.g1).toBe(nG);
    expect(win.h).toBe(nG);
    const tr = makeTransform(v, pw, ph);
    for (const r of [0, 5, 17, 23]) {
      const [, y, , h] = cellRect(win, tr, 0, r);
      // the bitmap cell centre and the label (yOf(r+0.5)) must coincide
      close(y + h / 2, tr.yOf(r + 0.5));
    }
  });

  test('column boundaries align with integer OG positions (fractional window)', () => {
    const nOG = 1000;
    const nG = 50;
    const pw = 1000;
    const ph = 500;
    const v: View = {cx: 250.37, cy: 10.42, zx: 0, zy: 0}; // 1 px per unit
    const win = computeWindow(v, nOG, nG, pw, ph);
    expect(Number.isInteger(win.og0)).toBe(true);
    const tr = makeTransform(v, pw, ph);
    // Just past the boundary between OG positions P and P+1.
    const P = win.og0 + 10;
    const hit = bufferAt(win, tr, tr.xOf(P + 1) + 0.01, tr.yOf(win.g0 + 3.5), nG, nOG)!;
    expect(hit.og).toBe(P + 1);
    expect(hit.col).toBe(hit.og - win.og0); // full detail on X
    const [x0] = cellRect(win, tr, hit.col, hit.row);
    close(x0, tr.xOf(P + 1));
  });

  test('buffer ranges partition the window exactly', () => {
    const win = {og0: 100, og1: 250, g0: 5, g1: 40, w: 75, h: 20};
    for (let i = 0; i < win.w; i++) {
      const [a, b] = colRange(win, i);
      if (i === 0) expect(a).toBe(win.og0);
      if (i === win.w - 1) expect(b).toBe(win.og1);
      const next = colRange(win, i + 1);
      if (i < win.w - 1) expect(b).toBe(next[0]);
    }
    for (let j = 0; j < win.h; j++) {
      const [a, b] = rowRange(win, j);
      if (j === 0) expect(a).toBe(win.g0);
      if (j === win.h - 1) expect(b).toBe(win.g1);
    }
  });

  test('bufferAt returns null outside the fetched half-open window (each edge)', () => {
    const win = {og0: 10, og1: 20, g0: 5, g1: 15, w: 10, h: 10};
    const tr = makeTransform({cx: 15, cy: 10, zx: 0, zy: 0}, 100, 100);
    const at = (x: number, y: number) => bufferAt(win, tr, x, y, 100, 100);
    // inside
    expect(at(50, 50)).not.toBeNull();
    // exactly on the inclusive low edges is inside
    expect(at(tr.xOf(10), 50)).not.toBeNull();
    expect(at(50, tr.yOf(5))).not.toBeNull();
    // left/right edges (half-open, right exclusive)
    expect(at(tr.xOf(10) - 0.5, 50)).toBeNull();
    expect(at(tr.xOf(20), 50)).toBeNull();
    expect(at(tr.xOf(20) - 0.5, 50)).not.toBeNull();
    // top/bottom edges
    expect(at(50, tr.yOf(5) - 0.5)).toBeNull();
    expect(at(50, tr.yOf(15))).toBeNull();
    expect(at(50, tr.yOf(15) - 0.5)).not.toBeNull();
    // corners
    expect(at(tr.xOf(10) - 1, tr.yOf(5) - 1)).toBeNull();
    expect(at(tr.xOf(20), tr.yOf(15))).toBeNull();
  });

  test('cellRect uses continuous equal-width texel boundaries', () => {
    // L=150 over w=40 -> 3.75 world units/texel; M=35 over h=20 -> 1.75
    const win = {og0: 100, og1: 250, g0: 5, g1: 40, w: 40, h: 20};
    const tr = makeTransform({cx: 0, cy: 0, zx: 0, zy: 0}, 1000, 1000);
    const r0 = cellRect(win, tr, 0, 0);
    const r1 = cellRect(win, tr, 1, 0);
    close(r0[2], tr.xOf(win.og0 + 3.75) - tr.xOf(win.og0)); // width = L/w
    close(r1[0], r0[0] + r0[2]); // contiguous
    close(r1[2], r0[2]); // equal widths
    close(r0[3], tr.yOf(win.g0 + 1.75) - tr.yOf(win.g0)); // height = M/h
    // the outline follows the bitmap texel edge, not the floor-rounded range
    expect(colRange(win, 1)[0]).toBe(win.og0 + 3);
    close(r1[0], tr.xOf(win.og0 + 3.75));
  });
});

describe('resize / fit', () => {
  test('fitView window covers the full range at several sizes', () => {
    const nOG = 50515;
    const nG = 2193;
    for (const [pw, ph] of [
      [1010, 792],
      [1800, 1100],
      [900, 500],
    ]) {
      const v = fitView(nOG, nG, pw, ph);
      const win = computeWindow(v, nOG, nG, pw, ph);
      expect(win.og0).toBe(0);
      expect(win.og1).toBe(nOG);
      expect(win.g0).toBe(0);
      expect(win.g1).toBe(nG);
    }
  });

  test('clampView is idempotent', () => {
    const v = fitView(100, 10, 1000, 500);
    const c = clampView(v, 100, 10, 1000, 500);
    expect(clampView(c, 100, 10, 1000, 500)).toEqual(c);
  });
});

describe('request cell cap', () => {
  test('keeps w*h <= 4M and preserves aspect ratio (3000x2000)', () => {
    const r = capCellCounts(3000, 2000);
    expect(r.w * r.h).toBeLessThanOrEqual(4_000_000);
    close(r.w / r.h, 3000 / 2000, 1);
  });

  test('no-op below the cap', () => {
    expect(capCellCounts(100, 50)).toEqual({w: 100, h: 50});
  });

  test('4096x4096 is reduced', () => {
    const r = capCellCounts(4096, 4096);
    expect(r.w * r.h).toBeLessThanOrEqual(4_000_000);
    close(r.w / r.h, 1, 1);
  });
});
