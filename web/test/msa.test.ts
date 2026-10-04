import {describe, expect, test} from 'bun:test';
import type {MsaRow} from '../src/api';
import {
  clampPxPerCol,
  clampScrollCol,
  clampScrollY,
  computeConservation,
  fitPxPerCol,
  GAP_COLOR,
  layoutMsaRows,
  LETTER_PX,
  MAX_PX_PER_COL,
  MIN_PX_PER_COL,
  residueColor,
  scrollFromBar,
  showLetters,
  visibleCols,
  zoomedScroll,
} from '../src/msa';

describe('MSA column math', () => {
  test('visible column range', () => {
    expect(visibleCols(0, 100, 10, 50)).toEqual([0, 10]);
    expect(visibleCols(5.5, 100, 10, 50)).toEqual([5, 16]);
    // clamped to length at the right edge
    expect(visibleCols(45, 100, 10, 50)).toEqual([45, 50]);
  });

  test('zoom clamps', () => {
    expect(clampPxPerCol(0)).toBe(MIN_PX_PER_COL);
    expect(clampPxPerCol(1000)).toBe(MAX_PX_PER_COL);
    expect(clampPxPerCol(4)).toBe(4);
  });

  test('letter threshold', () => {
    expect(showLetters(LETTER_PX - 0.01)).toBe(false);
    expect(showLetters(LETTER_PX)).toBe(true);
  });

  test('scroll clamps inside [0, length - viewport]', () => {
    expect(clampScrollCol(-5, 100, 10, 50)).toBe(0);
    expect(clampScrollCol(100, 100, 10, 50)).toBe(40);
    expect(clampScrollCol(10, 1000, 10, 50)).toBe(0); // viewport wider than content
  });

  test('vertical scroll clamps to [0, contentH - viewH]', () => {
    expect(clampScrollY(-3, 1000, 500)).toBe(0);
    expect(clampScrollY(800, 1000, 500)).toBe(500);
    expect(clampScrollY(100, 400, 600)).toBe(0); // content shorter than viewport
  });

  test('off-centre pointer-anchored zoom keeps the anchor pixel offset', () => {
    const width = 1000;
    const length = 200;
    const pxPerCol = 5;
    const scrollCol = 40;
    const anchorPx = 300; // pointer offset from the grid left edge
    const anchorCol = scrollCol + anchorPx / pxPerCol; // 100
    const newPx = 10;
    const next = zoomedScroll(anchorCol, anchorPx, newPx, width, length);
    expect(next).toBe(70); // anchorCol - anchorPx/newPx
    expect((anchorCol - next) * newPx).toBeCloseTo(anchorPx);
  });

  test('zero-delta zoom leaves the scroll unchanged', () => {
    const width = 1000;
    const length = 400;
    const pxPerCol = 10;
    const scrollCol = 40;
    const anchorPx = 300;
    const anchorCol = scrollCol + anchorPx / pxPerCol;
    expect(zoomedScroll(anchorCol, anchorPx, pxPerCol, width, length)).toBe(scrollCol);
  });

  test('fit-to-width caps at 12 px/col and zoom clamps apply', () => {
    expect(fitPxPerCol(1000, 200)).toBe(5);
    expect(fitPxPerCol(12000, 200)).toBe(12);
    expect(fitPxPerCol(100, 1000)).toBe(MIN_PX_PER_COL);
  });

  test('overview strip maps a fraction to a centred viewport', () => {
    expect(scrollFromBar(0.5, 1000, 10, 200)).toBe(50);
    expect(scrollFromBar(-1, 1000, 10, 200)).toBe(0);
    expect(scrollFromBar(2, 1000, 10, 200)).toBe(100);
  });
});

describe('MSA conservation', () => {
  const rows: MsaRow[] = [
    {genome: 'g0', row: 0, gene: 'a', raw: 'a', seq: 'AA-'},
    {genome: 'g0', row: 0, gene: 'b', raw: 'b', seq: 'AB-'},
    {genome: 'g1', row: 1, gene: 'c', raw: 'c', seq: 'A--'},
  ];
  test('per-column conservation and gap fraction', () => {
    const {conservation, gap} = computeConservation(rows, 3);
    expect(conservation[0]).toBe(1); // all A
    expect(gap[0]).toBe(0);
    expect(conservation[1]).toBe(0.5); // A,B (one gap) -> 1 of 2
    expect(gap[1]).toBeCloseTo(1 / 3);
    expect(conservation[2]).toBe(0); // all gap
    expect(gap[2]).toBe(1);
  });
});

describe('MSA residue colours', () => {
  test('Clustal-like classes and gaps', () => {
    const hydro = residueColor('A');
    const acid = residueColor('D');
    const basic = residueColor('K');
    const polar = residueColor('S');
    const other = residueColor('X');
    expect(new Set([hydro, acid, basic, polar, other]).size).toBe(5);
    expect(residueColor('-')).toBe(GAP_COLOR);
    expect(residueColor('a')).toBe(hydro); // case-insensitive
  });
});

describe('MSA row layout', () => {
  const rows: MsaRow[] = [
    {genome: 'g0', row: 0, gene: 'a', raw: 'a', seq: 'AA'},
    {genome: 'g0', row: 0, gene: 'b', raw: 'b', seq: 'AB'},
    {genome: 'g1', row: 1, gene: 'c', raw: 'c', seq: 'AC'},
  ];
  test('groups rows by genome and stacks group headers', () => {
    const l = layoutMsaRows(rows, 15, 22);
    expect(l.groups.length).toBe(2);
    expect(l.groups[0].rowStart).toBe(0);
    expect(l.groups[0].rowEnd).toBe(2);
    expect(l.groups[1].rowStart).toBe(2);
    expect(l.groups[1].rowEnd).toBe(3);
    expect(l.ys[0]).toBe(22); // first row below its group header
    expect(l.ys[1]).toBe(37);
    expect(l.ys[2]).toBe(22 + 2 * 15 + 22); // second group header
    expect(l.totalHeight).toBe(22 + 2 * 15 + 22 + 15);
  });
});
