// Alignment column math and residue colours (pure; the virtualised renderer
// lives in alignment.ts).

import type {MsaRow} from './api';

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
