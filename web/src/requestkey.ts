// Request dedupe keys for the matrix/tracks window. Extracted so the
// "reset on cleared state" rule is testable.

import type {Window} from './transform';

export interface RequestKeys {
  matrix: string;
  tracks: string;
}

/** Keys describing a matrix window (order/value/window) and its track layers. */
export function requestKeys(
  order: string,
  value: string,
  win: Window,
  tracks: string[],
  annos: string[],
): RequestKeys {
  const matrix = `${order}|${value}|${win.og0}|${win.og1}|${win.g0}|${win.g1}|${win.w}|${win.h}`;
  return {matrix, tracks: `${matrix}|${tracks.join(',')}|${annos.join(',')}`};
}

export interface FetchDecision {
  matrix: boolean;
  tracks: boolean;
}

/**
 * Decide what to fetch. `hasMatrix` is false when the displayed matrix state
 * was cleared (e.g. an order change), which forces a matrix fetch even if the
 * key happens to match a stale one.
 */
export function needsFetch(
  keys: RequestKeys,
  lastMatrixKey: string | null,
  lastTracksKey: string | null,
  hasMatrix: boolean,
): FetchDecision {
  return {
    matrix: !hasMatrix || keys.matrix !== lastMatrixKey,
    tracks: keys.tracks !== lastTracksKey,
  };
}
