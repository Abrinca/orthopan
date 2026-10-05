import {describe, expect, test} from 'bun:test';
import {needsFetch, requestKeys} from '../src/requestkey';

const win = {og0: 0, og1: 10, g0: 0, g1: 5, w: 10, h: 5};
const tree = requestKeys('tree', 'presence', win, [], []);
const freq = requestKeys('frequency', 'presence', win, [], []);

describe('request dedupe keys', () => {
  test('unchanged state fetches nothing', () => {
    expect(needsFetch(tree, tree.matrix, tree.tracks, true)).toEqual({matrix: false, tracks: false});
  });

  test('a different order changes the matrix key', () => {
    expect(freq.matrix).not.toBe(tree.matrix);
    expect(needsFetch(freq, tree.matrix, tree.tracks, true).matrix).toBe(true);
  });

  test('toggling a track only changes the tracks key', () => {
    const withTrack = requestKeys('tree', 'presence', win, ['copy_number_mean'], []);
    expect(withTrack.matrix).toBe(tree.matrix);
    expect(withTrack.tracks).not.toBe(tree.tracks);
    const d = needsFetch(withTrack, tree.matrix, tree.tracks, true);
    expect(d).toEqual({matrix: false, tracks: true});
  });

  test('cleared displayed state forces a refetch even when the key matches', () => {
    // A -> B -> A: B was fetched, then setOrder() cleared dataOrder and the keys,
    // so returning to A must fetch again rather than dedupe against stale A.
    expect(needsFetch(freq, tree.matrix, tree.tracks, true).matrix).toBe(true);
    const reset = needsFetch(tree, null, null, false);
    expect(reset).toEqual({matrix: true, tracks: true});
  });

  test('hasMatrix=false forces the matrix fetch but may skip tracks', () => {
    expect(needsFetch(tree, tree.matrix, tree.tracks, false)).toEqual({matrix: true, tracks: false});
  });
});
