import {describe, expect, test} from 'bun:test';
import {binRange, makeTransform} from '../src/transform';
import {Tracks, type CategoryStrip} from '../src/tracks';

const strip: CategoryStrip = {
  kind: 'category',
  name: 'cog',
  categories: ['A', 'B', 'C'],
  result: {name: 'cog', data: new Uint16Array([1, 2, 3, 0, 3, 2, 1]), bins: 7},
  world: {og0: 10, og1: 25},
};

// sx = 1 px/unit, so ogAt(px) = px - 50; the strip spans px [60, 75].
const tr = makeTransform({cx: 0, cy: 0, zx: 0, zy: 0}, 100, 100);
const tracks = new Tracks({} as HTMLCanvasElement);

describe('strip hit-testing', () => {
  test('returns null outside the fetched orthogroup span', () => {
    expect(tracks.hit(tr, [strip], 55, 10, 100)).toBeNull(); // op 5 < 10
    expect(tracks.hit(tr, [strip], 75, 10, 100)).toBeNull(); // op 25 == og1
    expect(tracks.hit(tr, [strip], 80, 10, 100)).toBeNull(); // stale window edge
  });

  test('returns null outside the canvas', () => {
    expect(tracks.hit(tr, [strip], 60, -1, 100)).toBeNull();
    expect(tracks.hit(tr, [strip], 60, 100, 100)).toBeNull();
  });

  test('maps a pointer to the right bin', () => {
    const first = tracks.hit(tr, [strip], 60, 10, 100);
    expect(first?.index).toBe(0);
    expect(first?.value).toBe(1);
    const last = tracks.hit(tr, [strip], 74.5, 10, 100);
    expect(last?.index).toBe(6);
    expect(last?.value).toBe(1);
  });
});

describe('binRange', () => {
  test('partitions the span with floor-rounded boundaries (non-even window)', () => {
    const bins = 7;
    const start = 10;
    const end = 25;
    const ranges = Array.from({length: bins}, (_, i) => binRange(start, end, bins, i));
    expect(ranges[0]).toEqual([10, 12]);
    expect(ranges[bins - 1]).toEqual([22, 25]);
    for (let i = 0; i < bins; i++) {
      expect(ranges[i][0]).toBeLessThanOrEqual(ranges[i][1]);
      if (i > 0) expect(ranges[i][0]).toBe(ranges[i - 1][1]); // contiguous
    }
  });

  test('a full-detail window maps one bin to one OG', () => {
    expect(binRange(100, 150, 50, 7)).toEqual([107, 108]);
  });
});
