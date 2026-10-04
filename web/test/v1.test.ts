import {describe, expect, test} from 'bun:test';
import {CATEGORY_COLORS, CATEGORY_GREY, annotationStripColor, categoryColor} from '../src/annotations';
import {copiesBucket, copiesMean, decodeFloat32LE, decodeUint16LE, formatCopiesExact, formatCopiesMean, COPIES_SATURATED} from '../src/decode';
import {hitViewport} from '../src/search';
import type {View} from '../src/transform';

const close = (a: number, b: number, digits = 6) => expect(a).toBeCloseTo(b, digits);

describe('little-endian decoders', () => {
  test('decodeUint16LE', () => {
    const buf = new ArrayBuffer(6);
    const dv = new DataView(buf);
    dv.setUint16(0, 1, true);
    dv.setUint16(2, 258, true);
    dv.setUint16(4, 65535, true);
    expect([...decodeUint16LE(buf)]).toEqual([1, 258, 65535]);
  });

  test('decodeUint16LE ignores a trailing odd byte', () => {
    const buf = new ArrayBuffer(3);
    new DataView(buf).setUint16(0, 7, true);
    expect([...decodeUint16LE(buf)]).toEqual([7]);
  });

  test('decodeFloat32LE', () => {
    const buf = new ArrayBuffer(8);
    const dv = new DataView(buf);
    dv.setFloat32(0, 1.5, true);
    dv.setFloat32(4, -2.25, true);
    const out = decodeFloat32LE(buf);
    expect(out.length).toBe(2);
    close(out[0], 1.5);
    close(out[1], -2.25);
  });
});

describe('copy-number values', () => {
  test('raw byte is 10x the mean', () => {
    expect(copiesMean(0)).toBe(0);
    expect(copiesMean(10)).toBe(1);
    close(copiesMean(25), 2.5);
  });

  test('discrete buckets: 0 absent, 1..5, 5+ clamps', () => {
    expect(copiesBucket(0)).toBe(0);
    expect(copiesBucket(3)).toBe(1); // 0.3 -> clamped to 1
    expect(copiesBucket(10)).toBe(1);
    expect(copiesBucket(25)).toBe(3); // round(2.5)
    expect(copiesBucket(49)).toBe(5);
    expect(copiesBucket(90)).toBe(5); // 9 -> 5+
    expect(copiesBucket(255)).toBe(5);
  });

  test('saturated copy bytes are shown as a lower bound', () => {
    expect(COPIES_SATURATED).toBe(255);
    expect(formatCopiesMean(255)).toBe('≥25.5');
    expect(formatCopiesMean(0)).toBe('0');
    expect(formatCopiesMean(25)).toBe('2.50');
    expect(formatCopiesExact(255)).toBe('≥25.5');
    expect(formatCopiesExact(20)).toBe('2');
    expect(formatCopiesExact(0)).toBe('0');
  });
});

describe('categorical palette', () => {
  test('top-N categories are distinct colours', () => {
    const seen = new Set<string>();
    for (let i = 0; i < CATEGORY_COLORS.length; i++) {
      const c = categoryColor(i);
      expect(c).not.toBe(CATEGORY_GREY);
      expect(seen.has(c)).toBe(false);
      seen.add(c);
    }
  });

  test('categories beyond the palette are grey', () => {
    expect(categoryColor(-1)).toBe(CATEGORY_GREY);
    expect(categoryColor(CATEGORY_COLORS.length)).toBe(CATEGORY_GREY);
    expect(categoryColor(1000)).toBe(CATEGORY_GREY);
  });

  test('annotation values map: 0 missing, else 1 + category index', () => {
    expect(annotationStripColor(0)).not.toBe(categoryColor(0));
    expect(annotationStripColor(1)).toBe(categoryColor(0));
    expect(annotationStripColor(2)).toBe(categoryColor(1));
  });
});

describe('search-hit viewport', () => {
  const base: View = {cx: 100, cy: 10, zx: 1, zy: 3};
  test('centres the hit and sets a label-visible column width', () => {
    const v = hitViewport(1234, base);
    close(v.cx, 1234.5);
    close(2 ** v.zx, 40);
    // other axes untouched
    close(v.cy, base.cy);
    close(v.zy, base.zy);
  });

  test('respects a custom target px/OG', () => {
    const v = hitViewport(7, base, 24);
    close(2 ** v.zx, 24);
    close(v.cx, 7.5);
  });
});
