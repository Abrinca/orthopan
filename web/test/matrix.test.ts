import {describe, expect, test} from 'bun:test';
import {
  FLAT_TILE,
  isShaderError,
  TILE_ABSENT_DARKEN,
  TILE_FULL_PX,
  TILE_GAP_MAX,
  TILE_GAP_MIN,
  TILE_MIN_PX,
  TILE_RADIUS_FRAC,
  tileFade,
  tileGap,
  tileParams,
  tileRadius,
} from '../src/matrix';

describe('tile activation fade', () => {
  test('flat below the threshold, fully tiled above it', () => {
    expect(tileFade(0)).toBe(0);
    expect(tileFade(TILE_MIN_PX)).toBe(0);
    expect(tileFade(TILE_MIN_PX + 0.001)).toBeGreaterThan(0);
    expect(tileFade(TILE_FULL_PX)).toBe(1);
    expect(tileFade(100)).toBe(1);
  });

  test('smoothstep midpoint', () => {
    expect(tileFade((TILE_MIN_PX + TILE_FULL_PX) / 2)).toBeCloseTo(0.5, 6);
  });

  test('monotonic', () => {
    let prev = -1;
    for (let px = 0; px <= 14; px += 0.5) {
      const f = tileFade(px);
      expect(f).toBeGreaterThanOrEqual(prev);
      prev = f;
    }
  });
});

describe('tile gap and radius', () => {
  test('gap is ~12% clamped to [1, 3]', () => {
    expect(tileGap(2)).toBe(TILE_GAP_MIN);
    expect(tileGap(10)).toBeCloseTo(1.2, 6);
    expect(tileGap(40)).toBe(TILE_GAP_MAX);
  });

  test('radius is ~25% of the smaller dimension', () => {
    expect(tileRadius(20, 40, 2, 2)).toBeCloseTo(20 * TILE_RADIUS_FRAC, 6);
  });

  test('radius never collapses the tile into the gap', () => {
    const w = 8;
    const h = 8;
    const g = tileGap(8);
    const r = tileRadius(w, h, g, g);
    expect(r).toBeLessThanOrEqual(w / 2 - g / 2);
  });
});

describe('tile params', () => {
  test('binned (not full detail) is flat regardless of cell size', () => {
    expect(tileParams(40, 40, [10, 10], false, 1)).toEqual({...FLAT_TILE, cellCount: [10, 10]});
  });

  test('small cells are flat', () => {
    const p = tileParams(5, 5, [100, 100], true, 1);
    expect(p.mix).toBe(0);
    expect(p.gap).toEqual([0, 0]);
  });

  test('large cells are tiled and scaled by dpr', () => {
    const p = tileParams(16, 16, [100, 100], true, 2);
    expect(p.mix).toBe(1);
    expect(p.cellCount).toEqual([100, 100]);
    expect(p.cellSize).toEqual([32, 32]);
    expect(p.gap).toEqual([tileGap(16) * 2, tileGap(16) * 2]);
    expect(p.radius).toBeCloseTo(tileRadius(16, 16, tileGap(16), tileGap(16)) * 2, 6);
    expect(p.absentDarken).toBe(TILE_ABSENT_DARKEN);
  });

  test('partial fade between 6 and 10 px', () => {
    const p = tileParams(8, 8, [50, 50], true, 1);
    expect(p.mix).toBeCloseTo(0.5, 6);
    expect(p.mix).toBeGreaterThan(0);
  });
});

describe('shader failure detection', () => {
  test('matches compile/link errors, ignores others', () => {
    expect(isShaderError('Fragment shader is not compiled.')).toBe(true);
    expect(isShaderError('link error during link-error')).toBe(true);
    expect(isShaderError('GLSL program link failed')).toBe(true);
    expect(isShaderError('network request failed')).toBe(false);
  });
});
