import {describe, expect, test} from 'bun:test';
import {copiesLegend, colorForValue, displayT} from '../src/matrix';
import {categoryColor} from '../src/annotations';
import {residueColor} from '../src/msa';
import {LIGHT_PALETTE, DARK_PALETTE, paletteFor, resolveTheme} from '../src/theme';

const sum = (c: [number, number, number]) => c[0] + c[1] + c[2];

describe('theme resolution', () => {
  test('auto follows the OS preference', () => {
    expect(resolveTheme('auto', true)).toBe('dark');
    expect(resolveTheme('auto', false)).toBe('light');
  });

  test('explicit light/dark ignore the OS preference', () => {
    expect(resolveTheme('light', true)).toBe('light');
    expect(resolveTheme('light', false)).toBe('light');
    expect(resolveTheme('dark', true)).toBe('dark');
    expect(resolveTheme('dark', false)).toBe('dark');
  });

  test('an OS preference change only affects auto', () => {
    // same setting, flipped preference
    expect(resolveTheme('auto', false)).not.toBe(resolveTheme('auto', true));
    expect(resolveTheme('dark', false)).toBe(resolveTheme('dark', true));
    expect(resolveTheme('light', false)).toBe(resolveTheme('light', true));
  });

  test('paletteFor maps themes', () => {
    expect(paletteFor('light')).toBe(LIGHT_PALETTE);
    expect(paletteFor('dark')).toBe(DARK_PALETTE);
  });
});

describe('colormap endpoints per theme', () => {
  test('presence: light absent is near-white, present is dark', () => {
    expect(sum(colorForValue(0, 'light'))).toBeGreaterThan(720);
    expect(sum(colorForValue(255, 'light'))).toBeLessThan(300);
  });

  test('presence: dark absent is near the background, present is bright', () => {
    expect(sum(colorForValue(0, 'dark'))).toBeLessThan(140);
    expect(sum(colorForValue(255, 'dark'))).toBeGreaterThan(560);
    // value 0 is the ramp start in both themes
    expect(colorForValue(0, 'dark')).toEqual(DARK_PALETTE.presenceStops[0]);
    expect(colorForValue(255, 'dark')).toEqual(DARK_PALETTE.presenceStops[DARK_PALETTE.presenceStops.length - 1]);
  });

  test('gamma still lifts sparse presence on both themes', () => {
    expect(displayT(1)).toBeGreaterThan(0);
    expect(sum(colorForValue(1, 'dark'))).toBeGreaterThan(sum(colorForValue(0, 'dark')));
  });

  test('copies legend has 6 levels; absent follows the theme', () => {
    const light = copiesLegend('light');
    const dark = copiesLegend('dark');
    expect(light).toHaveLength(6);
    expect(dark).toHaveLength(6);
    expect(light[0].color).toBe('rgb(255,255,255)');
    // dark absent is near the dark background, not white
    expect(dark[0].color).toMatch(/^rgb\(\d+,\d+,\d+\)$/);
    expect(dark[0].color).not.toBe('rgb(255,255,255)');
  });
});

describe('categorical + residue palettes', () => {
  test('dark categories stay distinct', () => {
    const cols = DARK_PALETTE.categories;
    expect(cols.length).toBeGreaterThanOrEqual(12);
    expect(new Set(cols).size).toBe(cols.length);
    expect(categoryColor(0, 'dark')).not.toBe(categoryColor(0, 'light'));
  });

  test('residue hues differ per theme and letters stay distinguishable', () => {
    expect(residueColor('A', 'dark')).not.toBe(residueColor('A', 'light'));
    expect(residueColor('D', 'dark')).not.toBe(residueColor('A', 'dark'));
    expect(residueColor('K', 'dark')).not.toBe(residueColor('D', 'dark'));
    // gaps take the theme background-ish gap colour
    expect(residueColor('-', 'dark')).toBe(DARK_PALETTE.gapResidue);
    expect(residueColor('-', 'light')).toBe(LIGHT_PALETTE.gapResidue);
  });

  test('residue letters contrast with their fills', () => {
    // dark letters on the light-ish residue fills
    const ratio = (a: [number, number, number], b: string) => {
      const lum = (rgb: [number, number, number]) => {
        const f = (v: number) => {
          const x = v / 255;
          return x <= 0.03928 ? x / 12.92 : Math.pow((x + 0.055) / 1.055, 2.4);
        };
        return 0.2126 * f(rgb[0]) + 0.7152 * f(rgb[1]) + 0.0722 * f(rgb[2]);
      };
      const hex = (h: string): [number, number, number] => [
        parseInt(h.slice(1, 3), 16),
        parseInt(h.slice(3, 5), 16),
        parseInt(h.slice(5, 7), 16),
      ];
      const la = lum(a);
      const lb = lum(hex(b));
      return (Math.max(la, lb) + 0.05) / (Math.min(la, lb) + 0.05);
    };
    for (const theme of ['light', 'dark'] as const) {
      const pal = paletteFor(theme);
      const letter = hexToRgb(pal.residueLetter);
      for (const fill of Object.values(pal.residue)) {
        expect(ratio(letter, fill)).toBeGreaterThanOrEqual(4.5);
      }
    }
  });
});

function hexToRgb(h: string): [number, number, number] {
  return [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
}

function luminance(c: [number, number, number]): number {
  const f = (v: number) => {
    const x = v / 255;
    return x <= 0.03928 ? x / 12.92 : Math.pow((x + 0.055) / 1.055, 2.4);
  };
  return 0.2126 * f(c[0]) + 0.7152 * f(c[1]) + 0.0722 * f(c[2]);
}

function assertMonotonic(vals: number[]): void {
  const dir = Math.sign(vals[1] - vals[0]);
  expect(dir).not.toBe(0);
  for (let i = 1; i < vals.length; i++) {
    expect(Math.sign(vals[i] - vals[i - 1])).toBe(dir);
    expect(Math.abs(vals[i] - vals[i - 1])).toBeGreaterThan(0);
  }
}

describe('colormap luminance is monotonic', () => {
  for (const theme of ['light', 'dark'] as const) {
    test(`presence ${theme}`, () => {
      const pal = paletteFor(theme);
      assertMonotonic(pal.presenceStops.map(luminance));
      const samples = [0, 1, 4, 16, 64, 128, 192, 255].map((v) => luminance(colorForValue(v, theme)));
      assertMonotonic(samples);
    });

    test(`copies ${theme}`, () => {
      assertMonotonic(paletteFor(theme).copiesColors.map(luminance));
    });
  }
});
