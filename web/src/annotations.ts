// Categorical annotation palette and per-value colour mapping.
// Top-N categories get distinct colours; the rest fall back to grey.

import type {Theme} from './theme';
import {paletteFor} from './theme';

/** Categorical palette for a theme (index 0..N-1). */
export function categoryColors(theme: Theme = 'light'): string[] {
  return paletteFor(theme).categories;
}

export function categoryGrey(theme: Theme = 'light'): string {
  return paletteFor(theme).categoryGrey;
}

export function categoryMissing(theme: Theme = 'light'): string {
  return paletteFor(theme).categoryMissing;
}

/** Legacy light-theme palette (kept for compatibility). */
export const CATEGORY_COLORS = categoryColors('light');
export const CATEGORY_GREY = categoryGrey('light');
export const CATEGORY_MISSING = categoryMissing('light');

/** Colour for a category index (0-based); beyond the palette it is grey. */
export function categoryColor(index: number, theme: Theme = 'light'): string {
  const cols = categoryColors(theme);
  return index >= 0 && index < cols.length ? cols[index] : categoryGrey(theme);
}

/** Colour for an annotation-track value: 0 = missing, else 1 + category index. */
export function annotationStripColor(value: number, theme: Theme = 'light'): string {
  return value <= 0 ? categoryMissing(theme) : categoryColor(value - 1, theme);
}

/** Category label for a track value (categories are ordered by frequency). */
export function annotationCategory(value: number, categories: string[]): string {
  if (value <= 0) return '(missing)';
  return categories[value - 1] ?? String(value);
}
