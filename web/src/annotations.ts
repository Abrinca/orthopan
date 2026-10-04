// Categorical annotation palette and per-value colour mapping.
// Top-N categories get distinct colours; the rest fall back to grey.

export const CATEGORY_COLORS = [
  '#4e79a7',
  '#f28e2b',
  '#e15759',
  '#76b7b2',
  '#59a14f',
  '#edc948',
  '#b07aa1',
  '#ff9da7',
  '#9c755f',
  '#8c6d31',
  '#86bcb6',
  '#d37295',
];

export const CATEGORY_GREY = '#c7ccd4';
export const CATEGORY_MISSING = '#f2f4f7';

/** Colour for a category index (0-based); beyond the palette it is grey. */
export function categoryColor(index: number): string {
  return index >= 0 && index < CATEGORY_COLORS.length ? CATEGORY_COLORS[index] : CATEGORY_GREY;
}

/** Colour for an annotation-track value: 0 = missing, else 1 + category index. */
export function annotationStripColor(value: number): string {
  return value <= 0 ? CATEGORY_MISSING : categoryColor(value - 1);
}

/** Category label for a track value (categories are ordered by frequency). */
export function annotationCategory(value: number, categories: string[]): string {
  if (value <= 0) return '(missing)';
  return categories[value - 1] ?? String(value);
}
