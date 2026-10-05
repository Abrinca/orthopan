// Theme resolution and the single source of colour for the whole UI.
//
// The resolved theme is one of 'light' | 'dark'. The user setting is
// auto | light | dark; 'auto' follows `prefers-color-scheme` and reacts live.
// CSS custom properties are driven by the `data-theme` attribute, and the
// canvas/deck.gl code reads `palette()` at draw time.

export type Theme = 'light' | 'dark';
export type ThemeSetting = 'auto' | 'light' | 'dark';

/** Resolve the user setting against the OS preference. */
export function resolveTheme(setting: ThemeSetting, prefersDark: boolean): Theme {
  if (setting === 'light') return 'light';
  if (setting === 'dark') return 'dark';
  return prefersDark ? 'dark' : 'light';
}

export type Rgb = [number, number, number];

export interface Palette {
  theme: Theme;
  // UI chrome
  bg: string;
  panel: string;
  surface: string;
  border: string;
  borderSoft: string;
  text: string;
  dim: string;
  faint: string;
  accent: string;
  accentBg: string;
  accentText: string;
  inputBg: string;
  hover: string;
  tooltipBg: string;
  tooltipText: string;
  tooltipAccent: string;
  tooltipDim: string;
  tooltipBorder: string;
  notice: string;
  disabled: string;
  disabledBorder: string;
  badgeBg: string;
  badgeText: string;
  legendText: string;
  // shared canvas chrome
  treeLine: string;
  treeLineSelected: string;
  treePath: string;
  cladeTint: string;
  genomeLabel: string;
  ogLabel: string;
  trackBgA: string;
  trackBgB: string;
  trackBorder: string;
  trackBars: string[];
  trackLabel: string;
  metaLabel: string;
  // heatmap
  presenceStops: Rgb[];
  copiesColors: Rgb[];
  absentDarken: number;
  // highlights / selection
  selectionTint: string;
  hoverTint: string;
  hoverOutline: string;
  searchFill: string;
  searchStroke: string;
  // alignment view
  alnBg: string;
  alnBarBg: string;
  alnBorder: string;
  alnBarLabel: string;
  alnName: string;
  alnMsaBg: string;
  alnCons: Rgb;
  alnGapBar: Rgb;
  alnOverviewTrack: string;
  alnOverviewHandle: string;
  alnDupFill: string;
  alnDupStroke: string;
  missing: string;
  // categorical + numeric
  categories: string[];
  categoryGrey: string;
  categoryMissing: string;
  numericStops: Rgb[];
  gapResidue: string;
  residueLetter: string;
  residue: {hydro: string; acid: string; basic: string; polar: string; other: string};
}

const PRESENCE_STOPS_LIGHT: Rgb[] = [
  [247, 251, 255],
  [198, 219, 239],
  [107, 174, 214],
  [33, 113, 181],
  [8, 48, 107],
];

// dark navy -> cyan -> pale yellow: bright present cells on a dark background
const PRESENCE_STOPS_DARK: Rgb[] = [
  [22, 34, 54],
  [40, 92, 140],
  [66, 145, 190],
  [120, 205, 225],
  [235, 246, 175],
];

const COPIES_LIGHT: Rgb[] = [
  [255, 255, 255],
  [254, 227, 145],
  [254, 196, 79],
  [254, 153, 41],
  [217, 95, 14],
  [153, 52, 4],
];

const COPIES_DARK: Rgb[] = [
  [26, 36, 56],
  [47, 104, 168],
  [60, 160, 180],
  [140, 205, 160],
  [236, 222, 120],
  [255, 248, 200],
];

const CATEGORIES_LIGHT = [
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

// brighter / more saturated so the hues stay distinct on a dark background
const CATEGORIES_DARK = [
  '#6ea8fe',
  '#f2a65a',
  '#e86a6a',
  '#63c7c0',
  '#7cc46a',
  '#e6d357',
  '#c79be0',
  '#ff9dc0',
  '#c9a27e',
  '#d4b06a',
  '#8fd0c8',
  '#e08fb0',
];

export const LIGHT_PALETTE: Palette = {
  theme: 'light',
  bg: '#ffffff',
  panel: '#f8fafc',
  surface: '#ffffff',
  border: '#d0d5dd',
  borderSoft: '#eef2f6',
  text: '#1f2937',
  dim: '#475467',
  faint: '#667085',
  accent: '#2563eb',
  accentBg: '#2563eb',
  accentText: '#ffffff',
  inputBg: '#ffffff',
  hover: '#eff4ff',
  tooltipBg: 'rgba(17,24,39,0.92)',
  tooltipText: '#f9fafb',
  tooltipAccent: '#93c5fd',
  tooltipDim: '#cbd5e1',
  tooltipBorder: 'transparent',
  notice: '#b54708',
  disabled: '#98a2b3',
  disabledBorder: '#e4e7ec',
  badgeBg: '#fef0c7',
  badgeText: '#b54708',
  legendText: '#1f2937',
  treeLine: '#475569',
  treeLineSelected: '#334155',
  treePath: '#2563eb',
  cladeTint: 'rgba(37,99,235,0.12)',
  genomeLabel: '#344054',
  ogLabel: '#475467',
  trackBgA: '#f5f8fc',
  trackBgB: '#fbfcfe',
  trackBorder: '#e4e7ec',
  trackBars: ['#2563eb', '#059669', '#d97706'],
  trackLabel: '#344054',
  metaLabel: '#475467',
  presenceStops: PRESENCE_STOPS_LIGHT,
  copiesColors: COPIES_LIGHT,
  absentDarken: 0.955,
  selectionTint: 'rgba(245,158,11,0.14)',
  hoverTint: 'rgba(37,99,235,0.10)',
  hoverOutline: 'rgba(17,24,39,0.7)',
  searchFill: 'rgba(37,99,235,0.20)',
  searchStroke: 'rgba(37,99,235,0.9)',
  alnBg: '#ffffff',
  alnBarBg: '#f8fafc',
  alnBorder: '#e4e7ec',
  alnBarLabel: '#475467',
  alnName: '#1f2937',
  alnMsaBg: '#f3f4f6',
  alnCons: [56, 161, 105],
  alnGapBar: [180, 83, 9],
  alnOverviewTrack: '#eef2f6',
  alnOverviewHandle: '#94a3b8',
  alnDupFill: '#ffffff',
  alnDupStroke: '#111827',
  missing: '#98a2b3',
  categories: CATEGORIES_LIGHT,
  categoryGrey: '#c7ccd4',
  categoryMissing: '#f2f4f7',
  numericStops: PRESENCE_STOPS_LIGHT,
  gapResidue: '#eceff3',
  residueLetter: '#0b1220',
  residue: {hydro: '#7fa8d9', acid: '#e08bb0', basic: '#e06666', polar: '#8fce8f', other: '#f0c869'},
};

export const DARK_PALETTE: Palette = {
  theme: 'dark',
  bg: '#0b1220',
  panel: '#111827',
  surface: '#0f172a',
  border: '#334155',
  borderSoft: '#1f2937',
  text: '#e5e7eb',
  dim: '#cbd5e1',
  faint: '#94a3b8',
  accent: '#60a5fa',
  accentBg: '#2563eb',
  accentText: '#ffffff',
  inputBg: '#1f2937',
  hover: '#1e293b',
  tooltipBg: 'rgba(2,6,23,0.95)',
  tooltipText: '#f1f5f9',
  tooltipAccent: '#93c5fd',
  tooltipDim: '#94a3b8',
  tooltipBorder: '#334155',
  notice: '#fbbf24',
  disabled: '#64748b',
  disabledBorder: '#334155',
  badgeBg: '#422006',
  badgeText: '#fbbf24',
  legendText: '#e5e7eb',
  treeLine: '#94a3b8',
  treeLineSelected: '#cbd5e1',
  treePath: '#60a5fa',
  cladeTint: 'rgba(96,165,250,0.18)',
  genomeLabel: '#cbd5e1',
  ogLabel: '#94a3b8',
  trackBgA: '#0f172a',
  trackBgB: '#111827',
  trackBorder: '#1f2937',
  trackBars: ['#60a5fa', '#34d399', '#fbbf24'],
  trackLabel: '#cbd5e1',
  metaLabel: '#94a3b8',
  presenceStops: PRESENCE_STOPS_DARK,
  copiesColors: COPIES_DARK,
  absentDarken: 1.0,
  selectionTint: 'rgba(250,204,21,0.18)',
  hoverTint: 'rgba(96,165,250,0.18)',
  hoverOutline: 'rgba(226,232,240,0.8)',
  searchFill: 'rgba(96,165,250,0.25)',
  searchStroke: 'rgba(96,165,250,0.9)',
  alnBg: '#0b1220',
  alnBarBg: '#0f172a',
  alnBorder: '#1f2937',
  alnBarLabel: '#94a3b8',
  alnName: '#e5e7eb',
  alnMsaBg: '#111827',
  alnCons: [74, 222, 128],
  alnGapBar: [251, 146, 60],
  alnOverviewTrack: '#1f2937',
  alnOverviewHandle: '#64748b',
  alnDupFill: '#0b1220',
  alnDupStroke: '#e5e7eb',
  missing: '#64748b',
  categories: CATEGORIES_DARK,
  categoryGrey: '#4b5563',
  categoryMissing: '#2b3647',
  numericStops: PRESENCE_STOPS_DARK,
  gapResidue: '#111827',
  residueLetter: '#0b1220',
  residue: {hydro: '#5b9bd5', acid: '#e07aa8', basic: '#ef6f6f', polar: '#5fc27e', other: '#e6b84c'},
};

export function paletteFor(theme: Theme): Palette {
  return theme === 'dark' ? DARK_PALETTE : LIGHT_PALETTE;
}

// --- live theme store -------------------------------------------------------

const STORAGE_KEY = 'orthopan.theme';
const prefersDark = () => typeof window !== 'undefined' && window.matchMedia?.('(prefers-color-scheme: dark)').matches === true;

let setting: ThemeSetting = 'auto';
let resolved: Theme = 'light';
const listeners = new Set<(theme: Theme) => void>();
let media: MediaQueryList | null = null;

export function currentTheme(): Theme {
  return resolved;
}

/** The current palette (canvas + deck.gl colours). */
export function palette(): Palette {
  return paletteFor(resolved);
}

export function getThemeSetting(): ThemeSetting {
  return setting;
}

export function onThemeChange(cb: (theme: Theme) => void): () => void {
  listeners.add(cb);
  return () => listeners.delete(cb);
}

function readSetting(): ThemeSetting {
  try {
    const s = localStorage.getItem(STORAGE_KEY);
    if (s === 'light' || s === 'dark' || s === 'auto') return s;
  } catch {
    /* storage unavailable */
  }
  return 'auto';
}

function notify(): void {
  for (const cb of listeners) cb(resolved);
}

/** Apply the setting, update the DOM attribute and notify listeners if changed. */
export function setThemeSetting(next: ThemeSetting, opts: {persist?: boolean} = {}): void {
  setting = next;
  if (opts.persist !== false) {
    try {
      localStorage.setItem(STORAGE_KEY, next);
    } catch {
      /* storage unavailable */
    }
  }
  const r = resolveTheme(setting, prefersDark());
  applyResolved(r);
}

function applyResolved(r: Theme): void {
  const changed = r !== resolved;
  resolved = r;
  if (typeof document !== 'undefined') document.documentElement.dataset.theme = r;
  if (changed) notify();
}

function onMediaChange(): void {
  if (setting === 'auto') applyResolved(resolveTheme('auto', prefersDark()));
}

/** Initialise from storage + OS preference and start listening for OS changes. */
export function initTheme(): void {
  setting = readSetting();
  resolved = resolveTheme(setting, prefersDark());
  if (typeof document !== 'undefined') document.documentElement.dataset.theme = resolved;
  if (typeof window !== 'undefined' && window.matchMedia) {
    media = window.matchMedia('(prefers-color-scheme: dark)');
    media.addEventListener?.('change', onMediaChange);
  }
}
