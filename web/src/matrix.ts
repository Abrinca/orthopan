// Presence/absence heatmap: display colormap and the deck.gl BitmapLayer.
// Viewport math lives in ./transform so every panel shares one mapping.

import {Deck, OrthographicView} from '@deck.gl/core';
import {BitmapLayer} from '@deck.gl/layers';
import {copiesBucket} from './decode';
import type {Theme} from './theme';
import {paletteFor} from './theme';
import type {View} from './transform';

// Sparse windows produce tiny non-zero fractions; a gamma lift keeps them
// visible while keeping 0 pure white/dark (see the LUT below).
const GAMMA = 0.45;
const NONZERO_FLOOR = 0.12;

function stopColorAt(stops: [number, number, number][], t: number): [number, number, number] {
  const u = Math.min(1, Math.max(0, t)) * (stops.length - 1);
  const a = Math.min(stops.length - 2, Math.floor(u));
  const f = u - a;
  return [
    Math.round(stops[a][0] * (1 - f) + stops[a + 1][0] * f),
    Math.round(stops[a][1] * (1 - f) + stops[a + 1][1] * f),
    Math.round(stops[a][2] * (1 - f) + stops[a + 1][2] * f),
  ];
}

/** Map a raw 0..255 fraction to the displayed ramp position (theme-independent). */
export function displayT(v: number): number {
  if (v <= 0) return 0;
  return NONZERO_FLOOR + (1 - NONZERO_FLOOR) * Math.pow(v / 255, GAMMA);
}

export function colorForValue(v: number, theme: Theme = 'light'): [number, number, number] {
  return stopColorAt(paletteFor(theme).presenceStops, displayT(v));
}

const LUT_CACHE = new Map<Theme, Uint8Array>();
function presenceLut(theme: Theme): Uint8Array {
  const cached = LUT_CACHE.get(theme);
  if (cached) return cached;
  const stops = paletteFor(theme).presenceStops;
  const lut = new Uint8Array(256 * 3);
  for (let v = 0; v < 256; v++) {
    const [r, g, b] = stopColorAt(stops, displayT(v));
    lut[v * 3] = r;
    lut[v * 3 + 1] = g;
    lut[v * 3 + 2] = b;
  }
  LUT_CACHE.set(theme, lut);
  return lut;
}

/** Expand one-channel bin values into the RGBA image the BitmapLayer uploads. */
export function toImageData(values: Uint8Array, w: number, h: number, theme: Theme = 'light'): ImageData {
  const lut = presenceLut(theme);
  const rgba = new Uint8ClampedArray(w * h * 4);
  const n = Math.min(values.length, w * h);
  for (let i = 0, j = 0; i < n; i++, j += 4) {
    const c = values[i] * 3;
    rgba[j] = lut[c];
    rgba[j + 1] = lut[c + 1];
    rgba[j + 2] = lut[c + 2];
    // alpha is only a presence mask for the tiled shader; the flat path ignores
    // it, so the rendering below the activation threshold is unchanged.
    rgba[j + 3] = values[i] > 0 ? 255 : 0;
  }
  return new ImageData(rgba, w, h);
}

// --- copy-number discrete colormap -----------------------------------------

export function copiesLegend(theme: Theme = 'light'): {level: number; label: string; color: string}[] {
  return paletteFor(theme).copiesColors.map((c, i) => ({
    level: i,
    label: i === 0 ? 'absent' : i === 5 ? '5+' : String(i),
    color: `rgb(${c[0]},${c[1]},${c[2]})`,
  }));
}

/** Expand copy-number bin bytes (raw = round(10*mean)) into an RGBA image. */
export function toCopiesImageData(values: Uint8Array, w: number, h: number, theme: Theme = 'light'): ImageData {
  const colors = paletteFor(theme).copiesColors;
  const rgba = new Uint8ClampedArray(w * h * 4);
  const n = Math.min(values.length, w * h);
  for (let i = 0, j = 0; i < n; i++, j += 4) {
    const bucket = copiesBucket(values[i]);
    const c = colors[bucket];
    rgba[j] = c[0];
    rgba[j + 1] = c[1];
    rgba[j + 2] = c[2];
    // presence mask only (see toImageData); copies 0 = absent
    rgba[j + 3] = bucket > 0 ? 255 : 0;
  }
  return new ImageData(rgba, w, h);
}

export function toImageDataFor(
  mode: 'presence' | 'copies',
  values: Uint8Array,
  w: number,
  h: number,
  theme: Theme = 'light',
): ImageData {
  return mode === 'copies' ? toCopiesImageData(values, w, h, theme) : toImageData(values, w, h, theme);
}

// --- rounded-tile style -----------------------------------------------------

/** Below this CSS-px cell size the flat rendering is used unchanged. */
export const TILE_MIN_PX = 6;
/** Fully tiled at/above this CSS-px cell size; smoothstep in between. */
export const TILE_FULL_PX = 10;
/** Gap is ~12% of the cell, clamped to [1, 3] CSS px (gap between neighbours). */
export const TILE_GAP_FRAC = 0.12;
export const TILE_GAP_MIN = 1;
export const TILE_GAP_MAX = 3;
/** Corner radius is ~25% of the smaller cell dimension. */
export const TILE_RADIUS_FRAC = 0.25;
/** Absent tiles multiply the sampled colour by this (slightly darker tile). */
export const TILE_ABSENT_DARKEN = 0.955;

/** Smoothstep fade for the tiled style from TILE_MIN_PX to TILE_FULL_PX. */
export function tileFade(minCellPx: number): number {
  if (minCellPx <= TILE_MIN_PX) return 0;
  if (minCellPx >= TILE_FULL_PX) return 1;
  const t = (minCellPx - TILE_MIN_PX) / (TILE_FULL_PX - TILE_MIN_PX);
  return t * t * (3 - 2 * t);
}

/** Visible gap between neighbouring tiles, in the same units as `cellPx`. */
export function tileGap(cellPx: number): number {
  return Math.min(TILE_GAP_MAX, Math.max(TILE_GAP_MIN, TILE_GAP_FRAC * cellPx));
}

/** Corner radius, never so large that the tile collapses into the gap. */
export function tileRadius(cellW: number, cellH: number, gapX: number, gapY: number): number {
  const minDim = Math.min(cellW, cellH);
  const maxR = Math.max(0, minDim / 2 - Math.max(gapX, gapY) / 2 - 0.5);
  return Math.min(TILE_RADIUS_FRAC * minDim, maxR);
}

export interface TileParams {
  /** texture size in cells (x = orthogroup columns, y = genome rows) */
  cellCount: [number, number];
  /** framebuffer pixels per cell */
  cellSize: [number, number];
  /** framebuffer-pixel gap between neighbours */
  gap: [number, number];
  /** framebuffer-pixel corner radius */
  radius: number;
  /** 0 = flat (today's look), 1 = fully tiled */
  mix: number;
  /** absent tile colour multiplier */
  absentDarken: number;
}

export const FLAT_TILE: TileParams = {
  cellCount: [1, 1],
  cellSize: [0, 0],
  gap: [0, 0],
  radius: 0,
  mix: 0,
  absentDarken: TILE_ABSENT_DARKEN,
};

/**
 * Tile style for a viewport. `cellW`/`cellH` are CSS px per cell, `fullDetail`
 * is true only when the window is not binned on either axis, and `dpr` converts
 * to framebuffer pixels.
 */
export function tileParams(
  cellW: number,
  cellH: number,
  cellCount: [number, number],
  fullDetail: boolean,
  dpr: number,
  absentDarken: number = TILE_ABSENT_DARKEN,
): TileParams {
  const mix = fullDetail ? tileFade(Math.min(cellW, cellH)) : 0;
  if (mix <= 0) return {...FLAT_TILE, cellCount, absentDarken};
  const gapX = tileGap(cellW);
  const gapY = tileGap(cellH);
  return {
    cellCount,
    cellSize: [cellW * dpr, cellH * dpr],
    gap: [gapX * dpr, gapY * dpr],
    radius: tileRadius(cellW, cellH, gapX, gapY) * dpr,
    mix,
    absentDarken,
  };
}

// --- deck.gl heatmap --------------------------------------------------------

// Custom uniforms for the rounded-tile fragment shader. Layout mirrors the
// bitmap module's uniform block (luma.gl computes std140 offsets from the types).
const TILE_UNIFORM_GLSL = `\
layout(std140) uniform tileUniforms {
  vec2 cellCount;
  vec2 cellSize;
  vec2 gap;
  float radius;
  float mixAmount;
  float absentDarken;
} tile;
`;

const TILE_UNIFORM_WGSL = `\
struct TileUniforms {
  cellCount: vec2<f32>,
  cellSize: vec2<f32>,
  gap: vec2<f32>,
  radius: f32,
  mixAmount: f32,
  absentDarken: f32,
};
@group(0) @binding(auto) var<uniform> tile: TileUniforms;
`;

export const tileModule = {
  name: 'tile',
  source: TILE_UNIFORM_WGSL,
  vs: TILE_UNIFORM_GLSL,
  fs: TILE_UNIFORM_GLSL,
  uniformTypes: {
    cellCount: 'vec2<f32>',
    cellSize: 'vec2<f32>',
    gap: 'vec2<f32>',
    radius: 'f32',
    mixAmount: 'f32',
    absentDarken: 'f32',
  },
} as const;

// BitmapLayer's fragment shader with a per-cell rounded tile applied in
// fragment space. `mixAmount = 0` reproduces the flat rendering exactly.
const TILE_FS = `\
#version 300 es
#define SHADER_NAME tiled-bitmap-fragment-shader

#ifdef GL_ES
precision highp float;
#endif

uniform sampler2D bitmapTexture;

in vec2 vTexCoord;
in vec2 vTexPos;

out vec4 fragColor;

vec3 color_desaturate(vec3 color) {
  float luminance = (color.r + color.g + color.b) * 0.333333333;
  return mix(color, vec3(luminance), bitmap.desaturate);
}

vec3 color_tint(vec3 color) {
  return color * bitmap.tintColor;
}

vec4 apply_opacity(vec3 color, float alpha) {
  if (bitmap.transparentColor.a == 0.0) {
    return vec4(color, alpha);
  }
  float blendedAlpha = alpha + bitmap.transparentColor.a * (1.0 - alpha);
  float highLightRatio = alpha / blendedAlpha;
  vec3 blendedRGB = mix(bitmap.transparentColor.rgb, color, highLightRatio);
  return vec4(blendedRGB, blendedAlpha);
}

float sd_round_box(vec2 p, vec2 b, float r) {
  vec2 q = abs(p) - b + r;
  return min(max(q.x, q.y), 0.0) + length(max(q, 0.0)) - r;
}

void main(void) {
  vec2 uv = vTexCoord;
  vec4 bitmapColor = texture(bitmapTexture, uv);
  vec3 base = color_tint(color_desaturate(bitmapColor.rgb));

  geometry.uv = uv;

  // Below the activation threshold keep today's flat rendering exactly.
  if (tile.mixAmount <= 0.0) {
    fragColor = apply_opacity(base, layer.opacity);
    DECKGL_FILTER_COLOR(fragColor, geometry);
    return;
  }

  // alpha is a presence mask (see toImageData)
  float present = step(0.5, bitmapColor.a);
  // absent tiles darken smoothly as the tile style fades in
  float darken = mix(1.0, tile.absentDarken, tile.mixAmount);
  vec3 tileRgb = base * mix(darken, 1.0, present);

  // Fade the gap width (not the gap colour): gaps always show the background.
  vec2 halfCell = tile.cellSize * 0.5;
  vec2 b = halfCell - tile.gap * 0.5 * tile.mixAmount;
  float r = clamp(tile.radius * tile.mixAmount, 0.0, max(0.0, min(b.x, b.y)));
  vec2 cell = fract(uv * tile.cellCount);
  vec2 p = (cell - 0.5) * tile.cellSize;
  float sd = sd_round_box(p, b, r);
  float aa = clamp(fwidth(sd), 0.75, 2.0);
  float cov = clamp(0.5 - sd / aa, 0.0, 1.0);

  fragColor = vec4(tileRgb, cov * layer.opacity);
  DECKGL_FILTER_COLOR(fragColor, geometry);
}
`;

type TileFn = () => TileParams;

class TiledBitmapLayer extends BitmapLayer<{tileParams?: TileFn}> {
  static layerName = 'TiledBitmapLayer';
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  static defaultProps: any = {
    ...BitmapLayer.defaultProps,
    tileParams: {type: 'function', value: () => FLAT_TILE},
  };

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  getShaders(): any {
    const shaders = super.getShaders();
    const modules = [...(shaders.modules ?? []), tileModule];
    return {...shaders, fs: TILE_FS, modules};
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  draw(opts: any): void {
    const tp = this.props.tileParams?.() ?? FLAT_TILE;
    (this as any).setShaderModuleProps({
      tile: {
        cellCount: tp.cellCount,
        cellSize: tp.cellSize,
        gap: tp.gap,
        radius: tp.radius,
        mixAmount: tp.mix,
        absentDarken: tp.absentDarken,
      },
    });
    super.draw(opts);
  }
}

/** True when an error looks like a shader compile/link failure. */
export function isShaderError(message: string): boolean {
  return /shader|compile|link|program|glsl|wgsl/i.test(message);
}

declare global {
  interface Window {
    /** installed by the inline head script; fires on a shader compile/link error */
    __onShaderError?: ((message: string) => void) | null;
  }
}

/** Copy of the texture with opaque alpha (the flat rendering's expectation). */
function opaqueImage(img: ImageData): ImageData {
  const d = new Uint8ClampedArray(img.data);
  for (let i = 3; i < d.length; i += 4) d[i] = 255;
  return new ImageData(d, img.width, img.height);
}

export class Heatmap {
  readonly deck: Deck<OrthographicView>;
  /** set once the tiled shader fails, so we stay on the plain BitmapLayer path */
  private flat = false;
  private lastImage: ImageData | null = null;
  private lastBounds: [number, number, number, number] | null = null;
  private lastTileParams?: TileFn;

  constructor(
    container: HTMLDivElement,
    private constrain: (v: View) => View,
    private onView: (v: View) => void,
  ) {
    // the inline head script installs the console.error hook (before luma loads)
    if (typeof window !== 'undefined') window.__onShaderError = (msg: string) => this.handleError(new Error(msg));
    this.deck = new Deck<OrthographicView>({
      parent: container,
      views: new OrthographicView({id: 'ortho'}),
      viewState: viewStateOf({cx: 0, cy: 0, zx: 0, zy: 0}),
      controller: {dragPan: true, scrollZoom: false, doubleClickZoom: false, touchZoom: false},
      onError: (error: Error) => this.handleError(error),
      onViewStateChange: (params) => {
        const vs = params.viewState;
        const v: View = {
          cx: vs.target?.[0] ?? 0,
          cy: vs.target?.[1] ?? 0,
          zx: vs.zoomX ?? (typeof vs.zoom === 'number' ? vs.zoom : 0),
          zy: vs.zoomY ?? (typeof vs.zoom === 'number' ? vs.zoom : 0),
        };
        const clamped = this.constrain(v);
        this.onView(clamped);
        // controlled mode: push the clamped state back so the controller continues from it
        this.deck.setProps({viewState: viewStateOf(clamped)});
      },
      layers: [],
    });
  }

  /** True once the fallback flat rendering is active. */
  get usingFlat(): boolean {
    return this.flat;
  }

  /** Fall back permanently to the plain BitmapLayer if the tiled shader fails. */
  handleError(error: Error): void {
    const msg = `${error?.message ?? ''} ${String((error as {cause?: unknown})?.cause ?? '')}`;
    if (this.flat || !isShaderError(msg)) return;
    this.flat = true;
    // eslint-disable-next-line no-console
    console.warn('OrthoPan: tiled matrix shader failed, falling back to flat rendering —', msg.slice(0, 200));
    this.applyLayer();
  }

  setView(v: View): void {
    this.deck.setProps({viewState: viewStateOf(v)});
  }

  /** Force the WebGL canvas to re-read devicePixelRatio after a DPR-only change. */
  resyncDpr(): void {
    this.deck.setProps({useDevicePixels: window.devicePixelRatio || 1});
    this.deck.redraw('dpr-change');
  }

  setData(image: ImageData, bounds: [number, number, number, number], tileParams?: TileFn): void {
    this.lastImage = image;
    this.lastBounds = bounds;
    this.lastTileParams = tileParams;
    this.applyLayer();
  }

  private applyLayer(): void {
    if (!this.lastImage || !this.lastBounds) return;
    const base = {
      id: 'matrix',
      bounds: this.lastBounds,
      pickable: false,
      textureParameters: {magFilter: 'nearest' as const, minFilter: 'nearest' as const},
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      onError: ((e: Error) => this.handleError(e)) as any,
    };
    const layer = this.flat
      ? new BitmapLayer({...base, image: opaqueImage(this.lastImage)})
      : new TiledBitmapLayer({...base, image: this.lastImage, tileParams: this.lastTileParams});
    this.deck.setProps({layers: [layer]});
  }
}

function viewStateOf(v: View) {
  return {target: [v.cx, v.cy, 0] as [number, number, number], zoomX: v.zx, zoomY: v.zy};
}
