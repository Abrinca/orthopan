// Presence/absence heatmap: display colormap and the deck.gl BitmapLayer.
// Viewport math lives in ./transform so every panel shares one mapping.

import {Deck, OrthographicView} from '@deck.gl/core';
import {BitmapLayer} from '@deck.gl/layers';
import {copiesBucket} from './decode';
import type {View} from './transform';

// Sparse windows produce tiny non-zero fractions; a gamma lift keeps them
// visible while keeping 0 pure white and 255 dark (see COLOR_LUT below).
const GAMMA = 0.45;
const NONZERO_FLOOR = 0.12;

const STOPS: [number, number, number][] = [
  [247, 251, 255],
  [198, 219, 239],
  [107, 174, 214],
  [33, 113, 181],
  [8, 48, 107],
];

function colorAt(t: number): [number, number, number] {
  const u = Math.min(1, Math.max(0, t)) * (STOPS.length - 1);
  const a = Math.min(STOPS.length - 2, Math.floor(u));
  const f = u - a;
  return [
    Math.round(STOPS[a][0] * (1 - f) + STOPS[a + 1][0] * f),
    Math.round(STOPS[a][1] * (1 - f) + STOPS[a + 1][1] * f),
    Math.round(STOPS[a][2] * (1 - f) + STOPS[a + 1][2] * f),
  ];
}

/** Map a raw 0..255 fraction to the displayed ramp position. */
export function displayT(v: number): number {
  if (v <= 0) return 0;
  return NONZERO_FLOOR + (1 - NONZERO_FLOOR) * Math.pow(v / 255, GAMMA);
}

export function colorForValue(v: number): [number, number, number] {
  return colorAt(displayT(v));
}

const COLOR_LUT = (() => {
  const lut = new Uint8Array(256 * 3);
  for (let v = 0; v < 256; v++) {
    const [r, g, b] = colorAt(displayT(v));
    lut[v * 3] = r;
    lut[v * 3 + 1] = g;
    lut[v * 3 + 2] = b;
  }
  return lut;
})();

/** Expand one-channel bin values into the RGBA image the BitmapLayer uploads. */
export function toImageData(values: Uint8Array, w: number, h: number): ImageData {
  const rgba = new Uint8ClampedArray(w * h * 4);
  const n = Math.min(values.length, w * h);
  for (let i = 0, j = 0; i < n; i++, j += 4) {
    const c = values[i] * 3;
    rgba[j] = COLOR_LUT[c];
    rgba[j + 1] = COLOR_LUT[c + 1];
    rgba[j + 2] = COLOR_LUT[c + 2];
    rgba[j + 3] = 255;
  }
  return new ImageData(rgba, w, h);
}

// --- copy-number discrete colormap -----------------------------------------

// bucket 0 (absent) is white; 1..5 are "1, 2, 3, 4, 5+" copies
const COPIES_COLORS: [number, number, number][] = [
  [255, 255, 255],
  [254, 227, 145],
  [254, 196, 79],
  [254, 153, 41],
  [217, 95, 14],
  [153, 52, 4],
];

export function copiesLegend(): {level: number; label: string; color: string}[] {
  return COPIES_COLORS.map((c, i) => ({
    level: i,
    label: i === 0 ? 'absent' : i === 5 ? '5+' : String(i),
    color: `rgb(${c[0]},${c[1]},${c[2]})`,
  }));
}

/** Expand copy-number bin bytes (raw = round(10*mean)) into an RGBA image. */
export function toCopiesImageData(values: Uint8Array, w: number, h: number): ImageData {
  const rgba = new Uint8ClampedArray(w * h * 4);
  const n = Math.min(values.length, w * h);
  for (let i = 0, j = 0; i < n; i++, j += 4) {
    const c = COPIES_COLORS[copiesBucket(values[i])];
    rgba[j] = c[0];
    rgba[j + 1] = c[1];
    rgba[j + 2] = c[2];
    rgba[j + 3] = 255;
  }
  return new ImageData(rgba, w, h);
}

export function toImageDataFor(
  mode: 'presence' | 'copies',
  values: Uint8Array,
  w: number,
  h: number,
): ImageData {
  return mode === 'copies' ? toCopiesImageData(values, w, h) : toImageData(values, w, h);
}

// --- deck.gl heatmap --------------------------------------------------------

export class Heatmap {
  readonly deck: Deck<OrthographicView>;

  constructor(
    container: HTMLDivElement,
    private constrain: (v: View) => View,
    private onView: (v: View) => void,
  ) {
    this.deck = new Deck<OrthographicView>({
      parent: container,
      views: new OrthographicView({id: 'ortho'}),
      viewState: viewStateOf({cx: 0, cy: 0, zx: 0, zy: 0}),
      controller: {dragPan: true, scrollZoom: false, doubleClickZoom: false, touchZoom: false},
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

  setView(v: View): void {
    this.deck.setProps({viewState: viewStateOf(v)});
  }

  /** Force the WebGL canvas to re-read devicePixelRatio after a DPR-only change. */
  resyncDpr(): void {
    this.deck.setProps({useDevicePixels: window.devicePixelRatio || 1});
    this.deck.redraw('dpr-change');
  }

  setData(image: ImageData, bounds: [number, number, number, number]): void {
    this.deck.setProps({
      layers: [
        new BitmapLayer({
          id: 'matrix',
          image,
          bounds,
          pickable: false,
          textureParameters: {magFilter: 'nearest', minFilter: 'nearest'},
        }),
      ],
    });
  }
}

function viewStateOf(v: View) {
  return {target: [v.cx, v.cy, 0] as [number, number, number], zoomX: v.zx, zoomY: v.zy};
}
