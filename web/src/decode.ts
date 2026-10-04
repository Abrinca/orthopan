// Little-endian decoders and copy-number value math (docs/api.md v1).

/** Decode an ArrayBuffer of little-endian uint16 values. */
export function decodeUint16LE(buf: ArrayBuffer): Uint16Array {
  const n = Math.floor(buf.byteLength / 2);
  const out = new Uint16Array(n);
  const dv = new DataView(buf);
  for (let i = 0; i < n; i++) out[i] = dv.getUint16(i * 2, true);
  return out;
}

/** Decode an ArrayBuffer of little-endian float32 values. */
export function decodeFloat32LE(buf: ArrayBuffer): Float32Array {
  const n = Math.floor(buf.byteLength / 4);
  const out = new Float32Array(n);
  const dv = new DataView(buf);
  for (let i = 0; i < n; i++) out[i] = dv.getFloat32(i * 4, true);
  return out;
}

/**
 * Copy-number matrix contract: 0 = no present cell, otherwise the raw byte is
 * `round(10 * mean copies over present cells)`. These helpers keep the decode
 * in one place for the image, legend and tooltip.
 */
export function copiesMean(v: number): number {
  return v / 10;
}

/** Discrete display bucket: 0 absent, 1..5 where 5 means "5 or more". */
export function copiesBucket(v: number): number {
  if (v <= 0) return 0;
  return Math.max(1, Math.min(5, Math.round(v / 10)));
}

// A raw copy byte is saturated at 255 (>= 25.5 mean copies): report it as a
// lower bound rather than an exact value.
export const COPIES_SATURATED = 255;

/** Mean copies for a binned cell (2 decimals), or "≥25.5" when saturated. */
export function formatCopiesMean(v: number): string {
  if (v <= 0) return '0';
  if (v >= COPIES_SATURATED) return '≥25.5';
  return copiesMean(v).toFixed(2);
}

/** Integer copies for a single cell, or "≥25.5" when saturated. */
export function formatCopiesExact(v: number): string {
  if (v <= 0) return '0';
  if (v >= COPIES_SATURATED) return '≥25.5';
  return String(Math.round(copiesMean(v)));
}
