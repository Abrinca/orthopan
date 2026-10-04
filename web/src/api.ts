// Typed client for the OrthoPan viewer API (docs/api.md).

export interface TrackInfo {
  name: string;
  min: number;
  max: number;
  missing: number;
}

export interface Meta {
  n_genomes: number;
  n_orthogroups: number;
  tree_source: string;
  orders: string[];
  tracks: TrackInfo[];
}

export interface Genome {
  row: number;
  genome_idx: number;
  name: string;
  n_orthogroups: number;
  n_genes: number;
}

export interface TreeArrays {
  parent: number[];
  length: number[];
  label: string[];
  row: number[];
}

export interface MatrixResult {
  data: Uint8Array;
  width: number;
  height: number;
}

export interface TracksResult {
  data: Float32Array;
  /** track names, in the order the server sent them (X-Tracks header) */
  names: string[];
  bins: number;
}

export interface OrthogroupDetails {
  orthogroup: string;
  og_idx: number;
  n_genomes: number;
  n_genes: number;
  mrca_label: string | null;
  purity: number | null;
  tracks: Record<string, number | null>;
  members: {genome: string; genes: string[]}[];
}

export class HttpError extends Error {
  constructor(
    public status: number,
    message: string,
  ) {
    super(message);
    this.name = 'HttpError';
  }
}

const BASE = '';

async function getJSON<T>(path: string, signal?: AbortSignal): Promise<T> {
  const res = await fetch(BASE + path, {signal});
  if (!res.ok) throw new HttpError(res.status, `${path}: HTTP ${res.status}`);
  return (await res.json()) as T;
}

export function getMeta(): Promise<Meta> {
  return getJSON<Meta>('/api/meta');
}

export function getGenomes(): Promise<Genome[]> {
  return getJSON<Genome[]>('/api/genomes');
}

export function getTree(): Promise<TreeArrays> {
  return getJSON<TreeArrays>('/api/tree');
}

export interface MatrixParams {
  order: string;
  og0: number;
  og1: number;
  g0: number;
  g1: number;
  w: number;
  h: number;
}

export async function getMatrix(p: MatrixParams, signal?: AbortSignal): Promise<MatrixResult> {
  const q = new URLSearchParams({
    order: p.order,
    og0: String(p.og0),
    og1: String(p.og1),
    g0: String(p.g0),
    g1: String(p.g1),
    w: String(p.w),
    h: String(p.h),
  });
  const res = await fetch(`/api/matrix?${q}`, {signal});
  if (!res.ok) throw new HttpError(res.status, `/api/matrix: HTTP ${res.status}`);
  const width = Number(res.headers.get('X-Width') ?? p.w);
  const height = Number(res.headers.get('X-Height') ?? p.h);
  const buf = await res.arrayBuffer();
  if (buf.byteLength !== width * height) {
    throw new Error(`matrix payload ${buf.byteLength} != ${width}*${height}`);
  }
  return {data: new Uint8Array(buf), width, height};
}

export interface TracksParams {
  order: string;
  og0: number;
  og1: number;
  w: number;
  names: string[];
}

export async function getTracks(p: TracksParams, signal?: AbortSignal): Promise<TracksResult> {
  const q = new URLSearchParams({
    order: p.order,
    og0: String(p.og0),
    og1: String(p.og1),
    w: String(p.w),
    names: p.names.join(','),
  });
  const res = await fetch(`/api/tracks?${q}`, {signal});
  if (!res.ok) throw new HttpError(res.status, `/api/tracks: HTTP ${res.status}`);
  const names = (res.headers.get('X-Tracks') ?? p.names.join(',')).split(',').filter(Boolean);
  const buf = await res.arrayBuffer();
  const count = names.length * p.w;
  if (buf.byteLength !== count * 4) {
    throw new Error(`tracks payload ${buf.byteLength} != ${count}*4`);
  }
  const dv = new DataView(buf);
  const data = new Float32Array(count);
  for (let i = 0; i < count; i++) data[i] = dv.getFloat32(i * 4, true);
  return {data, names, bins: p.w};
}

export async function getOrthogroups(
  order: string,
  og0: number,
  og1: number,
  signal?: AbortSignal,
): Promise<string[]> {
  if (og1 - og0 > 5000) throw new Error('orthogroup range exceeds 5000');
  const q = new URLSearchParams({order, og0: String(og0), og1: String(og1)});
  return getJSON<string[]>(`/api/orthogroups?${q}`, signal);
}

export function getOrthogroup(og: string, signal?: AbortSignal): Promise<OrthogroupDetails> {
  return getJSON<OrthogroupDetails>(`/api/orthogroup/${encodeURIComponent(og)}`, signal);
}

export async function postOrder(body: {
  by: string;
  genomes?: string[];
}): Promise<{order: string}> {
  const res = await fetch('/api/orders', {
    method: 'POST',
    headers: {'content-type': 'application/json'},
    body: JSON.stringify(body),
  });
  if (!res.ok) throw new HttpError(res.status, `/api/orders: HTTP ${res.status}`);
  return (await res.json()) as {order: string};
}
