// Typed client for the OrthoPan viewer API (docs/api.md).

export interface TrackInfo {
  name: string;
  min: number;
  max: number;
  missing: number;
}

export interface AnnotationInfo {
  name: string;
  kind: 'text' | 'categorical';
  categories?: string[];
}

export interface GenomeMetadataInfo {
  name: string;
  kind: 'categorical' | 'numeric';
  categories?: string[];
  min?: number;
  max?: number;
}

export interface Meta {
  n_genomes: number;
  n_orthogroups: number;
  tree_source: string;
  orders: string[];
  tracks: TrackInfo[];
  annotations?: AnnotationInfo[];
  genome_metadata?: GenomeMetadataInfo[];
  detail?: Record<string, boolean>;
}

export interface Genome {
  row: number;
  genome_idx: number;
  name: string;
  n_orthogroups: number;
  n_genes: number;
  metadata?: Record<string, string | number | null>;
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
  /** echo of the requested value mode */
  value: 'presence' | 'copies';
}

export interface TracksResult {
  data: Float32Array;
  /** track names, in the order the server sent them (X-Tracks header) */
  names: string[];
  bins: number;
}

export interface AnnotationTrackResult {
  name: string;
  data: Uint16Array;
  bins: number;
}

export interface SearchHit {
  orthogroup: string;
  og_idx: number;
  position: number;
  field: string;
  text: string;
}

export interface SearchResult {
  hits: SearchHit[];
  total: number;
}

export interface PatternResult {
  pattern_id: number;
  total: number;
  orthogroups: {orthogroup: string; og_idx: number; annotations: Record<string, unknown>}[];
}

export interface MemberGene {
  id: string;
  raw: string;
}

export interface OrthogroupDetails {
  orthogroup: string;
  og_idx: number;
  n_genomes: number;
  n_genes: number;
  mrca_label: string | null;
  purity: number | null;
  tracks: Record<string, number | null>;
  annotations?: Record<string, unknown>;
  pattern_id?: number;
  n_same_pattern?: number;
  members: {genome: string; genes: MemberGene[]}[];
}

import {decodeFloat32LE, decodeUint16LE} from './decode';
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
  value?: 'presence' | 'copies';
}

export async function getMatrix(p: MatrixParams, signal?: AbortSignal): Promise<MatrixResult> {
  const value = p.value ?? 'presence';
  const q = new URLSearchParams({
    order: p.order,
    og0: String(p.og0),
    og1: String(p.og1),
    g0: String(p.g0),
    g1: String(p.g1),
    w: String(p.w),
    h: String(p.h),
    value,
  });
  const res = await fetch(`/api/matrix?${q}`, {signal});
  if (!res.ok) throw new HttpError(res.status, `/api/matrix: HTTP ${res.status}`);
  const width = Number(res.headers.get('X-Width') ?? p.w);
  const height = Number(res.headers.get('X-Height') ?? p.h);
  const buf = await res.arrayBuffer();
  if (buf.byteLength !== width * height) {
    throw new Error(`matrix payload ${buf.byteLength} != ${width}*${height}`);
  }
  return {data: new Uint8Array(buf), width, height, value};
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
  return {data: decodeFloat32LE(buf), names, bins: p.w};
}

export interface AnnotationTrackParams {
  order: string;
  og0: number;
  og1: number;
  w: number;
  name: string;
}

export async function getAnnotationTrack(
  p: AnnotationTrackParams,
  signal?: AbortSignal,
): Promise<AnnotationTrackResult> {
  const q = new URLSearchParams({
    order: p.order,
    og0: String(p.og0),
    og1: String(p.og1),
    w: String(p.w),
    name: p.name,
  });
  const res = await fetch(`/api/annotation-track?${q}`, {signal});
  if (!res.ok) throw new HttpError(res.status, `/api/annotation-track: HTTP ${res.status}`);
  const bins = p.w;
  const buf = await res.arrayBuffer();
  if (buf.byteLength !== bins * 2) {
    throw new Error(`annotation-track payload ${buf.byteLength} != ${bins}*2`);
  }
  return {name: p.name, data: decodeUint16LE(buf), bins};
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

export function getPattern(og: string, limit = 200, signal?: AbortSignal): Promise<PatternResult> {
  const q = new URLSearchParams({limit: String(limit)});
  return getJSON<PatternResult>(`/api/orthogroup/${encodeURIComponent(og)}/pattern?${q}`, signal);
}

export async function search(q: string, order: string, limit = 50, signal?: AbortSignal): Promise<SearchResult> {
  const params = new URLSearchParams({q, order, limit: String(limit)});
  return getJSON<SearchResult>(`/api/search?${params}`, signal);
}

export async function postOrder(body: {
  by: string;
  genomes?: string[];
  q?: string;
  og?: string;
}): Promise<{order: string}> {
  const res = await fetch('/api/orders', {
    method: 'POST',
    headers: {'content-type': 'application/json'},
    body: JSON.stringify(body),
  });
  if (!res.ok) throw new HttpError(res.status, `/api/orders: HTTP ${res.status}`);
  return (await res.json()) as {order: string};
}

// --- v2: single-orthogroup detail view -------------------------------------

export interface MsaRow {
  genome: string;
  row: number;
  gene: string;
  raw: string;
  seq: string;
}

export interface MsaResult {
  orthogroup: string;
  source: 'orthofinder' | 'famsa';
  length: number;
  trimmed: boolean | null;
  rows: MsaRow[];
}

export interface GeneTreeArrays {
  orthogroup: string;
  parent: number[];
  length: number[];
  label: string[];
  gene: (string | null)[];
  genome_row: number[];
  duplication: boolean[];
  support: (number | null)[];
}

export type IdMode = 'display' | 'raw';

/** Read a JSON error body and surface its `detail` as the error message. */
async function getJSONDetail<T>(path: string, signal?: AbortSignal): Promise<T> {
  const res = await fetch(BASE + path, {signal});
  if (!res.ok) {
    let detail = `${path}: HTTP ${res.status}`;
    try {
      const body = (await res.json()) as {detail?: unknown};
      if (body && body.detail !== undefined) detail = String(body.detail);
    } catch {
      /* not JSON */
    }
    throw new HttpError(res.status, detail);
  }
  return (await res.json()) as T;
}

export function getMsa(og: string, signal?: AbortSignal): Promise<MsaResult> {
  return getJSONDetail<MsaResult>(`/api/orthogroup/${encodeURIComponent(og)}/msa`, signal);
}

export function getGeneTree(og: string, signal?: AbortSignal): Promise<GeneTreeArrays> {
  return getJSONDetail<GeneTreeArrays>(`/api/orthogroup/${encodeURIComponent(og)}/gene-tree`, signal);
}

export function fastaUrl(og: string, aligned: boolean, ids: IdMode): string {
  return `/api/orthogroup/${encodeURIComponent(og)}/fasta?aligned=${aligned ? 1 : 0}&ids=${ids}`;
}

export function geneTreeNwkUrl(og: string, ids: IdMode): string {
  return `/api/orthogroup/${encodeURIComponent(og)}/gene-tree.nwk?ids=${ids}`;
}
