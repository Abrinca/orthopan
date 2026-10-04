// Synthetic OrthoPan API server for frontend development without the backend.
// Implements docs/api.md on 50 genomes x 2000 orthogroups with clade structure.
// Also serves web/dist so `bun run build` output can be checked headlessly.
//
//   bun run web/mock/server.ts            # port 8000, serves dist if present
//   bun run web/mock/server.ts --port 9000
import {existsSync} from 'node:fs';
import {join} from 'node:path';

const args = process.argv.slice(2);
const PORT = Number(args.indexOf('--port') >= 0 ? args[args.indexOf('--port') + 1] : 8000);
const N_G = 50;
const N_OG = 2000;
const CLADES = 5;
const CLADE_SIZE = N_G / CLADES;
const DIST = join(import.meta.dir, '..', 'dist');
const MAX_CELLS = 4096;
const MAX_CELLS_TOTAL = 4_000_000;
const MAX_LIST = 5000;
const MAX_BODY = 64 * 1024;

// --- deterministic PRNG -----------------------------------------------------
function mulberry32(seed: number) {
  let a = seed >>> 0;
  return () => {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const rng = mulberry32(42);

// --- synthetic tree ---------------------------------------------------------
interface Node {
  id: number;
  label: string;
  length: number;
  children: number[];
  parent: number;
  row: number;
  lo: number;
  hi: number;
  depth: number;
  clade: number;
}

const nodes: Node[] = [];
let nodeId = 0;
function mkNode(label: string, clade: number): Node {
  const n: Node = {id: nodeId++, label, length: 0.01 + rng() * 0.3, children: [], parent: -1, row: -1, lo: 0, hi: 0, depth: 0, clade};
  nodes.push(n);
  return n;
}
const subId: number[] = [];

function buildRange(rows: number[], clade: number, depth: number): number {
  if (rows.length === 1) {
    const leaf = mkNode(`GENOME_${String(rows[0]).padStart(3, '0')}`, clade);
    return leaf.id;
  }
  const label = depth === 0 ? `clade${clade}` : `N${subId[clade]++}`;
  const node = mkNode(label, clade);
  const mid = Math.ceil(rows.length / 2);
  for (const half of [rows.slice(0, mid), rows.slice(mid)]) {
    const child = buildRange(half, clade, depth + 1);
    nodes[child].parent = node.id;
    node.children.push(child);
  }
  return node.id;
}

for (let c = 0; c < CLADES; c++) {
  subId.push(1);
  const rows = Array.from({length: CLADE_SIZE}, (_, i) => c * CLADE_SIZE + i);
  buildRange(rows, c, 0);
}

// root spanning all clades
const root = mkNode('root', -1);
root.children = nodes.filter((n) => n.parent === -1 && n.id !== root.id).map((n) => n.id);
for (const c of root.children) nodes[c].parent = root.id;

// pre-order arrays, leaf rows in DFS order
const parentArr: number[] = [];
const lengthArr: number[] = [];
const labelArr: string[] = [];
const rowArr: number[] = [];
let leaf = 0;
const preorder: number[] = [];
const preIdx = new Map<number, number>();
const stack = [root.id];
while (stack.length) {
  const id = stack.pop()!;
  preIdx.set(id, preorder.length);
  preorder.push(id);
  const node = nodes[id];
  const pi = node.parent < 0 ? -1 : (preIdx.get(node.parent) ?? -1);
  parentArr.push(pi);
  lengthArr.push(pi < 0 ? 0 : node.length);
  labelArr.push(node.label);
  if (node.children.length === 0) {
    node.row = leaf++;
    rowArr.push(node.row);
  } else {
    rowArr.push(-1);
    for (let k = node.children.length - 1; k >= 0; k--) stack.push(node.children[k]);
  }
}
// depth / lo / hi in preorder index space
const idx = preIdx;
const childrenIdx: number[][] = preorder.map((id) => nodes[id].children.map((c) => idx.get(c)!));
const depthArr: number[] = new Array(preorder.length).fill(0);
for (let i = 0; i < preorder.length; i++) {
  for (const c of childrenIdx[i]) depthArr[c] = depthArr[i] + lengthArr[c];
}
const loArr: number[] = new Array(preorder.length).fill(0);
const hiArr: number[] = new Array(preorder.length).fill(0);
for (let i = preorder.length - 1; i >= 0; i--) {
  if (childrenIdx[i].length === 0) {
    loArr[i] = hiArr[i] = rowArr[i];
  } else {
    loArr[i] = loArr[childrenIdx[i][0]];
    hiArr[i] = hiArr[childrenIdx[i][childrenIdx[i].length - 1]];
  }
}

const genomeClade: number[] = new Array(N_G);
preorder.forEach((id, i) => {
  if (nodes[id].children.length === 0) genomeClade[rowArr[i]] = nodes[id].clade;
});
const cladeRootIdx: number[] = [];
for (let c = 0; c < CLADES; c++) {
  const node = nodes.find((n) => n.label === `clade${c}`);
  cladeRootIdx.push(node ? idx.get(node.id)! : 0);
}
const genomeNames: string[] = new Array(N_G);
preorder.forEach((id, i) => {
  if (rowArr[i] >= 0) genomeNames[rowArr[i]] = nodes[id].label;
});
// node whose range does not cross a clade boundary (mrca surrogate)
function mrcaIndex(minRow: number, maxRow: number): number {
  let best = 0;
  for (let i = 0; i < preorder.length; i++) {
    if (loArr[i] <= minRow && hiArr[i] >= maxRow && depthArr[i] >= depthArr[best]) best = i;
  }
  return best;
}
function cladeOfRow(r: number): number {
  return genomeClade[r];
}

// --- presence + tracks ------------------------------------------------------
const presence = new Uint8Array(N_G * N_OG); // 0/1
const copies = new Uint8Array(N_G * N_OG); // 0..3
const geom = new Float64Array(N_OG);
const func = new Float64Array(N_OG);
const cnMean = new Float64Array(N_OG);
const ogNGenomes = new Int32Array(N_OG);
const ogNGenes = new Int32Array(N_OG);
const ogMrca: number[] = new Array(N_OG);
const ogPurity = new Float64Array(N_OG);
const ogBlockKey = new Float64Array(N_OG);

for (let og = 0; og < N_OG; og++) {
  const cladeRestricted = og < CLADES * 300;
  const clade = cladeRestricted ? Math.floor(og / 300) : -1;
  let minR = N_G;
  let maxR = -1;
  let nG = 0;
  let nGenes = 0;
  for (let r = 0; r < N_G; r++) {
    const p = cladeRestricted
      ? cladeOfRow(r) === clade
        ? 0.9
        : 0.05
      : 0.15;
    if (rng() < p) {
      const c = rng() < 0.15 ? 2 : 1;
      presence[r * N_OG + og] = 1;
      copies[r * N_OG + og] = c;
      nG++;
      nGenes += c;
      if (r < minR) minR = r;
      if (r > maxR) maxR = r;
    }
  }
  ogNGenomes[og] = nG;
  ogNGenes[og] = nGenes;
  cnMean[og] = nG ? nGenes / nG : NaN;
  const m = cladeRestricted ? cladeRootIdx[clade] : mrcaIndex(minR, maxR);
  ogMrca[og] = m;
  ogPurity[og] = nG ? nG / (hiArr[m] - loArr[m] + 1) : NaN;
  ogBlockKey[og] = cladeRestricted ? cladeRootIdx[clade] : 1e9;
  geom[og] = rng() < 0.1 ? NaN : (cladeRestricted ? 0.45 : 0.25) + 0.6 * rng();
  func[og] = rng() < 0.1 ? NaN : (cladeRestricted ? 0.4 : 0.2) + 0.55 * rng();
}

const trackData: Record<string, Float64Array> = {
  copy_number_mean: cnMean,
  geometric_homogeneity: geom,
  functional_homogeneity: func,
};

function trackInfo(name: string) {
  const arr = trackData[name];
  let min = Infinity;
  let max = -Infinity;
  let missing = 0;
  for (const v of arr) {
    if (Number.isNaN(v)) missing++;
    else {
      if (v < min) min = v;
      if (v > max) max = v;
    }
  }
  return {name, min: Number.isFinite(min) ? min : 0, max: Number.isFinite(max) ? max : 1, missing};
}

// --- orderings --------------------------------------------------------------
function ranks(n: number, key: (i: number) => number | string): number[] {
  const perm = Array.from({length: n}, (_, i) => i).sort((a, b) => {
    const ka = key(a);
    const kb = key(b);
    if (ka < kb) return -1;
    if (ka > kb) return 1;
    return a - b;
  });
  return perm; // perm[pos] = og index
}
const orderTree = ranks(N_OG, (i) => ogBlockKey[i] * 1e6 + (1 - ogPurity[i]) * 1000 + (N_OG - ogNGenomes[i]));
const orderFreq = ranks(N_OG, (i) => N_OG - ogNGenomes[i]);
const orders: Record<string, number[]> = {tree: orderTree, frequency: orderFreq};
let dynamicId = 0;

function trackPerm(spec: string): number[] | null {
  const desc = spec.endsWith(':desc');
  const name = desc ? spec.slice('track:'.length, -':desc'.length) : spec.slice('track:'.length);
  if (!(name in trackData)) return null;
  const arr = trackData[name];
  const present: number[] = [];
  const missing: number[] = [];
  for (let i = 0; i < N_OG; i++) (Number.isNaN(arr[i]) ? missing : present).push(i);
  present.sort((a, b) => (desc ? arr[b] - arr[a] : arr[a] - arr[b]) || a - b);
  return [...present, ...missing]; // missing values last
}

/** Resolve an order id to a position -> og-index permutation, or null if unknown. */
function resolveOrder(order: string): number[] | null {
  if (order in orders) return orders[order];
  if (order.startsWith('track:')) return trackPerm(order);
  return null;
}

// --- API helpers ------------------------------------------------------------
function json(obj: unknown): Response {
  return new Response(JSON.stringify(obj), {headers: {'content-type': 'application/json'}});
}
function bad(status: number, detail: string): Response {
  return new Response(JSON.stringify({detail}), {status, headers: {'content-type': 'application/json'}});
}
function intParam(q: URLSearchParams, key: string, def: number): number {
  const v = q.get(key);
  return v === null ? def : Number(v);
}
function validateWindow(og0: number, og1: number, g0: number, g1: number, w: number, h: number): string | null {
  const ints = [og0, og1, g0, g1, w, h];
  if (!ints.every((x) => Number.isInteger(x))) return 'window out of range';
  if (!(0 <= og0 && og0 < og1 && og1 <= N_OG) || !(0 <= g0 && g0 < g1 && g1 <= N_G)) {
    return 'window out of range';
  }
  if (!(1 <= w && w <= MAX_CELLS) || !(1 <= h && h <= MAX_CELLS)) return `w and h must be in 1..${MAX_CELLS}`;
  if (w > og1 - og0 || h > g1 - g0) return 'w/h exceed window size';
  if (w * h > MAX_CELLS_TOTAL) return `w*h exceeds ${MAX_CELLS_TOTAL} output cells`;
  return null;
}
function matrixResponse(perm: number[], og0: number, og1: number, g0: number, g1: number, w: number, h: number): Response {
  const L = og1 - og0;
  const M = g1 - g0;
  const outW = Math.min(w, L);
  const outH = Math.min(h, M);
  const out = new Uint8Array(outW * outH);
  for (let j = 0; j < outH; j++) {
    const r0 = g0 + Math.floor((j * M) / outH);
    const r1 = g0 + Math.floor(((j + 1) * M) / outH);
    for (let i = 0; i < outW; i++) {
      const p0 = og0 + Math.floor((i * L) / outW);
      const p1 = og0 + Math.floor(((i + 1) * L) / outW);
      let present = 0;
      let total = 0;
      for (let r = r0; r < r1; r++) {
        const base = r * N_OG;
        for (let p = p0; p < p1; p++) {
          total++;
          present += presence[base + perm[p]];
        }
      }
      out[j * outW + i] = total ? Math.round((present / total) * 255) : 0;
    }
  }
  return new Response(out, {
    headers: {
      'content-type': 'application/octet-stream',
      'X-Width': String(outW),
      'X-Height': String(outH),
    },
  });
}

function tracksResponse(perm: number[], og0: number, og1: number, w: number, names: string[]): Response {
  const L = og1 - og0;
  const outW = Math.min(w, L);
  const data = new Float32Array(names.length * outW);
  names.forEach((name, t) => {
    const arr = trackData[name];
    for (let i = 0; i < outW; i++) {
      const p0 = og0 + Math.floor((i * L) / outW);
      const p1 = og0 + Math.floor(((i + 1) * L) / outW);
      let sum = 0;
      let n = 0;
      for (let p = p0; p < p1; p++) {
        const v = arr[perm[p]];
        if (!Number.isNaN(v)) {
          sum += v;
          n++;
        }
      }
      data[t * outW + i] = n ? sum / n : NaN;
    }
  });
  return new Response(data, {
    headers: {
      'content-type': 'application/octet-stream',
      'X-Tracks': names.join(','),
    },
  });
}

function orthogroupDetails(ogIdx: number): Response {
  const og = `OG${String(ogIdx).padStart(7, '0')}`;
  const members: {genome: string; genes: string[]}[] = [];
  for (let r = 0; r < N_G; r++) {
    if (presence[r * N_OG + ogIdx]) {
      const c = copies[r * N_OG + ogIdx];
      members.push({
        genome: genomeNames[r],
        genes: Array.from({length: c}, (_, k) => `gnl|${genomeNames[r]}_${String(ogIdx).padStart(6, '0')}_${k + 1}`),
      });
    }
  }
  return json({
    orthogroup: og,
    og_idx: ogIdx,
    n_genomes: ogNGenomes[ogIdx],
    n_genes: ogNGenes[ogIdx],
    mrca_label: labelArr[ogMrca[ogIdx]],
    purity: ogPurity[ogIdx],
    tracks: {
      copy_number_mean: Number.isNaN(cnMean[ogIdx]) ? null : cnMean[ogIdx],
      geometric_homogeneity: Number.isNaN(geom[ogIdx]) ? null : geom[ogIdx],
    },
    members,
  });
}

// --- routing ----------------------------------------------------------------
const server = Bun.serve({
  port: PORT,
  async fetch(req) {
    const url = new URL(req.url);
    const p = url.pathname;
    if (!p.startsWith('/api/')) {
      return serveStatic(p);
    }
    const q = url.searchParams;
    if (p === '/api/meta') {
      return json({
        n_genomes: N_G,
        n_orthogroups: N_OG,
        tree_source: 'mock',
        orders: Object.keys(orders),
        tracks: Object.keys(trackData).map(trackInfo),
      });
    }
    if (p === '/api/genomes') {
      const list = genomeNames.map((name, row) => {
        let nOG = 0;
        let nGenes = 0;
        for (let og = 0; og < N_OG; og++) {
          const c = copies[row * N_OG + og];
          if (c) {
            nOG++;
            nGenes += c;
          }
        }
        return {row, genome_idx: row, name, n_orthogroups: nOG, n_genes: nGenes};
      });
      return json(list);
    }
    if (p === '/api/tree') {
      return json({parent: parentArr, length: lengthArr, label: labelArr, row: rowArr});
    }
    if (p === '/api/matrix') {
      const order = q.get('order') ?? 'tree';
      const perm = resolveOrder(order);
      if (!perm) return bad(404, `unknown order: ${order}`);
      const og0 = intParam(q, 'og0', 0);
      const og1 = intParam(q, 'og1', N_OG);
      const g0 = intParam(q, 'g0', 0);
      const g1 = intParam(q, 'g1', N_G);
      const w = intParam(q, 'w', 0);
      const h = intParam(q, 'h', 0);
      const e = validateWindow(og0, og1, g0, g1, w, h);
      if (e) return bad(400, e);
      return matrixResponse(perm, og0, og1, g0, g1, w, h);
    }
    if (p === '/api/tracks') {
      const order = q.get('order') ?? 'tree';
      const perm = resolveOrder(order);
      if (!perm) return bad(404, `unknown order: ${order}`);
      const requested = (q.get('names') ?? '').split(',').filter(Boolean);
      const names = requested.length ? requested : Object.keys(trackData);
      for (const n of names) if (!(n in trackData)) return bad(404, `unknown track: ${n}`);
      const og0 = intParam(q, 'og0', 0);
      const og1 = intParam(q, 'og1', N_OG);
      const w = intParam(q, 'w', 0);
      if (!(Number.isInteger(og0) && Number.isInteger(og1) && 0 <= og0 && og0 < og1 && og1 <= N_OG)) {
        return bad(400, 'window out of range');
      }
      if (!(Number.isInteger(w) && 1 <= w && w <= MAX_CELLS)) return bad(400, `w must be in 1..${MAX_CELLS}`);
      if (w > og1 - og0) return bad(400, 'w exceeds window size');
      return tracksResponse(perm, og0, og1, w, names);
    }
    if (p === '/api/orthogroups') {
      const order = q.get('order') ?? 'tree';
      const perm = resolveOrder(order);
      if (!perm) return bad(404, `unknown order: ${order}`);
      const og0 = intParam(q, 'og0', 0);
      const og1 = intParam(q, 'og1', N_OG);
      if (!(Number.isInteger(og0) && Number.isInteger(og1) && 0 <= og0 && og0 < og1 && og1 <= N_OG)) {
        return bad(400, 'window out of range');
      }
      if (og1 - og0 > MAX_LIST) return bad(400, `at most ${MAX_LIST} per request`);
      const ids = [];
      for (let pos = og0; pos < og1; pos++) ids.push(`OG${String(perm[pos]).padStart(7, '0')}`);
      return json(ids);
    }
    if (p.startsWith('/api/orthogroup/')) {
      const id = decodeURIComponent(p.slice('/api/orthogroup/'.length));
      const m = /^OG(\d+)$/.exec(id);
      if (!m || Number(m[1]) >= N_OG) return bad(404, `unknown orthogroup: ${id}`);
      return orthogroupDetails(Number(m[1]));
    }
    if (p === '/api/orders' && req.method === 'POST') {
      const length = req.headers.get('content-length');
      if (length !== null && length !== '' && Number(length) > MAX_BODY) return bad(413, 'request body too large');
      let raw: string;
      try {
        raw = await req.text();
      } catch {
        return bad(400, 'invalid JSON body');
      }
      if (raw.length > MAX_BODY) return bad(413, 'request body too large');
      let body: {by?: unknown; genomes?: unknown};
      try {
        body = JSON.parse(raw) as {by?: unknown; genomes?: unknown};
      } catch {
        return bad(400, 'invalid JSON body');
      }
      if (!body || typeof body !== 'object' || body.by !== 'contrast') {
        return bad(400, "only 'contrast' orders are supported");
      }
      const names = body.genomes;
      if (!Array.isArray(names) || names.length === 0 || !names.every((n) => typeof n === 'string')) {
        return bad(400, "'genomes' must be a non-empty list of strings");
      }
      const unique = [...new Set(names as string[])];
      if (unique.length > N_G) return bad(400, 'more genomes than exist');
      for (const n of unique) if (!genomeNames.includes(n)) return bad(400, `unknown genome: ${n}`);
      if (unique.length === N_G) return bad(400, 'contrast needs 1..n-1 genomes');
      const insideRows = new Set<number>();
      genomeNames.forEach((n, r) => {
        if (unique.includes(n)) insideRows.add(r);
      });
      const score = new Float64Array(N_OG);
      for (let og = 0; og < N_OG; og++) {
        let inP = 0;
        let outP = 0;
        let inN = 0;
        let outN = 0;
        for (let r = 0; r < N_G; r++) {
          const v = presence[r * N_OG + og];
          if (insideRows.has(r)) {
            inP += v;
            inN++;
          } else {
            outP += v;
            outN++;
          }
        }
        score[og] = (inN ? inP / inN : 0) - (outN ? outP / outN : 0);
      }
      const perm = ranks(N_OG, (i) => -score[i]);
      const id = `contrast-${++dynamicId}`;
      orders[id] = perm;
      return json({order: id});
    }
    return bad(404, 'not found');
  },
});

async function serveStatic(pathname: string): Promise<Response> {
  let file = join(DIST, pathname === '/' ? 'index.html' : pathname);
  if (!existsSync(file)) {
    if (pathname.includes('.')) return new Response('not found', {status: 404});
    file = join(DIST, 'index.html');
  }
  if (!existsSync(file)) return new Response('frontend not built (run: bun run build)', {status: 404});
  return new Response(Bun.file(file));
}

console.log(`mock OrthoPan API on http://127.0.0.1:${server.port} (dist: ${existsSync(join(DIST, 'index.html')) ? 'yes' : 'no'})`);
