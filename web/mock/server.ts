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

// --- v1: annotations, genome metadata and presence patterns ----------------

const COG_NAMES = ['J', 'K', 'L', 'D', 'O', 'M', 'G', 'E', 'F', 'H'];
const ogDescription: string[] = new Array(N_OG);
const ogCog: string[] = new Array(N_OG);
const cogCount = new Map<string, number>();
const DESC_WORDS = ['transporter', 'kinase', 'synthase', 'regulator', 'hydrolase', 'binding', 'ribosomal', 'oxidase', 'reductase', 'permease'];
for (let og = 0; og < N_OG; og++) {
  ogDescription[og] = `${DESC_WORDS[og % DESC_WORDS.length]} family OG${String(og).padStart(7, '0')}`;
  const cog = COG_NAMES[Math.floor(rng() * COG_NAMES.length)];
  ogCog[og] = cog;
  cogCount.set(cog, (cogCount.get(cog) ?? 0) + 1);
}
// categories ordered by frequency, most common first (ties keep the fixed order)
const cogCategories = [...COG_NAMES].sort(
  (a, b) => (cogCount.get(b)! - cogCount.get(a)!) || COG_NAMES.indexOf(a) - COG_NAMES.indexOf(b),
);
const cogIndex = new Map(cogCategories.map((c, i) => [c, i]));

const genomeCladeName: string[] = new Array(N_G);
const genomeGc: number[] = new Array(N_G);
for (let r = 0; r < N_G; r++) {
  genomeCladeName[r] = `clade${genomeClade[r]}`;
  genomeGc[r] = Number((0.35 + 0.06 * rng()).toFixed(4));
}
const metaAnnotations = [
  {name: 'description', kind: 'text' as const},
  {name: 'cog_category', kind: 'categorical' as const, categories: cogCategories},
];
const metaGenomeMetadata = [
  {name: 'clade', kind: 'categorical' as const, categories: [...new Set(genomeCladeName)]},
  {name: 'gc', kind: 'numeric' as const, min: 0.35, max: 0.41},
];

// pattern ids: identical presence columns share one
const ogPatternId = new Int32Array(N_OG);
const patternGroups = new Map<string, number[]>();
function ogSignature(og: number): string {
  let s = '';
  for (let r = 0; r < N_G; r++) s += presence[r * N_OG + og];
  return s;
}
for (let og = 0; og < N_OG; og++) {
  const sig = ogSignature(og);
  let members = patternGroups.get(sig);
  if (!members) {
    members = [];
    patternGroups.set(sig, members);
    ogPatternId[og] = patternGroups.size - 1;
  } else {
    ogPatternId[og] = ogPatternId[members[0]];
  }
  members.push(og);
}
const patternSize = (og: number) => patternGroups.get(ogSignature(og))!.length;

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
function matrixResponse(
  perm: number[],
  og0: number,
  og1: number,
  g0: number,
  g1: number,
  w: number,
  h: number,
  value: 'presence' | 'copies',
): Response {
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
      let copySum = 0;
      for (let r = r0; r < r1; r++) {
        const base = r * N_OG;
        for (let p = p0; p < p1; p++) {
          const og = perm[p];
          total++;
          present += presence[base + og];
          copySum += copies[base + og];
        }
      }
      if (value === 'copies') {
        out[j * outW + i] = present ? Math.min(255, Math.round((10 * copySum) / present)) : 0;
      } else {
        out[j * outW + i] = total ? Math.round((present / total) * 255) : 0;
      }
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

function annotationTrackResponse(
  perm: number[],
  og0: number,
  og1: number,
  w: number,
  name: string,
): Response {
  const L = og1 - og0;
  const outW = Math.min(w, L);
  const out = new Uint16Array(outW);
  for (let i = 0; i < outW; i++) {
    const p0 = og0 + Math.floor((i * L) / outW);
    const p1 = og0 + Math.floor(((i + 1) * L) / outW);
    const counts = new Map<string, number>();
    for (let p = p0; p < p1; p++) {
      const c = ogCog[perm[p]];
      counts.set(c, (counts.get(c) ?? 0) + 1);
    }
    let best = '';
    let bestN = -1;
    for (const c of cogCategories) {
      const n = counts.get(c) ?? 0;
      if (n > bestN) {
        bestN = n;
        best = c;
      }
    }
    out[i] = bestN > 0 ? 1 + cogIndex.get(best)! : 0;
  }
  return new Response(out, {
    headers: {'content-type': 'application/octet-stream', 'X-Width': String(outW)},
  });
}

function computeHits(term: string, perm: number[]): {orthogroup: string; og_idx: number; position: number; field: string; text: string}[] {
  const posOf = new Int32Array(N_OG);
  perm.forEach((og, pos) => (posOf[og] = pos));
  const needle = term.toLowerCase();
  const hits: {orthogroup: string; og_idx: number; position: number; field: string; text: string}[] = [];
  for (let og = 0; og < N_OG; og++) {
    const id = `OG${String(og).padStart(7, '0')}`;
    let field = '';
    let text = '';
    if (id.toLowerCase().includes(needle)) {
      field = 'orthogroup';
      text = id;
    } else if (ogDescription[og].toLowerCase().includes(needle)) {
      field = 'description';
      text = ogDescription[og];
    } else if (ogCog[og].toLowerCase().includes(needle)) {
      field = 'cog_category';
      text = ogCog[og];
    }
    if (field) hits.push({orthogroup: id, og_idx: og, position: posOf[og], field, text});
  }
  hits.sort((a, b) => a.position - b.position);
  return hits;
}

function searchResponse(q: string, order: string, limit: number): Response {
  const perm = resolveOrder(order);
  if (!perm) return bad(404, `unknown order: ${order}`);
  const hits = computeHits(q, perm);
  return json({hits: hits.slice(0, limit), total: hits.length});
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
  const members: {genome: string; genes: {id: string; raw: string}[]}[] = [];
  for (let r = 0; r < N_G; r++) {
    if (presence[r * N_OG + ogIdx]) {
      const c = copies[r * N_OG + ogIdx];
      members.push({
        genome: genomeNames[r],
        genes: Array.from({length: c}, (_, k) => {
          const raw = `gnl|extdb|${genomeNames[r]}_${String(ogIdx).padStart(6, '0')}_${k + 1}`;
          return {id: raw.replace(/^gnl\|[^|]+\|/, ''), raw};
        }),
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
    annotations: {description: ogDescription[ogIdx], cog_category: ogCog[ogIdx]},
    pattern_id: ogPatternId[ogIdx],
    n_same_pattern: patternSize(ogIdx),
    members,
  });
}

function patternResponse(ogIdx: number): Response {
  const members = patternGroups.get(ogSignature(ogIdx))!;
  return json({
    pattern_id: ogPatternId[ogIdx],
    total: members.length,
    orthogroups: members.map((og) => ({
      orthogroup: `OG${String(og).padStart(7, '0')}`,
      og_idx: og,
      annotations: {description: ogDescription[og], cog_category: ogCog[og]},
    })),
  });
}

// --- v2: MSA / gene tree / downloads ---------------------------------------

const AMINO = 'ACDEFGHIKLMNPQRSTVWY';
const msaBusyAt = new Map<number, number>();

interface GeneRef {
  genome: string;
  row: number;
  gene: string;
  raw: string;
}

function genesOf(ogIdx: number): GeneRef[] {
  const out: GeneRef[] = [];
  for (let r = 0; r < N_G; r++) {
    if (!presence[r * N_OG + ogIdx]) continue;
    const c = copies[r * N_OG + ogIdx];
    for (let k = 0; k < c; k++) {
      const raw = `gnl|extdb|${genomeNames[r]}_${String(ogIdx).padStart(6, '0')}_${k + 1}`;
      out.push({genome: genomeNames[r], row: r, gene: raw.replace(/^gnl\|[^|]+\|/, ''), raw});
    }
  }
  out.sort((a, b) => a.row - b.row || (a.gene < b.gene ? -1 : 1));
  return out;
}

const hasMsa = (og: number) => og % 97 !== 42; // some OGs have no alignment
const msaSource = (og: number): 'orthofinder' | 'famsa' => (og % 11 === 3 ? 'famsa' : 'orthofinder');

function makeMsa(ogIdx: number) {
  const genes = genesOf(ogIdx);
  const length = 24 + (ogIdx % 40);
  const rows = genes.map((g, gi) => {
    let seq = '';
    for (let c = 0; c < length; c++) {
      const base = AMINO[(c * 7 + genomeClade[g.row] * 3) % AMINO.length];
      const m = (c * 13 + gi * 5) % 23;
      seq += m === 0 ? AMINO[(c * 3 + gi) % AMINO.length] : m === 1 ? '-' : base;
    }
    return {genome: g.genome, row: g.row, gene: g.gene, raw: g.raw, seq};
  });
  return {orthogroup: `OG${String(ogIdx).padStart(7, '0')}`, source: msaSource(ogIdx), length, trimmed: ogIdx % 3 === 0 ? true : ogIdx % 3 === 1 ? false : null, rows};
}

const hasTree = (og: number) => genesOf(og).length >= 4 && og % 101 !== 42;

function makeGeneTree(ogIdx: number) {
  const genes = genesOf(ogIdx);
  const parent: number[] = [];
  const length: number[] = [];
  const label: string[] = [];
  const gene: (string | null)[] = [];
  const genomeRow: number[] = [];
  const duplication: boolean[] = [];
  const support: (number | null)[] = [];
  let id = 0;
  const build = (lo: number, hi: number, par: number): void => {
    const my = id++;
    parent[my] = par;
    length[my] = 0.05 + ((my * 7) % 10) / 100;
    support[my] = par < 0 ? null : Number((0.5 + (my % 5) / 10).toFixed(2));
    if (hi - lo === 1) {
      const g = genes[lo];
      label[my] = g.gene;
      gene[my] = g.gene;
      genomeRow[my] = g.row;
      duplication[my] = false;
    } else {
      label[my] = `n${my}`;
      gene[my] = null;
      genomeRow[my] = -1;
      duplication[my] = my % 3 === 0;
      const mid = Math.floor((lo + hi) / 2);
      build(lo, mid, my);
      build(mid, hi, my);
    }
  };
  build(0, genes.length, -1);
  return {orthogroup: `OG${String(ogIdx).padStart(7, '0')}`, parent, length, label, gene, genome_row: genomeRow, duplication, support};
}

function newickOf(tree: ReturnType<typeof makeGeneTree>, ids: 'display' | 'raw', ogIdx: number): string {
  const children: number[][] = tree.parent.map(() => []);
  tree.parent.forEach((p, i) => {
    if (p >= 0) children[p].push(i);
  });
  const genes = new Map(genesOf(ogIdx).map((g) => [g.gene, g]));
  const name = (i: number) => {
    const leafGene = tree.gene[i];
    if (leafGene) return (ids === 'raw' ? genes.get(leafGene)?.raw : leafGene) ?? leafGene;
    return tree.label[i];
  };
  const rec = (i: number): string => {
    const kids = children[i];
    const bl = tree.length[i] ? `:${tree.length[i].toFixed(4)}` : '';
    if (!kids.length) return `${name(i)}${bl}`;
    return `(${kids.map(rec).join(',')})${name(i)}${bl}`;
  };
  return rec(0) + ';';
}

function download(content: string, filename: string): Response {
  return new Response(content, {
    headers: {
      'content-type': 'text/plain; charset=utf-8',
      'content-disposition': `attachment; filename="${filename}"`,
    },
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
        annotations: metaAnnotations,
        genome_metadata: metaGenomeMetadata,
        detail: {msa: true, gene_trees: true, sequences: true, famsa: true},
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
        return {
          row,
          genome_idx: row,
          name,
          n_orthogroups: nOG,
          n_genes: nGenes,
          metadata: {clade: genomeCladeName[row], gc: genomeGc[row]},
        };
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
      const value = q.get('value') === 'copies' ? 'copies' : 'presence';
      const og0 = intParam(q, 'og0', 0);
      const og1 = intParam(q, 'og1', N_OG);
      const g0 = intParam(q, 'g0', 0);
      const g1 = intParam(q, 'g1', N_G);
      const w = intParam(q, 'w', 0);
      const h = intParam(q, 'h', 0);
      const e = validateWindow(og0, og1, g0, g1, w, h);
      if (e) return bad(400, e);
      return matrixResponse(perm, og0, og1, g0, g1, w, h, value);
    }
    if (p === '/api/annotation-track') {
      const order = q.get('order') ?? 'tree';
      const perm = resolveOrder(order);
      if (!perm) return bad(404, `unknown order: ${order}`);
      const name = q.get('name') ?? '';
      const info = metaAnnotations.find((a) => a.name === name);
      if (!info) return bad(404, `unknown annotation: ${name}`);
      if (info.kind !== 'categorical') return bad(400, `not a categorical annotation: ${name}`);
      const og0 = intParam(q, 'og0', 0);
      const og1 = intParam(q, 'og1', N_OG);
      const w = intParam(q, 'w', 0);
      if (!(Number.isInteger(og0) && Number.isInteger(og1) && 0 <= og0 && og0 < og1 && og1 <= N_OG)) {
        return bad(400, 'window out of range');
      }
      if (!(Number.isInteger(w) && 1 <= w && w <= MAX_CELLS)) return bad(400, `w must be in 1..${MAX_CELLS}`);
      if (w > og1 - og0) return bad(400, 'w exceeds window size');
      return annotationTrackResponse(perm, og0, og1, w, name);
    }
    if (p === '/api/search') {
      const term = (q.get('q') ?? '').trim();
      if (!term) return bad(400, 'q is required');
      const order = q.get('order') ?? 'tree';
      const limit = Math.min(500, Math.max(1, intParam(q, 'limit', 50)));
      return searchResponse(term, order, limit);
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
      const rest = p.slice('/api/orthogroup/'.length);
      const suffix = ['/pattern', '/msa', '/gene-tree.nwk', '/gene-tree', '/fasta'].find((s) => rest.endsWith(s));
      const id = decodeURIComponent(suffix ? rest.slice(0, -suffix.length) : rest);
      const m = /^OG(\d+)$/.exec(id);
      if (!m || Number(m[1]) >= N_OG) return bad(404, `unknown orthogroup: ${id}`);
      const ogIdx = Number(m[1]);
      const ids: 'display' | 'raw' = q.get('ids') === 'raw' ? 'raw' : 'display';
      if (suffix === '/msa') {
        if (!hasMsa(ogIdx)) return bad(404, 'no MSA for this orthogroup (no sequences available)');
        if (ogIdx === 1500 && (msaBusyAt.get(ogIdx) ?? 0) < Date.now() - 3000) {
          msaBusyAt.set(ogIdx, Date.now());
          return bad(503, 'alignment limit busy, try again');
        }
        if (msaSource(ogIdx) === 'famsa') await Bun.sleep(250);
        return json(makeMsa(ogIdx));
      }
      if (suffix === '/gene-tree') {
        if (!hasTree(ogIdx)) return bad(404, 'no gene tree for this orthogroup (<4 genes)');
        return json(makeGeneTree(ogIdx));
      }
      if (suffix === '/gene-tree.nwk') {
        if (!hasTree(ogIdx)) return bad(404, 'no gene tree for this orthogroup (<4 genes)');
        return download(newickOf(makeGeneTree(ogIdx), ids, ogIdx), `OG${String(ogIdx).padStart(7, '0')}.nwk`);
      }
      if (suffix === '/fasta') {
        const aligned = q.get('aligned') === '1';
        if (aligned && !hasMsa(ogIdx)) return bad(404, 'no MSA for this orthogroup');
        const rows = aligned
          ? makeMsa(ogIdx).rows
          : genesOf(ogIdx).map((g) => ({...g, seq: ''}));
        const body = rows
          .map((r) => `>${ids === 'raw' ? (r as GeneRef).raw : (r as GeneRef).gene}\n${aligned ? (r as {seq: string}).seq : 'MKT'}`)
          .join('\n');
        const name = `OG${String(ogIdx).padStart(7, '0')}${aligned ? '.aligned' : ''}.fa`;
        return download(body + '\n', name);
      }
      return suffix === '/pattern' ? patternResponse(ogIdx) : orthogroupDetails(ogIdx);
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
      let body: {by?: unknown; genomes?: unknown; q?: unknown; og?: unknown};
      try {
        body = JSON.parse(raw) as {by?: unknown; genomes?: unknown; q?: unknown; og?: unknown};
      } catch {
        return bad(400, 'invalid JSON body');
      }
      if (!body || typeof body !== 'object') return bad(400, 'invalid JSON body');
      const by = body.by;
      if (by === 'search') {
        const term = body.q;
        if (typeof term !== 'string' || !term.trim()) return bad(400, "'q' must be a non-empty string");
        const hitSet = new Set(computeHits(term, orders.tree).map((h) => h.og_idx));
        const perm = [...orders.tree.filter((og) => hitSet.has(og)), ...orders.tree.filter((og) => !hitSet.has(og))];
        const id = `search-${++dynamicId}`;
        orders[id] = perm;
        return json({order: id});
      }
      if (by === 'pattern') {
        const idm = typeof body.og === 'string' ? /^OG(\d+)$/.exec(body.og) : null;
        if (!idm || Number(idm[1]) >= N_OG) return bad(400, "'og' must be a known orthogroup id");
        const set = new Set(patternGroups.get(ogSignature(Number(idm[1])))!);
        const perm = [...orders.tree.filter((og) => set.has(og)), ...orders.tree.filter((og) => !set.has(og))];
        const id = `pattern-${++dynamicId}`;
        orders[id] = perm;
        return json({order: id});
      }
      if (by !== 'contrast') {
        return bad(400, "unsupported order; expected 'contrast', 'search' or 'pattern'");
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
