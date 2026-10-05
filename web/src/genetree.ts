// Gene-tree layout and leaf colours (pure; the alignment view renders it).

import type {GeneTreeArrays, Genome, GenomeMetadataInfo} from './api';
import {categoryColor} from './annotations';
import {numericColor} from './metadata';

// --- pure layout ------------------------------------------------------------

export interface GeneTreeLayout {
  n: number;
  children: number[][];
  depth: number[];
  /** number of edges from the root (for cladogram layouts) */
  edges: number[];
  lo: number[];
  hi: number[];
  leafIndex: number[];
  nLeaves: number;
  maxDepth: number;
  maxEdges: number;
}

/** Pre-order node arrays -> children/depth/edge-depth/leaf ranges (see tree.ts). */
export function layoutGeneTree(a: GeneTreeArrays): GeneTreeLayout {
  const n = a.parent.length;
  const children: number[][] = Array.from({length: n}, () => []);
  const depth = new Array<number>(n).fill(0);
  const edges = new Array<number>(n).fill(0);
  for (let i = 0; i < n; i++) {
    const p = a.parent[i];
    if (p >= 0 && p < i) {
      children[p].push(i);
      depth[i] = depth[p] + (a.length[i] || 0);
      edges[i] = edges[p] + 1;
    }
  }
  const lo = new Array<number>(n).fill(0);
  const hi = new Array<number>(n).fill(0);
  const leafIndex = new Array<number>(n).fill(-1);
  let leaf = 0;
  for (let i = 0; i < n; i++) {
    if (children[i].length === 0) {
      lo[i] = hi[i] = leaf;
      leafIndex[i] = leaf++;
    }
  }
  for (let i = n - 1; i >= 0; i--) {
    const kids = children[i];
    if (kids.length) {
      lo[i] = lo[kids[0]];
      hi[i] = hi[kids[kids.length - 1]];
    }
  }
  return {
    n,
    children,
    depth,
    edges,
    lo,
    hi,
    leafIndex,
    nLeaves: leaf,
    maxDepth: Math.max(1e-9, ...depth),
    maxEdges: Math.max(1, ...edges),
  };
}

// --- colours ----------------------------------------------------------------

function hashString(s: string): number {
  let h = 2166136261;
  for (let i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  return h >>> 0;
}

export const MISSING_COLOR = '#98a2b3';

/** Leaf colour from a genome metadata field, else a stable genome identity colour. */
export function leafColor(genome: Genome | undefined, spec: GenomeMetadataInfo | undefined): string {
  if (!genome) return MISSING_COLOR;
  if (!spec) return categoryColor(hashString(genome.name) % 12);
  const value = genome.metadata?.[spec.name];
  if (value === null || value === undefined || value === '') return MISSING_COLOR;
  if (spec.kind === 'categorical') {
    const idx = (spec.categories ?? []).indexOf(String(value));
    return idx < 0 ? MISSING_COLOR : categoryColor(idx);
  }
  const min = spec.min ?? 0;
  const max = spec.max ?? 1;
  return numericColor((Number(value) - min) / (max - min || 1));
}

// --- pure view helpers ------------------------------------------------------

/** X of a node in a phylogram (branch length) or cladogram (edge count, tips aligned). */
export function treeNodeX(
  mode: 'phylogram' | 'cladogram',
  depth: number,
  maxDepth: number,
  edges: number,
  maxEdges: number,
  isLeaf: boolean,
  width: number,
  pad: number,
): number {
  const span = width - 2 * pad;
  if (mode === 'cladogram') return pad + ((isLeaf ? maxEdges : edges) / maxEdges) * span;
  return pad + (depth / maxDepth) * span;
}
