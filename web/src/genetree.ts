// Gene-tree tab: pure pre-order layout plus a rectangular Canvas2D renderer.

import type {GeneTreeArrays, Genome, GenomeMetadataInfo} from './api';
import {categoryColor} from './annotations';
import {numericColor} from './metadata';
import {prepareCanvas} from './transform';

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

const ROW_H = 15;
export const MIN_TREE_ROW = 4;
export const HIDE_LABEL_ROW = 9;

/** Fit all leaves into viewH, but never below MIN_TREE_ROW px per leaf. */
export function fitTreeScale(nLeaves: number, viewH: number, rowH = ROW_H): number {
  if (nLeaves <= 0) return 1;
  return Math.min(1, Math.max(MIN_TREE_ROW / rowH, viewH / (nLeaves * rowH)));
}

/** Hide leaf labels once rows are too short to read. */
export const hideLeafLabels = (rowH: number) => rowH < HIDE_LABEL_ROW;

export function clampTreeScroll(scrollY: number, contentH: number, viewH: number): number {
  return Math.max(0, Math.min(scrollY, Math.max(0, contentH - viewH)));
}

/** X of a node in a phylogram (branch length) or cladogram (edge count, tips aligned). */
export function treeNodeX(mode: 'phylogram' | 'cladogram', depth: number, maxDepth: number, edges: number, maxEdges: number, isLeaf: boolean, width: number, pad: number): number {
  const span = width - 2 * pad;
  if (mode === 'cladogram') return pad + ((isLeaf ? maxEdges : edges) / maxEdges) * span;
  return pad + (depth / maxDepth) * span;
}

// --- viewer -----------------------------------------------------------------

const PAD = 10;
const LABEL_W = 240;

export interface GeneTreeHover {
  node: number;
  gene: string | null;
  genome: string | null;
  duplication: boolean;
  support: number | null;
}

export class GeneTreeView {
  private data: GeneTreeArrays | null = null;
  private layout: GeneTreeLayout | null = null;
  private genomes: Genome[] = [];
  private spec: GenomeMetadataInfo | undefined = undefined;
  private mode: 'phylogram' | 'cladogram' = 'phylogram';
  private scale = 1;
  private userZoomed = false;
  private scrollY = 0;
  private observer: ResizeObserver;
  private t = {xOf: (_: number) => 0, yOf: (_: number) => 0};

  constructor(
    private canvas: HTMLCanvasElement,
    private onHover: (h: GeneTreeHover | null, clientX: number, clientY: number) => void,
    private onLeafHover: (genomeRow: number | null) => void,
  ) {
    canvas.addEventListener('wheel', (e) => this.wheel(e), {passive: false});
    canvas.addEventListener('mousemove', (e) => this.mousemove(e));
    canvas.addEventListener('mouseleave', () => this.clearHover());
    this.observer = new ResizeObserver(() => this.onResize());
    this.observer.observe(canvas.parentElement ?? canvas);
  }

  private clearHover(): void {
    this.onHover(null, 0, 0);
    this.onLeafHover(null);
  }

  setData(data: GeneTreeArrays | null, genomes: Genome[], spec: GenomeMetadataInfo | undefined): void {
    this.data = data;
    this.layout = data ? layoutGeneTree(data) : null;
    this.genomes = genomes;
    this.spec = spec;
    this.userZoomed = false;
    this.scrollY = 0;
    this.scale = this.layout ? fitTreeScale(this.layout.nLeaves, this.viewH()) : 1;
    this.draw();
  }

  setSpec(spec: GenomeMetadataInfo | undefined): void {
    this.spec = spec;
    this.draw();
  }

  setMode(mode: 'phylogram' | 'cladogram'): void {
    this.mode = mode;
    this.draw();
  }

  private viewH(): number {
    return Math.max(1, this.canvas.clientHeight - 16);
  }

  private contentH(): number {
    return this.layout ? 16 + this.layout.nLeaves * ROW_H * this.scale : 0;
  }

  private onResize(): void {
    if (this.layout && !this.userZoomed) this.scale = fitTreeScale(this.layout.nLeaves, this.viewH());
    this.scrollY = clampTreeScroll(this.scrollY, this.contentH(), this.viewH());
    this.draw();
  }

  zoomBy(factor: number): void {
    this.scale = Math.min(6, Math.max(0.2, this.scale * factor));
    if (factor !== 1) this.userZoomed = true;
    this.scrollY = clampTreeScroll(this.scrollY, this.contentH(), this.viewH());
    this.draw();
  }

  private wheel(e: WheelEvent): void {
    if (!this.data) return;
    e.preventDefault();
    this.clearHover();
    if (e.ctrlKey || e.metaKey) this.zoomBy(2 ** (-e.deltaY * 0.0015));
    else {
      this.scrollY = clampTreeScroll(this.scrollY + e.deltaY, this.contentH(), this.viewH());
      this.draw();
    }
  }

  private colorFor(genomeRow: number): string {
    if (genomeRow < 0) return MISSING_COLOR;
    return leafColor(this.genomes[genomeRow], this.spec);
  }

  draw(): void {
    const cssW = this.canvas.clientWidth;
    const cssH = this.canvas.clientHeight;
    const ctx = prepareCanvas(this.canvas, cssW, cssH);
    if (!ctx) return;
    ctx.fillStyle = '#fff';
    ctx.fillRect(0, 0, cssW, cssH);
    if (!this.data || !this.layout) {
      ctx.fillStyle = '#98a2b3';
      ctx.font = '12px system-ui, sans-serif';
      ctx.textBaseline = 'top';
      ctx.fillText('no gene tree', 12, 12);
      return;
    }
    const a = this.data;
    const L = this.layout;
    const treeW = Math.max(1, cssW - LABEL_W);
    const rowH = ROW_H * this.scale;
    const xOf = (i: number) => treeNodeX(this.mode, L.depth[i], L.maxDepth, L.edges[i], L.maxEdges, L.leafIndex[i] >= 0, treeW, PAD);
    const yOf = (i: number) => 8 + (L.leafIndex[i] >= 0 ? L.leafIndex[i] : (L.lo[i] + L.hi[i]) / 2) * rowH - this.scrollY;
    this.t = {xOf, yOf};

    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, treeW, cssH);
    ctx.clip();
    ctx.lineWidth = 1;
    ctx.strokeStyle = '#475569';
    ctx.beginPath();
    for (let i = 0; i < L.n; i++) {
      const kids = L.children[i];
      if (kids.length) {
        const x = xOf(i);
        ctx.moveTo(x, yOf(kids[0]));
        ctx.lineTo(x, yOf(kids[kids.length - 1]));
      }
      const p = a.parent[i];
      if (p >= 0 && p < i) {
        const y = yOf(i);
        ctx.moveTo(xOf(p), y);
        ctx.lineTo(xOf(i), y);
      }
    }
    ctx.stroke();
    // dotted leader lines from each tip to the label column
    ctx.save();
    ctx.setLineDash([2, 3]);
    ctx.strokeStyle = '#cbd5e1';
    ctx.beginPath();
    for (let i = 0; i < L.n; i++) {
      if (L.leafIndex[i] < 0) continue;
      const y = yOf(i);
      if (y < -rowH || y > cssH + rowH) continue;
      ctx.moveTo(xOf(i), y);
      ctx.lineTo(treeW, y);
    }
    ctx.stroke();
    ctx.restore();
    // duplication markers: hollow black diamond
    for (let i = 0; i < L.n; i++) {
      if (!a.duplication[i]) continue;
      const x = xOf(i);
      const y = yOf(i);
      if (y < -6 || y > cssH + 6) continue;
      ctx.beginPath();
      ctx.moveTo(x, y - 4);
      ctx.lineTo(x + 4, y);
      ctx.lineTo(x, y + 4);
      ctx.lineTo(x - 4, y);
      ctx.closePath();
      ctx.fillStyle = '#fff';
      ctx.fill();
      ctx.strokeStyle = '#111827';
      ctx.lineWidth = 1.4;
      ctx.stroke();
    }
    ctx.restore();

    // leaf labels, coloured by the metadata field
    const hideLabels = hideLeafLabels(rowH);
    ctx.save();
    ctx.beginPath();
    ctx.rect(treeW, 0, LABEL_W, cssH);
    ctx.clip();
    ctx.font = '11px system-ui, sans-serif';
    ctx.textBaseline = 'middle';
    for (let i = 0; i < L.n; i++) {
      if (L.leafIndex[i] < 0) continue;
      const y = yOf(i);
      if (y < -rowH || y > cssH + rowH) continue;
      const gr = a.genome_row[i];
      ctx.fillStyle = this.colorFor(gr);
      ctx.fillRect(treeW + 2, y - Math.min(6, rowH / 2 - 1), 8, Math.min(12, Math.max(3, rowH - 2)));
      if (hideLabels) continue;
      ctx.fillStyle = '#1f2937';
      const gene = a.gene[i] ?? a.label[i];
      ctx.fillText(gene ?? '', treeW + 14, y);
    }
    ctx.restore();
  }

  private mousemove(e: MouseEvent): void {
    if (!this.data || !this.layout) return;
    const rect = this.canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    const {xOf, yOf} = this.t;
    let best = -1;
    let bestD = 12;
    for (let i = 0; i < this.layout.n; i++) {
      const d = Math.hypot(x - xOf(i), y - yOf(i));
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    }
    if (best < 0) {
      this.clearHover();
      return;
    }
    const gr = this.data.genome_row[best];
    this.onHover(
      {node: best, gene: this.data.gene[best], genome: gr >= 0 ? this.genomes[gr]?.name ?? null : null, duplication: !!this.data.duplication[best], support: this.data.support[best]},
      e.clientX,
      e.clientY,
    );
    this.onLeafHover(this.layout.leafIndex[best] >= 0 && gr >= 0 ? gr : null);
  }
}
