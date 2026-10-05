// Genome dendrogram (rectangular cladogram) on a Canvas2D aligned with matrix rows.

import {getTree, type TreeArrays} from './api';
import type {Theme} from './theme';
import {paletteFor} from './theme';
import {prepareCanvas, type Transform} from './transform';

const PAD = 8;
const HIT_RADIUS = 7; // px

export class Tree {
  private arrays: TreeArrays | null = null;
  private children: number[][] = [];
  private lo: number[] = [];
  private hi: number[] = [];
  private depth: number[] = [];
  private edges: number[] = [];
  private maxDepth = 1;
  private maxEdges = 1;
  /** phylogram uses branch lengths; cladogram uses edge counts (always visible). */
  mode: 'phylogram' | 'cladogram' = 'phylogram';
  selected: number | null = null;

  // cached mapping from the last draw, for hit testing
  private t = {xOf: (_: number) => 0, nodeY: (_: number) => 0};

  constructor(private canvas: HTMLCanvasElement) {}

  setMode(mode: 'phylogram' | 'cladogram'): void {
    this.mode = mode;
  }

  async load(): Promise<void> {
    this.arrays = await getTree();
    const a = this.arrays;
    const n = a.parent.length;
    this.children = Array.from({length: n}, () => [] as number[]);
    this.depth = new Array<number>(n).fill(0);
    this.edges = new Array<number>(n).fill(0);
    this.lo = new Array<number>(n).fill(0);
    this.hi = new Array<number>(n).fill(0);
    for (let i = 0; i < n; i++) {
      const p = a.parent[i];
      // Guard against cycles / non-preorder parents: only a strictly earlier
      // node may be a parent.
      if (p >= 0 && p < i) {
        this.children[p].push(i);
        this.depth[i] = this.depth[p] + (a.length[i] || 0);
        this.edges[i] = this.edges[p] + 1;
      }
    }
    for (let i = n - 1; i >= 0; i--) {
      const kids = this.children[i];
      if (kids.length === 0) {
        this.lo[i] = this.hi[i] = a.row[i];
      } else {
        this.lo[i] = this.lo[kids[0]];
        this.hi[i] = this.hi[kids[kids.length - 1]];
      }
    }
    this.maxDepth = Math.max(1e-9, ...this.depth);
    this.maxEdges = Math.max(1, ...this.edges);
  }

  /**
   * Edge depth is used when the requested mode is cladogram, or when branch
   * lengths are degenerate (all zero / non-finite) so the topology stays visible.
   */
  private useEdges(): boolean {
    return this.mode === 'cladogram' || !(this.maxDepth > 1e-6 && Number.isFinite(this.maxDepth));
  }

  /** X pixel of node i's branch end for the current mode. */
  private xOfNode(i: number, pw: number): number {
    const span = pw - 2 * PAD;
    if (this.useEdges()) {
      const e = this.children[i].length === 0 ? this.maxEdges : this.edges[i];
      return PAD + (e / this.maxEdges) * span;
    }
    return PAD + ((this.depth[i] || 0) / this.maxDepth) * span;
  }

  /** Genome-row range covered by the selected node, inclusive. */
  get selectedRange(): [number, number] | null {
    if (this.selected === null) return null;
    return [this.lo[this.selected], this.hi[this.selected]];
  }

  draw(tr: Transform, pw: number, ph: number, theme: Theme = 'light'): void {
    const ctx = prepareCanvas(this.canvas, pw, ph);
    if (!ctx) return;
    if (!this.arrays) return;
    const pal = paletteFor(theme);

    const xOf = (i: number) => this.xOfNode(i, pw);
    const nodeY = (i: number) => tr.yOf((this.lo[i] + this.hi[i]) / 2 + 0.5);
    this.t = {xOf, nodeY};
    const a = this.arrays;

    // selected clade highlight
    const range = this.selectedRange;
    if (range) {
      ctx.fillStyle = pal.cladeTint;
      const y0 = tr.yOf(range[0]);
      const y1 = tr.yOf(range[1] + 1);
      ctx.fillRect(0, y0, pw, y1 - y0);
    }

    ctx.lineWidth = 1;
    ctx.strokeStyle = this.selected !== null ? pal.treeLineSelected : pal.treeLine;
    ctx.beginPath();
    for (let i = 0; i < a.parent.length; i++) {
      const kids = this.children[i];
      if (kids.length) {
        const x = xOf(i);
        ctx.moveTo(x, nodeY(kids[0]));
        ctx.lineTo(x, nodeY(kids[kids.length - 1]));
      }
      const p = a.parent[i];
      if (p >= 0 && p < i) {
        const y = nodeY(i);
        ctx.moveTo(xOf(p), y);
        ctx.lineTo(xOf(i), y);
      }
    }
    ctx.stroke();

    // selected path in accent colour (cycle-guarded)
    if (this.selected !== null) {
      ctx.strokeStyle = pal.treePath;
      ctx.lineWidth = 1.6;
      ctx.beginPath();
      const seen = new Set<number>();
      let i: number = this.selected;
      while (i >= 0 && i < a.parent.length && !seen.has(i)) {
        seen.add(i);
        const p = a.parent[i];
        if (p < 0 || p >= i) break;
        const y = nodeY(i);
        ctx.moveTo(xOf(p), y);
        ctx.lineTo(xOf(i), y);
        i = p;
      }
      ctx.stroke();
    }
  }

  /** Pick the tree node whose horizontal branch is nearest to a pixel position. */
  hitTest(px: number, py: number): number | null {
    if (!this.arrays) return null;
    const {xOf, nodeY} = this.t;
    const a = this.arrays;
    let best: number | null = null;
    let bestD = HIT_RADIUS;
    for (let i = 0; i < a.parent.length; i++) {
      const y = nodeY(i);
      const p = a.parent[i];
      // The horizontal branch of node i runs from its parent's branch end to its own.
      const x0 = p >= 0 && p < i ? xOf(p) : xOf(i) - 6;
      const x1 = xOf(i);
      const d = segmentDistance(px, py, x0, y, x1, y);
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    }
    return best;
  }
}

function segmentDistance(px: number, py: number, x0: number, y0: number, x1: number, y1: number): number {
  const dx = x1 - x0;
  const dy = y1 - y0;
  const len2 = dx * dx + dy * dy;
  const t = len2 === 0 ? 0 : Math.min(1, Math.max(0, ((px - x0) * dx + (py - y0) * dy) / len2));
  const cx = x0 + t * dx;
  const cy = y0 + t * dy;
  return Math.hypot(px - cx, py - cy);
}
