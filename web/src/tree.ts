// Genome dendrogram (rectangular cladogram) on a Canvas2D aligned with matrix rows.

import {getTree, type TreeArrays} from './api';
import {prepareCanvas, type Transform} from './transform';

const PAD = 8;
const HIT_RADIUS = 7; // px

export class Tree {
  private arrays: TreeArrays | null = null;
  private children: number[][] = [];
  private lo: number[] = [];
  private hi: number[] = [];
  private depth: number[] = [];
  private maxDepth = 1;
  selected: number | null = null;

  // cached mapping from the last draw, for hit testing
  private t = {xOf: (_: number) => 0, nodeY: (_: number) => 0};

  constructor(private canvas: HTMLCanvasElement) {}

  async load(): Promise<void> {
    this.arrays = await getTree();
    const a = this.arrays;
    const n = a.parent.length;
    this.children = Array.from({length: n}, () => [] as number[]);
    this.depth = new Array<number>(n).fill(0);
    this.lo = new Array<number>(n).fill(0);
    this.hi = new Array<number>(n).fill(0);
    for (let i = 0; i < n; i++) {
      const p = a.parent[i];
      // Guard against cycles / non-preorder parents: only a strictly earlier
      // node may be a parent.
      if (p >= 0 && p < i) {
        this.children[p].push(i);
        this.depth[i] = this.depth[p] + (a.length[i] || 0);
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
  }

  /** Genome-row range covered by the selected node, inclusive. */
  get selectedRange(): [number, number] | null {
    if (this.selected === null) return null;
    return [this.lo[this.selected], this.hi[this.selected]];
  }

  draw(tr: Transform, pw: number, ph: number): void {
    const ctx = prepareCanvas(this.canvas, pw, ph);
    if (!ctx) return;
    if (!this.arrays) return;

    const xOf = (d: number) => PAD + (d / this.maxDepth) * (pw - 2 * PAD);
    const nodeY = (i: number) => tr.yOf((this.lo[i] + this.hi[i]) / 2 + 0.5);
    this.t = {xOf, nodeY};
    const a = this.arrays;

    // selected clade highlight
    const range = this.selectedRange;
    if (range) {
      ctx.fillStyle = 'rgba(37,99,235,0.12)';
      const y0 = tr.yOf(range[0]);
      const y1 = tr.yOf(range[1] + 1);
      ctx.fillRect(0, y0, pw, y1 - y0);
    }

    ctx.lineWidth = 1;
    ctx.strokeStyle = this.selected !== null ? '#334155' : '#475569';
    ctx.beginPath();
    for (let i = 0; i < a.parent.length; i++) {
      const kids = this.children[i];
      if (kids.length) {
        const x = xOf(this.depth[i]);
        ctx.moveTo(x, nodeY(kids[0]));
        ctx.lineTo(x, nodeY(kids[kids.length - 1]));
      }
      const p = a.parent[i];
      if (p >= 0 && p < i) {
        const y = nodeY(i);
        ctx.moveTo(xOf(this.depth[p]), y);
        ctx.lineTo(xOf(this.depth[i]), y);
      }
    }
    ctx.stroke();

    // selected path in accent colour (cycle-guarded)
    if (this.selected !== null) {
      ctx.strokeStyle = '#2563eb';
      ctx.lineWidth = 1.6;
      ctx.beginPath();
      const seen = new Set<number>();
      let i: number = this.selected;
      while (i >= 0 && i < a.parent.length && !seen.has(i)) {
        seen.add(i);
        const p = a.parent[i];
        if (p < 0 || p >= i) break;
        const y = nodeY(i);
        ctx.moveTo(xOf(this.depth[p]), y);
        ctx.lineTo(xOf(this.depth[i]), y);
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
      // The horizontal branch of node i runs from its parent's depth to its own.
      const x0 = p >= 0 && p < i ? xOf(this.depth[p]) : xOf(this.depth[i]) - 6;
      const x1 = xOf(this.depth[i]);
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
