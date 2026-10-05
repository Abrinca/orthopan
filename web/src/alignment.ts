// Single-orthogroup alignment view: pure row-order / virtual-band / neighbour
// math plus a virtualised Canvas2D renderer (tree | colour | name | MSA).
//
// The view occupies the whole window. Vertical scrolling is the overlay's own
// scroll; only the MSA scrolls horizontally. A sticky canvas draws just the
// visible vertical band, so 2,000-row alignments never allocate a huge canvas.

import type {GeneTreeArrays, Genome, GenomeMetadataInfo, MsaResult} from './api';
import {layoutGeneTree, leafColor, MISSING_COLOR, treeNodeX} from './genetree';
import {clampPxPerCol, clampScrollCol, fitPxPerCol, residueColor, scrollFromBar, showLetters, visibleCols} from './msa';
import {prepareCanvas} from './transform';

// --- constants --------------------------------------------------------------

export const ALN_ROW_H = 16;
export const ALN_TREE_W = 300;
export const ALN_COLOR_W = 22;
export const ALN_NAME_W = 210;
export const ALN_MSA_X = ALN_TREE_W + ALN_COLOR_W + ALN_NAME_W;

export const ALN_BAR_H = 38;
const CONS_Y = 2;
const CONS_H = 12;
const GAP_Y = 16;
const GAP_H = 12;
const SB_Y = 30;
const SB_H = 7;

// --- pure row order ---------------------------------------------------------

export interface AlnRow {
  /** display gene id */
  gene: string | null;
  /** raw gene id */
  raw: string | null;
  /** genome row index, -1 when unknown */
  genomeRow: number;
  /** genome name (fallback when genomeRow < 0) */
  genome: string;
  /** aligned sequence, null when no MSA row matched */
  seq: string | null;
  /** index into msa.rows, -1 when none */
  msaIndex: number;
  /** gene-tree leaf node index, -1 when no tree */
  leaf: number;
}

export interface AlnOrder {
  rows: AlnRow[];
  hasTree: boolean;
}

/**
 * Row order: tree leaf order when a tree exists (so tree leaves, colour cells,
 * names and MSA rows line up exactly), otherwise the MSA/genome-tree order.
 */
export function orderAlignmentRows(msa: MsaResult | null, tree: GeneTreeArrays | null): AlnOrder {
  const msaRows = msa?.rows ?? [];
  const byKey = new Map<string, number>();
  msaRows.forEach((r, i) => {
    if (r.gene) byKey.set(r.gene, i);
    if (r.raw) byKey.set(r.raw, i);
  });

  if (tree && tree.parent.length) {
    const L = layoutGeneTree(tree);
    if (L.nLeaves > 0) {
      const rows: AlnRow[] = [];
      for (let i = 0; i < L.n; i++) {
        if (L.leafIndex[i] < 0) continue;
        let mi = -1;
        for (const key of [tree.gene[i], tree.label[i]]) {
          if (key === null || key === undefined) continue;
          const v = byKey.get(key);
          if (v !== undefined) {
            mi = v;
            break;
          }
        }
        const mr = mi >= 0 ? msaRows[mi] : null;
        rows.push({
          gene: mr?.gene ?? tree.gene[i] ?? tree.label[i] ?? null,
          raw: mr?.raw ?? tree.label[i] ?? null,
          genomeRow: tree.genome_row[i] >= 0 ? tree.genome_row[i] : (mr?.row ?? -1),
          genome: mr?.genome ?? '',
          seq: mr?.seq ?? null,
          msaIndex: mi,
          leaf: i,
        });
      }
      return {rows, hasTree: true};
    }
  }

  const rows: AlnRow[] = msaRows.map((r, i) => ({
    gene: r.gene,
    raw: r.raw,
    genomeRow: r.row,
    genome: r.genome,
    seq: r.seq,
    msaIndex: i,
    leaf: -1,
  }));
  return {rows, hasTree: false};
}

// --- pure virtual band ------------------------------------------------------

export interface Band {
  startRow: number;
  endRow: number;
  /** y of `startRow` relative to the viewport top */
  offset: number;
}

/** Rows intersecting [scrollTop, scrollTop+viewH), plus `overscan` on each side. */
export function visibleBand(
  scrollTop: number,
  viewH: number,
  nRows: number,
  rowH = ALN_ROW_H,
  overscan = 2,
): Band {
  const top = Math.max(0, scrollTop);
  const start = Math.max(0, Math.floor(top / rowH) - overscan);
  const end = Math.min(nRows, Math.ceil((top + Math.max(0, viewH)) / rowH) + overscan);
  return {startRow: start, endRow: Math.max(start, end), offset: start * rowH - scrollTop};
}

// --- pure neighbour logic ---------------------------------------------------/** Position of `og` within a known id window, or null when outside it. */
export function locateOg(ids: string[], og0: number, og: string): number | null {
  const i = ids.indexOf(og);
  return i < 0 ? null : og0 + i;
}

/** Neighbouring position within [0, n), or null at either end. */
export function stepPosition(pos: number, delta: number, n: number): number | null {
  const next = pos + delta;
  return next >= 0 && next < n ? next : null;
}

/** 1-based position of `col` among the non-gap residues, or -1 for a gap. */
export function ungappedPosition(seq: string, col: number): number {
  const here = seq[col];
  if (here === undefined || here === '-') return -1;
  let n = 1;
  for (let i = 0; i < col; i++) if (seq[i] !== '-') n++;
  return n;
}

/** Greedy word wrap for short canvas messages. */
function wrapText(ctx: CanvasRenderingContext2D, text: string, maxW: number): string[] {
  const words = text.split(/\s+/).filter(Boolean);
  const lines: string[] = [];
  let cur = '';
  for (const w of words) {
    const t = cur ? `${cur} ${w}` : w;
    if (cur && ctx.measureText(t).width > maxW) {
      lines.push(cur);
      cur = w;
    } else {
      cur = t;
    }
  }
  if (cur) lines.push(cur);
  return lines.length ? lines : [text];
}

// --- viewer -----------------------------------------------------------------

export interface AlnContent {
  og: string;
  rows: AlnRow[];
  hasTree: boolean;
  treeMessage: string | null;
  msaLength: number;
  conservation: Float32Array | null;
  gap: Float32Array | null;
  msaMessage: string | null;
  genomes: Genome[];
  spec: GenomeMetadataInfo | undefined;
  idMode: 'display' | 'raw';
  /** genome names to tint (the drawer's selected genome/range) */
  selectedGenomes: string[];
}

export interface AlnHover {
  row: number;
  col: number;
  gene: string;
  genome: string;
  residue: string;
  ungapped: number;
  conservation: number;
  gap: number;
}

export interface AlnBarHover {
  col: number;
  conservation: number;
  gap: number;
}

type Mode = 'phylogram' | 'cladogram';

export class AlignmentView {
  private content: AlnContent | null = null;
  private pxPerCol = 8;
  private scrollCol = 0;
  private userZoomed = false;
  private hoverRow: number | null = null;
  private dragging = false;
  private mode: Mode = 'phylogram';
  private observer: ResizeObserver;
  private raf = 0;

  constructor(
    private scroll: HTMLElement,
    private rowsEl: HTMLElement,
    private canvas: HTMLCanvasElement,
    private bar: HTMLCanvasElement,
    private onHover: (h: AlnHover | null, clientX: number, clientY: number) => void,
    private onBarHover: (h: AlnBarHover | null, clientX: number, clientY: number) => void,
  ) {
    this.scroll.addEventListener('scroll', () => this.scheduleDraw(), {passive: true});
    canvas.addEventListener('mousemove', (e) => this.mousemove(e));
    canvas.addEventListener('mouseleave', () => {
      this.hoverRow = null;
      this.onHover(null, 0, 0);
      this.scheduleDraw();
    });
    canvas.addEventListener('wheel', (e) => this.wheel(e), {passive: false});
    bar.addEventListener('wheel', (e) => this.wheel(e), {passive: false});
    bar.addEventListener('mousemove', (e) => this.barMousemove(e));
    bar.addEventListener('mouseleave', () => this.onBarHover(null, 0, 0));
    bar.addEventListener('pointerdown', (e) => this.barPointerdown(e));
    window.addEventListener('pointermove', (e) => this.pointermove(e));
    window.addEventListener('pointerup', () => (this.dragging = false));
    this.observer = new ResizeObserver(() => this.syncSize());
    this.observer.observe(scroll);
  }

  setData(c: AlnContent): void {
    this.content = c;
    this.hoverRow = null;
    this.scrollCol = 0;
    this.userZoomed = false;
    this.scroll.scrollTop = 0;
    this.rowsEl.style.height = `${c.rows.length * ALN_ROW_H}px`;
    this.syncSize();
  }

  /** The OG currently rendered, or null. */
  get contentOg(): string | null {
    return this.content?.og ?? null;
  }

  /** Drop the current alignment (on OG change) and clear hover immediately. */
  clear(): void {
    this.content = null;
    this.hoverRow = null;
    this.rowsEl.style.height = '0px';
    this.onHover(null, 0, 0);
    this.scheduleDraw();
  }

  setSpec(spec: GenomeMetadataInfo | undefined): void {
    if (!this.content) return;
    this.content.spec = spec;
    this.scheduleDraw();
  }

  setIdMode(idMode: 'display' | 'raw'): void {
    if (!this.content) return;
    this.content.idMode = idMode;
    this.scheduleDraw();
  }

  setMode(mode: Mode): void {
    this.mode = mode;
    this.scheduleDraw();
  }

  get modeValue(): Mode {
    return this.mode;
  }

  /** Scroll so that `row` sits near the top of the viewport. */
  scrollToRow(row: number): void {
    this.scroll.scrollTop = Math.max(0, row * ALN_ROW_H);
  }

  get pxPerColValue(): number {
    return this.pxPerCol;
  }

  zoomBy(factor: number): void {
    if (!this.content) return;
    const width = this.msaW();
    const anchor = this.scrollCol + width / this.pxPerCol / 2;
    this.pxPerCol = clampPxPerCol(this.pxPerCol * factor);
    this.scrollCol = clampScrollCol(anchor - width / this.pxPerCol / 2, width, this.pxPerCol, this.content.msaLength);
    if (factor !== 1) this.userZoomed = true;
    this.scheduleDraw();
  }

  /** Reset the overlay scroll to the top. */
  resetScroll(): void {
    this.scroll.scrollTop = 0;
  }

  private msaW(): number {
    return Math.max(1, this.canvas.clientWidth - ALN_MSA_X);
  }

  private syncSize(): void {
    const h = Math.max(1, this.scroll.clientHeight);
    const w = Math.max(1, this.scroll.clientWidth);
    this.canvas.style.height = `${h}px`;
    // Match the conservation bar exactly to the scrollable body width so the
    // MSA columns line up, even when a vertical scrollbar shrinks the body.
    this.bar.style.width = `${w}px`;
    if (!this.content) return;
    if (!this.userZoomed) this.pxPerCol = fitPxPerCol(this.msaW(), this.content.msaLength);
    this.scrollCol = clampScrollCol(this.scrollCol, this.msaW(), this.pxPerCol, this.content.msaLength);
    this.scheduleDraw();
  }

  private scheduleDraw(): void {
    if (this.raf) return;
    this.raf = requestAnimationFrame(() => {
      this.raf = 0;
      this.draw();
    });
  }

  private colAt(x: number): number {
    const c = this.content!;
    return Math.min(c.msaLength - 1, Math.max(0, Math.floor(this.scrollCol + (x - ALN_MSA_X) / this.pxPerCol)));
  }

  private wheel(e: WheelEvent): void {
    const c = this.content;
    if (!c) return;
    const rect = (e.currentTarget as HTMLCanvasElement).getBoundingClientRect();
    const x = e.clientX - rect.left;
    const overMsa = x >= ALN_MSA_X;
    if (e.ctrlKey || e.metaKey) {
      if (overMsa) {
        e.preventDefault();
        const anchorPx = x - ALN_MSA_X;
        const anchorCol = this.scrollCol + anchorPx / this.pxPerCol;
        const next = clampPxPerCol(this.pxPerCol * 2 ** (-e.deltaY * 0.0015));
        if (next !== this.pxPerCol) {
          this.pxPerCol = next;
          this.scrollCol = clampScrollCol(anchorCol - anchorPx / next, this.msaW(), next, c.msaLength);
          this.userZoomed = true;
          this.scheduleDraw();
        }
      }
      return;
    }
    const horizontal = e.shiftKey || Math.abs(e.deltaX) > Math.abs(e.deltaY);
    if (overMsa && horizontal) {
      e.preventDefault();
      const delta = (e.shiftKey ? e.deltaY : e.deltaX) / this.pxPerCol;
      this.scrollCol = clampScrollCol(this.scrollCol + delta, this.msaW(), this.pxPerCol, c.msaLength);
      this.scheduleDraw();
    }
    // vertical wheel is left to the page/overlay scroll
  }

  private mousemove(e: MouseEvent): void {
    const c = this.content;
    if (!c || !c.rows.length) return;
    const rect = this.canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    const row = Math.floor((this.scroll.scrollTop + y) / ALN_ROW_H);
    if (row < 0 || row >= c.rows.length) {
      this.onHover(null, e.clientX, e.clientY);
      return;
    }
    if (this.hoverRow !== row) {
      this.hoverRow = row;
      this.scheduleDraw();
    }
    const r = c.rows[row];
    const genome = r.genomeRow >= 0 ? c.genomes[r.genomeRow]?.name ?? r.genome : r.genome;
    if (x < ALN_MSA_X || !c.conservation) {
      this.onHover({row, col: -1, gene: r.gene ?? '', genome, residue: '', ungapped: -1, conservation: 0, gap: 0}, e.clientX, e.clientY);
      return;
    }
    const col = this.colAt(x);
    const residue = r.seq?.[col] ?? '-';
    const ungapped = r.seq ? ungappedPosition(r.seq, col) : -1;
    this.onHover({row, col, gene: r.gene ?? '', genome, residue, ungapped, conservation: c.conservation[col], gap: c.gap![col]}, e.clientX, e.clientY);
  }

  private barMousemove(e: MouseEvent): void {
    const c = this.content;
    if (!c || !c.conservation) return;
    const rect = this.bar.getBoundingClientRect();
    const x = e.clientX - rect.left;
    if (x < ALN_MSA_X) {
      this.onBarHover(null, 0, 0);
      return;
    }
    const col = this.colAt(x);
    this.onBarHover({col, conservation: c.conservation[col], gap: c.gap![col]}, e.clientX, e.clientY);
  }

  private barPointerdown(e: PointerEvent): void {
    const c = this.content;
    if (!c) return;
    const rect = this.bar.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    if (x < ALN_MSA_X || y < SB_Y || y > SB_Y + SB_H) return;
    this.dragging = true;
    this.scrollCol = scrollFromBar((x - ALN_MSA_X) / this.msaW(), this.msaW(), this.pxPerCol, c.msaLength);
    this.scheduleDraw();
  }

  private pointermove(e: PointerEvent): void {
    const c = this.content;
    if (!this.dragging || !c) return;
    const rect = this.bar.getBoundingClientRect();
    this.scrollCol = scrollFromBar((e.clientX - rect.left - ALN_MSA_X) / this.msaW(), this.msaW(), this.pxPerCol, c.msaLength);
    this.scheduleDraw();
  }

  draw(): void {
    this.drawBody();
    this.drawBar();
  }

  private drawBody(): void {
    const cssW = this.canvas.clientWidth;
    const cssH = this.canvas.clientHeight;
    const ctx = prepareCanvas(this.canvas, cssW, cssH);
    if (!ctx) return;
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, cssW, cssH);
    const c = this.content;
    if (!c) return;
    const S = this.scroll.scrollTop;
    const n = c.rows.length;
    const msaW = Math.max(1, cssW - ALN_MSA_X);
    const nameX = ALN_TREE_W + ALN_COLOR_W;

    if (!c.hasTree && c.treeMessage) {
      ctx.save();
      ctx.beginPath();
      ctx.rect(0, 0, ALN_TREE_W, cssH);
      ctx.clip();
      ctx.fillStyle = '#98a2b3';
      ctx.font = '11px system-ui, sans-serif';
      ctx.textBaseline = 'top';
      let ty = 6;
      for (const line of wrapText(ctx, c.treeMessage, ALN_TREE_W - 12)) {
        ctx.fillText(line, 6, ty);
        ty += 13;
      }
      ctx.restore();
    }
    if (c.msaMessage && !c.conservation) {
      ctx.save();
      ctx.beginPath();
      ctx.rect(ALN_MSA_X, 0, msaW, cssH);
      ctx.clip();
      ctx.fillStyle = '#98a2b3';
      ctx.font = '11px system-ui, sans-serif';
      ctx.textBaseline = 'top';
      let my = 6;
      for (const line of wrapText(ctx, c.msaMessage, msaW - 12)) {
        ctx.fillText(line, ALN_MSA_X + 6, my);
        my += 13;
      }
      ctx.restore();
    }

    const band = visibleBand(S, cssH, n, ALN_ROW_H, 2);
    const letters = showLetters(this.pxPerCol);
    const selected = new Set(c.selectedGenomes);
    // pass 1: full-width row chrome (selection tint, highlight, colour, name, separator)
    ctx.font = '11px system-ui, sans-serif';
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'left';
    for (let i = band.startRow; i < band.endRow; i++) {
      const y = i * ALN_ROW_H - S;
      const r = c.rows[i];
      const rowGenome = r.genomeRow >= 0 ? c.genomes[r.genomeRow]?.name ?? r.genome : r.genome;
      if (rowGenome && selected.has(rowGenome)) {
        ctx.fillStyle = 'rgba(245,158,11,0.14)';
        ctx.fillRect(0, y, cssW, ALN_ROW_H);
      }
      if (i === this.hoverRow) {
        ctx.fillStyle = 'rgba(37,99,235,0.10)';
        ctx.fillRect(0, y, cssW, ALN_ROW_H);
      }
      const gcol = r.genomeRow >= 0 ? leafColor(c.genomes[r.genomeRow], c.spec) : MISSING_COLOR;
      ctx.fillStyle = gcol;
      ctx.fillRect(ALN_TREE_W, y, ALN_COLOR_W - 2, ALN_ROW_H - 1);
      ctx.fillStyle = '#1f2937';
      const name = c.idMode === 'raw' ? r.raw : r.gene;
      if (name) ctx.fillText(name, nameX + 4, y + ALN_ROW_H / 2, ALN_NAME_W - 8);
      ctx.fillStyle = '#f3f4f6';
      ctx.fillRect(0, y + ALN_ROW_H - 1, cssW, 1);
    }
    // pass 2: MSA residues with a single clip (cheap for tall alignments)
    if (c.conservation) {
      const [c0, c1] = visibleCols(this.scrollCol, msaW, this.pxPerCol, c.msaLength);
      ctx.save();
      ctx.beginPath();
      ctx.rect(ALN_MSA_X, 0, msaW, cssH);
      ctx.clip();
      ctx.font = '11px ui-monospace, monospace';
      ctx.textAlign = letters ? 'center' : 'left';
      for (let i = band.startRow; i < band.endRow; i++) {
        const r = c.rows[i];
        if (!r.seq) continue;
        const y = i * ALN_ROW_H - S;
        const rowGenome = r.genomeRow >= 0 ? c.genomes[r.genomeRow]?.name ?? r.genome : r.genome;
        // leave the selection tint visible through the gaps
        if (!(rowGenome && selected.has(rowGenome))) {
          ctx.fillStyle = '#f3f4f6';
          ctx.fillRect(ALN_MSA_X, y, msaW, ALN_ROW_H - 1);
        }
        for (let col = c0; col < c1; col++) {
          const res = r.seq[col] ?? '-';
          if (res === '-') continue;
          const x = ALN_MSA_X + (col - this.scrollCol) * this.pxPerCol;
          const w = Math.max(1, this.pxPerCol);
          ctx.fillStyle = residueColor(res);
          ctx.fillRect(x, y, w, ALN_ROW_H - 1);
          if (letters) {
            ctx.fillStyle = '#1f2937';
            ctx.fillText(res, x + this.pxPerCol / 2, y + ALN_ROW_H / 2);
          }
        }
      }
      ctx.restore();
      ctx.textAlign = 'left';
    }

    if (c.hasTree) this.drawTree(ctx, cssH, S);
    this.drawGuides(ctx, cssH);
  }

  private drawTree(ctx: CanvasRenderingContext2D, cssH: number, S: number): void {
    const layout = this.treeLayout;
    if (!layout) return;
    const pad = 8;
    const treeW = ALN_TREE_W;
    const xOf = (i: number) =>
      treeNodeX(this.mode, layout.depth[i], layout.maxDepth, layout.edges[i], layout.maxEdges, layout.leafIndex[i] >= 0, treeW, pad);
    const yOf = (i: number) => {
      const leafY = layout.leafIndex[i] >= 0 ? layout.leafIndex[i] : (layout.lo[i] + layout.hi[i]) / 2;
      return leafY * ALN_ROW_H + ALN_ROW_H / 2 - S;
    };
    const margin = ALN_ROW_H * 2;
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, treeW, cssH);
    ctx.clip();
    ctx.lineWidth = 1;
    ctx.strokeStyle = '#475569';
    ctx.beginPath();
    for (let i = 0; i < layout.n; i++) {
      const top = layout.lo[i] * ALN_ROW_H - S;
      const bottom = (layout.hi[i] + 1) * ALN_ROW_H - S;
      if (bottom < -margin || top > cssH + margin) continue;
      const kids = layout.children[i];
      if (kids.length) {
        const x = xOf(i);
        ctx.moveTo(x, yOf(kids[0]));
        ctx.lineTo(x, yOf(kids[kids.length - 1]));
      }
      const p = this.treeData!.parent[i];
      if (p >= 0 && p < i) {
        const y = yOf(i);
        ctx.moveTo(xOf(p), y);
        ctx.lineTo(xOf(i), y);
      }
    }
    ctx.stroke();
    // dotted leaders from tip to the colour column
    ctx.save();
    ctx.setLineDash([2, 3]);
    ctx.strokeStyle = '#cbd5e1';
    ctx.beginPath();
    for (let i = 0; i < layout.n; i++) {
      if (layout.leafIndex[i] < 0) continue;
      const ly = layout.leafIndex[i];
      const y = ly * ALN_ROW_H + ALN_ROW_H / 2 - S;
      if (y < -ALN_ROW_H || y > cssH + ALN_ROW_H) continue;
      ctx.moveTo(xOf(i), y);
      ctx.lineTo(treeW, y);
    }
    ctx.stroke();
    ctx.restore();
    // duplication markers
    for (let i = 0; i < layout.n; i++) {
      if (!this.treeData!.duplication[i]) continue;
      const y = yOf(i);
      if (y < -8 || y > cssH + 8) continue;
      const x = xOf(i);
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
      ctx.lineWidth = 1;
    }
    ctx.restore();
  }

  private drawGuides(ctx: CanvasRenderingContext2D, cssH: number): void {
    ctx.strokeStyle = '#e4e7ec';
    ctx.beginPath();
    ctx.moveTo(ALN_TREE_W + 0.5, 0);
    ctx.lineTo(ALN_TREE_W + 0.5, cssH);
    ctx.moveTo(ALN_TREE_W + ALN_COLOR_W + 0.5, 0);
    ctx.lineTo(ALN_TREE_W + ALN_COLOR_W + 0.5, cssH);
    ctx.moveTo(ALN_MSA_X + 0.5, 0);
    ctx.lineTo(ALN_MSA_X + 0.5, cssH);
    ctx.stroke();
  }

  private drawBar(): void {
    const cssW = this.bar.clientWidth;
    const cssH = this.bar.clientHeight;
    const ctx = prepareCanvas(this.bar, cssW, cssH);
    if (!ctx) return;
    ctx.fillStyle = '#f8fafc';
    ctx.fillRect(0, 0, cssW, cssH);
    ctx.fillStyle = '#475467';
    ctx.font = '10px system-ui, sans-serif';
    ctx.textBaseline = 'middle';
    ctx.fillText('tree', 6, 12);
    ctx.fillText('gene', ALN_TREE_W + ALN_COLOR_W + 4, 12);
    const c = this.content;
    if (!c || !c.conservation) return;
    const msaW = Math.max(1, cssW - ALN_MSA_X);
    const [c0, c1] = visibleCols(this.scrollCol, msaW, this.pxPerCol, c.msaLength);
    const cons = c.conservation;
    const gap = c.gap!;
    for (let col = c0; col < c1; col++) {
      const x = ALN_MSA_X + (col - this.scrollCol) * this.pxPerCol;
      const w = Math.max(1, this.pxPerCol);
      const g = Math.round(60 + 170 * cons[col]);
      ctx.fillStyle = `rgb(${Math.round(g * 0.7)},${g},${Math.round(g * 0.7)})`;
      ctx.fillRect(x, CONS_Y, w, CONS_H);
      if (gap[col] > 0.005) {
        ctx.fillStyle = `rgba(146,64,14,${0.25 + 0.65 * Math.min(1, gap[col])})`;
        ctx.fillRect(x, GAP_Y, w, GAP_H);
      }
    }
    // overview strip
    ctx.fillStyle = '#eef2f6';
    ctx.fillRect(ALN_MSA_X, SB_Y, msaW, SB_H);
    const viewportCols = msaW / this.pxPerCol;
    const handleW = Math.max(8, msaW * Math.min(1, viewportCols / Math.max(1, c.msaLength)));
    const maxScroll = Math.max(0, c.msaLength - viewportCols);
    const frac = maxScroll > 0 ? this.scrollCol / maxScroll : 0;
    ctx.fillStyle = '#94a3b8';
    ctx.fillRect(ALN_MSA_X + Math.min(msaW - handleW, frac * (msaW - handleW)), SB_Y, handleW, SB_H);
    ctx.strokeStyle = '#e4e7ec';
    ctx.beginPath();
    ctx.moveTo(ALN_MSA_X + 0.5, 0);
    ctx.lineTo(ALN_MSA_X + 0.5, cssH);
    ctx.stroke();
  }

  // tree layout kept out of the render loop recomputation
  private treeLayout: ReturnType<typeof layoutGeneTree> | null = null;
  private treeData: GeneTreeArrays | null = null;

  setTree(data: GeneTreeArrays | null): void {
    this.treeData = data;
    this.treeLayout = data ? layoutGeneTree(data) : null;
  }
}
