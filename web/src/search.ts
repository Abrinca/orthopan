// Search box: debounced /api/search with a hits dropdown, plus the pure
// viewport math used to centre a hit.

import {search as apiSearch, type SearchHit} from './api';
import type {View} from './transform';

/**
 * View state that centres orthogroup position `pos` on the X axis at a column
 * width where its label is shown (default 40 CSS px/OG). The caller clamps it.
 */
export function hitViewport(pos: number, base: View, targetPxPerOG = 40): View {
  return {...base, cx: pos + 0.5, zx: Math.log2(targetPxPerOG)};
}

export interface SearchBoxOptions {
  input: HTMLInputElement;
  results: HTMLElement;
  orderButton: HTMLButtonElement;
  getOrder: () => string;
  onSelect: (hit: SearchHit) => void;
  onOrderByHits: (q: string) => void;
  /** called on user query changes (typing, clear, Escape), not on restore */
  onQueryChange?: () => void;
}

export class SearchBox {
  private timer = 0;
  private seq = 0;
  private hits: SearchHit[] = [];
  private controller: AbortController | null = null;

  constructor(private o: SearchBoxOptions) {
    o.input.addEventListener('input', () => {
      this.invalidate(); // query changed: abort in-flight, hide stale results
      this.schedule();
      this.o.onQueryChange?.();
    });
    o.input.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') {
        this.clear();
        this.o.onQueryChange?.();
      }
      if (e.key === 'Enter') this.run();
    });
    o.orderButton.addEventListener('click', () => {
      const q = this.query;
      if (q) this.o.onOrderByHits(q);
    });
    document.addEventListener('click', (e) => {
      if (!o.results.contains(e.target as Node) && e.target !== o.input) this.hide();
    });
  }

  get query(): string {
    return this.o.input.value.trim();
  }

  /** Restore from the URL hash: set the value and search without a hash echo. */
  setQuery(q: string): void {
    this.o.input.value = q;
    if (q) void this.run();
  }

  schedule(): void {
    window.clearTimeout(this.timer);
    this.timer = window.setTimeout(() => void this.run(), 220);
  }

  /**
   * Abort any in-flight/pending search and drop its results. Used on a query or
   * order change so a late response can never render or be selected.
   */
  invalidate(): void {
    window.clearTimeout(this.timer);
    this.seq++;
    this.controller?.abort();
    this.controller = null;
    this.hits = [];
    this.hide();
  }

  async run(): Promise<void> {
    window.clearTimeout(this.timer);
    this.controller?.abort();
    this.controller = null;
    const q = this.query;
    const my = ++this.seq;
    if (!q) {
      this.hits = [];
      this.hide();
      return;
    }
    const order = this.o.getOrder();
    const ac = new AbortController();
    this.controller = ac;
    try {
      const res = await apiSearch(q, order, 50, ac.signal);
      if (my !== this.seq || ac.signal.aborted) return;
      // query/order may have changed while the request was in flight
      if (q !== this.query || order !== this.o.getOrder()) return;
      this.hits = res.hits;
      this.render(res.total, my);
    } catch {
      /* aborted or best-effort */
    }
  }

  private render(total: number, my: number): void {
    const items = this.hits
      .slice(0, 50)
      .map((h, i) => {
        const text = h.text ? `<span class="snippet">${esc(h.text.slice(0, 80))}</span>` : '';
        return `<li data-i="${i}"><b>${esc(h.orthogroup)}</b> <span class="field">${esc(h.field)}</span> ${text}</li>`;
      })
      .join('');
    this.o.results.innerHTML = `<div class="total">${total} hit${total === 1 ? '' : 's'}</div><ul>${items || '<li class="empty">no hits</li>'}</ul>`;
    this.o.results.classList.add('visible');
    this.o.results.querySelectorAll('li[data-i]').forEach((el) => {
      el.addEventListener('click', () => {
        if (my !== this.seq) return; // stale result
        const hit = this.hits[Number((el as HTMLElement).dataset.i)];
        if (hit) {
          this.hide();
          this.o.onSelect(hit);
        }
      });
    });
  }

  hide(): void {
    this.o.results.classList.remove('visible');
  }

  clear(): void {
    this.invalidate();
    this.o.input.value = '';
  }
}

const esc = (s: string) =>
  s.replace(/[&<>"]/g, (c) => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c] as string);
