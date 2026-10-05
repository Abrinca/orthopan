// Toast notifications: a small transient message queue plus a bottom-centre
// container. The queue (dedupe, max-3 eviction, pause/resume timing) is pure and
// clock-injected so it can be unit-tested; the container does the DOM work.

export type ToastKind = 'hint' | 'notice' | 'error';

export interface ToastAction {
  label: string;
  onClick: () => void;
}

export interface ToastDetails {
  status?: number;
  url?: string;
  message?: string;
}

export interface ToastRecord {
  id: number;
  kind: ToastKind;
  text: string;
  action?: ToastAction;
  details?: ToastDetails;
  /** stable identity for error dedupe (context|status|message, or 'network') */
  key?: string;
  /** number of times this (deduped) toast fired; shown as "×N" when > 1 */
  count: number;
  /** ms left before auto-dismiss, or null when sticky */
  remaining: number | null;
  /** absolute time of auto-dismiss, or null when paused/sticky */
  deadline: number | null;
  paused: boolean;
}

/** Stable identity for an error toast. Network failures (no status) collapse. */
export function errorIdentity(context: string | undefined, status: number | undefined, message: string): string {
  return status === undefined ? 'network' : `${context ?? ''}|${status}|${message}`;
}

export const TOAST_DURATIONS: Record<ToastKind, number | null> = {
  hint: 3000,
  notice: 6000,
  error: null, // sticky until closed
};
export const MAX_TOASTS = 3;
/** When resuming after hover, leave at least this much time. */
export const RESUME_MIN_MS = 1500;

export function durationFor(kind: ToastKind): number | null {
  return TOAST_DURATIONS[kind];
}

/** Remaining time after a hover pause: the remainder, but at least 1.5 s. */
export function resumeDelay(remaining: number): number {
  return Math.max(remaining, RESUME_MIN_MS);
}

/** Index of the toast to evict: the oldest auto-dismissing one, else the oldest. */
export function pickEviction(items: readonly ToastRecord[]): number {
  const auto = items.findIndex((t) => t.remaining !== null);
  return auto >= 0 ? auto : 0;
}

export interface AddResult {
  id: number;
  /** true when an existing hint was restarted instead of adding a duplicate */
  restarted: boolean;
  /** ids of toasts evicted to keep the visible count at the cap */
  evicted: number[];
}

/** Pure toast queue with an injected clock. */
export class ToastQueue {
  private items: ToastRecord[] = [];
  private nextId = 1;

  constructor(
    private now: () => number = () => Date.now(),
    private maxVisible = MAX_TOASTS,
  ) {}

  list(): readonly ToastRecord[] {
    return this.items;
  }

  get(id: number): ToastRecord | undefined {
    return this.items.find((t) => t.id === id);
  }

  add(
    kind: ToastKind,
    text: string,
    opts: {action?: ToastAction; details?: ToastDetails; key?: string} = {},
  ): AddResult {
    const duration = durationFor(kind);
    if (kind === 'hint') {
      const existing = this.items.find((t) => t.kind === 'hint' && t.text === text);
      if (existing) {
        existing.remaining = duration;
        if (!existing.paused) existing.deadline = this.now() + (duration ?? 0);
        return {id: existing.id, restarted: true, evicted: []};
      }
    }
    if (kind === 'error' && opts.key) {
      const existing = this.items.find((t) => t.kind === 'error' && t.key === opts.key);
      if (existing) {
        // bump the count in place: keep the element (and its expanded details)
        existing.count += 1;
        if (opts.details) existing.details = opts.details;
        return {id: existing.id, restarted: true, evicted: []};
      }
    }
    const rec: ToastRecord = {
      id: this.nextId++,
      kind,
      text,
      action: opts.action,
      details: opts.details,
      key: opts.key,
      count: 1,
      remaining: duration,
      deadline: duration !== null ? this.now() + duration : null,
      paused: false,
    };
    this.items.push(rec);
    const evicted: number[] = [];
    while (this.items.length > this.maxVisible) {
      const idx = pickEviction(this.items);
      const [removed] = this.items.splice(idx, 1);
      evicted.push(removed.id);
    }
    return {id: rec.id, restarted: false, evicted};
  }

  /** Freeze the auto-dismiss countdown while hovered. */
  pause(id: number): void {
    const t = this.get(id);
    if (!t || t.remaining === null || t.paused) return;
    t.remaining = Math.max(0, (t.deadline ?? this.now()) - this.now());
    t.paused = true;
    t.deadline = null;
  }

  /** Resume the countdown, leaving at least RESUME_MIN_MS. */
  resume(id: number): void {
    const t = this.get(id);
    if (!t || t.remaining === null || !t.paused) return;
    t.remaining = resumeDelay(t.remaining);
    t.paused = false;
    t.deadline = this.now() + t.remaining;
  }

  remove(id: number): boolean {
    const i = this.items.findIndex((t) => t.id === id);
    if (i < 0) return false;
    this.items.splice(i, 1);
    return true;
  }

  /** Drop expired auto-dismissing toasts; returns their ids. */
  tick(): number[] {
    const now = this.now();
    const expired: number[] = [];
    this.items = this.items.filter((t) => {
      if (t.remaining === null || t.paused || t.deadline === null) return true;
      if (now >= t.deadline) {
        expired.push(t.id);
        return false;
      }
      return true;
    });
    return expired;
  }
}

// --- DOM container ----------------------------------------------------------

function detailsLines(d: ToastDetails): string {
  const lines: string[] = [];
  if (d.status !== undefined) lines.push(`status: ${d.status}`);
  if (d.url) lines.push(`url: ${d.url}`);
  if (d.message) lines.push(`message: ${d.message}`);
  return lines.join('\n');
}

/** Body text including the "×N" repeat count when deduped. */
function displayText(rec: ToastRecord): string {
  return rec.count > 1 ? `${rec.text} ×${rec.count}` : rec.text;
}

export class ToastContainer {
  private queue: ToastQueue;
  private els = new Map<number, HTMLElement>();
  private timer: number;

  constructor(
    private root: HTMLElement,
    queue?: ToastQueue,
  ) {
    this.queue = queue ?? new ToastQueue();
    this.timer = window.setInterval(() => this.tick(), 250);
  }

  hint(text: string): number {
    const res = this.queue.add('hint', text);
    this.apply(res);
    return res.id;
  }

  notice(text: string, opts: {action?: ToastAction; details?: ToastDetails} = {}): number {
    const res = this.queue.add('notice', text, opts);
    this.apply(res);
    return res.id;
  }

  error(text: string, opts: {details?: ToastDetails; key?: string} = {}): number {
    const res = this.queue.add('error', text, opts);
    this.apply(res);
    return res.id;
  }

  /** Dismiss a toast by id (no-op if already gone). */
  dismiss(id: number): void {
    this.queue.remove(id);
    this.removeEl(id);
  }

  /** Tear down the polling timer (mainly for tests). */
  destroy(): void {
    window.clearInterval(this.timer);
  }

  private apply(res: AddResult): void {
    for (const id of res.evicted) this.removeEl(id);
    const existing = this.els.get(res.id);
    if (existing) {
      // a deduped error (or restarted hint): update in place, keep expanded state
      const rec = this.queue.get(res.id);
      if (rec) this.updateEl(existing, rec);
      return;
    }
    const rec = this.queue.get(res.id);
    if (!rec) return;
    const el = this.render(rec);
    this.els.set(rec.id, el);
    this.root.appendChild(el);
  }

  private updateEl(el: HTMLElement, rec: ToastRecord): void {
    const body = el.querySelector('.toast-body');
    if (body) body.textContent = displayText(rec);
    if (rec.details) {
      const det = el.querySelector('.toast-details');
      if (det) det.textContent = detailsLines(rec.details);
    }
  }

  private removeEl(id: number): void {
    this.els.get(id)?.remove();
    this.els.delete(id);
  }

  private tick(): void {
    for (const id of this.queue.tick()) this.removeEl(id);
  }

  private render(rec: ToastRecord): HTMLElement {
    const el = document.createElement('div');
    el.className = `toast toast-${rec.kind}`;
    if (rec.kind === 'error') {
      el.setAttribute('role', 'alert');
    } else {
      el.setAttribute('role', 'status');
      el.setAttribute('aria-live', 'polite');
    }

    const row = document.createElement('div');
    row.className = 'toast-row';

    const body = document.createElement('div');
    body.className = 'toast-body';
    body.textContent = displayText(rec); // textContent escapes
    row.appendChild(body);

    if (rec.action) {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'toast-action';
      btn.textContent = rec.action.label;
      const action = rec.action;
      btn.addEventListener('click', () => {
        this.dismiss(rec.id);
        action.onClick();
      });
      row.appendChild(btn);
    }

    let det: HTMLDivElement | null = null;
    if (rec.kind === 'error' && rec.details) {
      const detId = `toast-details-${rec.id}`;
      const toggle = document.createElement('button');
      toggle.type = 'button';
      toggle.className = 'toast-details-toggle';
      toggle.textContent = 'details';
      toggle.setAttribute('aria-expanded', 'false');
      toggle.setAttribute('aria-controls', detId);
      row.appendChild(toggle);
      det = document.createElement('div');
      det.className = 'toast-details';
      det.id = detId;
      det.hidden = true;
      det.textContent = detailsLines(rec.details);
      toggle.addEventListener('click', () => {
        const open = det!.hidden;
        det!.hidden = !open;
        toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
      });
    }

    if (rec.kind === 'error') {
      const close = document.createElement('button');
      close.type = 'button';
      close.className = 'toast-close';
      close.setAttribute('aria-label', 'close');
      close.textContent = '×';
      close.addEventListener('click', (e) => {
        e.stopPropagation();
        this.dismiss(rec.id);
      });
      row.appendChild(close);
    }
    el.appendChild(row);
    if (det) el.appendChild(det);

    el.addEventListener('mouseenter', () => this.queue.pause(rec.id));
    el.addEventListener('mouseleave', () => this.queue.resume(rec.id));
    return el;
  }
}

let singleton: ToastContainer | null = null;

/** Lazily-created app-wide toast container (bottom centre, above overlays). */
export function toasts(): ToastContainer {
  if (!singleton) {
    const root = document.createElement('div');
    root.id = 'toast-container';
    document.body.appendChild(root);
    singleton = new ToastContainer(root);
  }
  return singleton;
}
