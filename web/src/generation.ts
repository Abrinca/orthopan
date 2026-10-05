// Drawer request lifecycle: a monotonically increasing generation plus an
// AbortController. Starting a new generation aborts the previous in-flight
// requests; a late response is only committed when its generation is current.

export interface GenerationTicket {
  gen: number;
  signal: AbortSignal;
}

export class Generation {
  private n = 0;
  private ctl: AbortController | null = null;

  /** Invalidate the current generation, abort its requests, return the new one. */
  begin(): GenerationTicket {
    this.ctl?.abort();
    this.n++;
    this.ctl = new AbortController();
    return {gen: this.n, signal: this.ctl.signal};
  }

  get current(): number {
    return this.n;
  }

  isCurrent(gen: number): boolean {
    return gen === this.n;
  }

  abort(): void {
    this.ctl?.abort();
    this.ctl = null;
  }
}

/**
 * A bare monotonic token for guarding async navigation (prev/next). Callers
 * capture `current()` before awaiting and only apply the result if
 * `isCurrent(token)` still holds; `invalidate()` on close/route/order changes.
 */
export class Sequence {
  private n = 0;

  current(): number {
    return this.n;
  }

  invalidate(): number {
    return ++this.n;
  }

  isCurrent(token: number): boolean {
    return token === this.n;
  }
}
