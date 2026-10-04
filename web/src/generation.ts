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
