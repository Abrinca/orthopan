// Tiny pending-aware memo cache: concurrent loads of the same key share one
// in-flight promise, so repeated events never fire duplicate requests.

export class PendingCache<V> {
  private cache = new Map<string, V>();
  private pending = new Map<string, Promise<V>>();

  constructor(private limit = 64) {}

  get(key: string): V | undefined {
    return this.cache.get(key);
  }

  load(key: string, loader: () => Promise<V>): Promise<V> {
    const cached = this.cache.get(key);
    if (cached !== undefined) return Promise.resolve(cached);
    const pending = this.pending.get(key);
    if (pending) return pending;
    const p: Promise<V> = loader().then(
      (v) => {
        if (this.pending.get(key) !== p) return v; // forgotten: don't cache
        this.cache.set(key, v);
        this.pending.delete(key);
        while (this.cache.size > this.limit) {
          const oldest = this.cache.keys().next().value as string | undefined;
          if (oldest === undefined) break;
          this.cache.delete(oldest);
        }
        return v;
      },
      (err) => {
        if (this.pending.get(key) === p) this.pending.delete(key);
        throw err;
      },
    );
    this.pending.set(key, p);
    return p;
  }

  /** Drop an in-flight entry (e.g. its request was aborted) so the next load starts afresh. */
  forget(key: string): void {
    this.pending.delete(key);
  }
}
