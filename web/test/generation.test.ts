import {describe, expect, test} from 'bun:test';
import {Generation} from '../src/generation';

describe('drawer generation / abort', () => {
  test('begin() increments and aborts the previous requests', () => {
    const g = new Generation();
    expect(g.current).toBe(0);
    const a = g.begin();
    expect(a.gen).toBe(1);
    expect(g.isCurrent(1)).toBe(true);
    expect(a.signal.aborted).toBe(false);
    const b = g.begin();
    expect(b.gen).toBe(2);
    expect(a.signal.aborted).toBe(true); // stale request aborted
    expect(b.signal.aborted).toBe(false);
    expect(g.isCurrent(1)).toBe(false);
    expect(g.isCurrent(2)).toBe(true);
  });

  test('abort() stops the current request without incrementing', () => {
    const g = new Generation();
    const a = g.begin();
    g.abort();
    expect(a.signal.aborted).toBe(true);
    expect(g.current).toBe(1);
  });
});
