import {describe, expect, test} from 'bun:test';
import {Generation, Sequence} from '../src/generation';

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

describe('navigation sequence', () => {
  test('a token is only current until the next invalidation', () => {
    const s = new Sequence();
    const a = s.current();
    expect(s.isCurrent(a)).toBe(true);
    const b = s.invalidate();
    expect(b).not.toBe(a);
    expect(s.isCurrent(a)).toBe(false); // late prev/next result discarded
    expect(s.isCurrent(b)).toBe(true);
  });

  test('invalidating on close discards an in-flight step', () => {
    const s = new Sequence();
    let applied = false;
    const token = s.current();
    const finish = () => {
      if (s.isCurrent(token)) applied = true;
    };
    s.invalidate(); // closeDrawer()
    finish();
    expect(applied).toBe(false);
  });
});
