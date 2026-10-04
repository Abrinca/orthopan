import {describe, expect, test} from 'bun:test';
import {PendingCache} from '../src/pending';

describe('PendingCache', () => {
  test('repeated loads of one key share a single in-flight request', async () => {
    const cache = new PendingCache<number>(4);
    let calls = 0;
    let resolve!: (v: number) => void;
    const loader = () => {
      calls++;
      return new Promise<number>((r) => (resolve = r));
    };

    // simulate repeated mousemoves over one uncached OG
    const p1 = cache.load('OG1', loader);
    const p2 = cache.load('OG1', loader);
    const p3 = cache.load('OG1', loader);
    expect(calls).toBe(1);

    resolve(7);
    expect(await p1).toBe(7);
    expect(await p2).toBe(7);
    expect(await p3).toBe(7);
    expect(calls).toBe(1);

    // after completion it is cached
    expect(cache.get('OG1')).toBe(7);
    expect(await cache.load('OG1', loader)).toBe(7);
    expect(calls).toBe(1);
  });

  test('a failed load is dropped so it can be retried', async () => {
    const cache = new PendingCache<number>(4);
    let calls = 0;
    await expect(
      cache.load('x', async () => {
        calls++;
        throw new Error('boom');
      }),
    ).rejects.toThrow('boom');
    expect(calls).toBe(1);
    expect(await cache.load('x', async () => 1)).toBe(1);
    expect(calls).toBe(1);
  });

  test('forget() lets A -> B -> A start a fresh request for A', async () => {
    const cache = new PendingCache<string>(4);
    const resolvers: Record<string, (v: string) => void> = {};
    let calls = 0;
    const loader = (k: string) => () => {
      calls++;
      return new Promise<string>((r) => (resolvers[k + calls] = r));
    };
    const first = cache.load('A', loader('A')); // A1, will be abandoned
    cache.forget('A');
    void cache.load('B', loader('B'));
    const again = cache.load('A', loader('A')); // must not reuse A1
    expect(again).not.toBe(first);
    expect(calls).toBe(3);
    resolvers['A1']('stale');
    await first;
    expect(cache.get('A')).toBeUndefined(); // forgotten result is not cached
    resolvers['A3']('fresh');
    expect(await again).toBe('fresh');
    expect(cache.get('A')).toBe('fresh');
  });
});

