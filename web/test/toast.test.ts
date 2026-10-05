import {describe, expect, test} from 'bun:test';
import {durationFor, errorIdentity, MAX_TOASTS, pickEviction, resumeDelay, RESUME_MIN_MS, ToastQueue, type ToastRecord} from '../src/toast';

describe('toast queue', () => {
  test('hint dedupes and restarts its timer', () => {
    let t = 0;
    const q = new ToastQueue(() => t);
    const a = q.add('hint', 'zoom in');
    expect(a.restarted).toBe(false);
    t = 1000;
    const b = q.add('hint', 'zoom in');
    expect(b.restarted).toBe(true);
    expect(b.id).toBe(a.id);
    expect(q.list()).toHaveLength(1);
    expect(q.get(a.id)!.deadline).toBe(4000); // restarted at t=1000
    // a different hint stacks
    q.add('hint', 'something else');
    expect(q.list()).toHaveLength(2);
  });

  test('at most 3 toasts; oldest auto-dismissing is evicted first', () => {
    let t = 0;
    const q = new ToastQueue(() => t);
    const a = q.add('notice', 'a');
    q.add('notice', 'b');
    q.add('notice', 'c');
    t = 100;
    const res = q.add('notice', 'd');
    expect(res.evicted).toEqual([a.id]);
    expect(q.list().map((x) => x.text)).toEqual(['b', 'c', 'd']);
    expect(q.list()).toHaveLength(MAX_TOASTS);
  });

  test('sticky errors survive while auto toasts exist; cap still enforced', () => {
    let t = 0;
    const q = new ToastQueue(() => t);
    const e1 = q.add('error', 'e1');
    const e2 = q.add('error', 'e2');
    const n1 = q.add('notice', 'n1');
    // adding a 4th: the oldest auto (n1) goes, not an error
    const res = q.add('notice', 'n2');
    expect(res.evicted).toEqual([n1.id]);
    expect(q.list().map((x) => x.text)).toEqual(['e1', 'e2', 'n2']);
    void e1;
    void e2;
  });

  test('all-sticky set evicts the oldest toast to stay at the cap', () => {
    const q = new ToastQueue();
    const e1 = q.add('error', 'e1');
    q.add('error', 'e2');
    q.add('error', 'e3');
    const res = q.add('error', 'e4');
    expect(res.evicted).toEqual([e1.id]);
    expect(q.list()).toHaveLength(MAX_TOASTS);
  });
});

describe('error dedupe and network collapse', () => {
  test('repeated errors bump the count without evicting or re-inserting', () => {
    const q = new ToastQueue(() => 0);
    const a = q.add('error', 'matrix failed: boom', {key: 'matrix|500|boom', details: {status: 500, message: 'boom'}});
    const b = q.add('error', 'matrix failed: boom', {key: 'matrix|500|boom', details: {status: 500, message: 'boom2'}});
    expect(b.restarted).toBe(true);
    expect(b.id).toBe(a.id);
    expect(b.evicted).toEqual([]);
    expect(q.list()).toHaveLength(1);
    expect(q.get(a.id)!.count).toBe(2);
    expect(q.get(a.id)!.details!.message).toBe('boom2'); // details refreshed in place
  });

  test('network failures from any endpoint collapse to one key', () => {
    const q = new ToastQueue(() => 0);
    const k1 = errorIdentity('matrix', undefined, 'Failed to fetch');
    const k2 = errorIdentity('search', undefined, 'Failed to fetch');
    expect(k1).toBe('network');
    expect(k2).toBe('network');
    q.add('error', 'server unreachable', {key: k1});
    q.add('error', 'server unreachable', {key: k2});
    expect(q.list()).toHaveLength(1);
    expect(q.list()[0].count).toBe(2);
  });

  test('http errors key on context + status + message', () => {
    expect(errorIdentity('matrix', 500, 'boom')).toBe('matrix|500|boom');
    expect(errorIdentity('matrix', 500, 'boom')).not.toBe(errorIdentity('search', 500, 'boom'));
    expect(errorIdentity('matrix', 503, 'busy')).not.toBe(errorIdentity('matrix', 500, 'boom'));
  });

  test('repeated errors do not evict other toasts', () => {
    const q = new ToastQueue(() => 0);
    q.add('error', 'a', {key: 'a'});
    q.add('error', 'b', {key: 'b'});
    q.add('error', 'c', {key: 'c'});
    const res = q.add('error', 'a', {key: 'a'});
    expect(res.evicted).toEqual([]);
    expect(q.list().map((x) => x.count)).toEqual([2, 1, 1]);
  });
});

describe('toast pause / resume timing', () => {
  test('pause freezes the countdown, resume keeps the remainder', () => {
    let t = 0;
    const q = new ToastQueue(() => t);
    const a = q.add('notice', 'x'); // deadline 6000
    t = 1000;
    q.pause(a.id);
    expect(q.get(a.id)!.remaining).toBe(5000);
    expect(q.get(a.id)!.paused).toBe(true);
    t = 5000;
    q.resume(a.id);
    expect(q.get(a.id)!.remaining).toBe(5000);
    expect(q.get(a.id)!.deadline).toBe(10000);
    t = 9999;
    expect(q.tick()).toEqual([]);
    t = 10000;
    expect(q.tick()).toEqual([a.id]);
  });

  test('resume leaves at least the 1.5 s floor', () => {
    let t = 0;
    const q = new ToastQueue(() => t);
    const a = q.add('notice', 'x'); // deadline 6000
    t = 5500;
    q.pause(a.id);
    expect(q.get(a.id)!.remaining).toBe(500);
    t = 6000;
    q.resume(a.id);
    expect(q.get(a.id)!.remaining).toBe(RESUME_MIN_MS);
    expect(q.get(a.id)!.deadline).toBe(7500);
  });

  test('restarting a paused hint keeps it paused', () => {
    let t = 0;
    const q = new ToastQueue(() => t);
    const a = q.add('hint', 'h'); // deadline 3000
    t = 1000;
    q.pause(a.id); // remaining 2000
    t = 1500;
    q.add('hint', 'h'); // restart
    expect(q.get(a.id)!.paused).toBe(true);
    expect(q.get(a.id)!.deadline).toBeNull();
    expect(q.get(a.id)!.remaining).toBe(3000);
    q.resume(a.id);
    expect(q.get(a.id)!.deadline).toBe(4500);
  });

  test('sticky errors never expire', () => {
    let t = 0;
    const q = new ToastQueue(() => t);
    const a = q.add('error', 'boom');
    t = 1_000_000;
    expect(q.tick()).toEqual([]);
    expect(q.list()).toHaveLength(1);
    expect(q.remove(a.id)).toBe(true);
    expect(q.list()).toHaveLength(0);
  });
});

describe('toast helpers', () => {
  test('durations and resume floor', () => {
    expect(durationFor('hint')).toBe(3000);
    expect(durationFor('notice')).toBe(6000);
    expect(durationFor('error')).toBeNull();
    expect(resumeDelay(100)).toBe(1500);
    expect(resumeDelay(4000)).toBe(4000);
  });

  test('pickEviction prefers the oldest auto toast', () => {
    const mk = (kind: 'hint' | 'notice' | 'error', id: number): ToastRecord => ({
      id,
      kind,
      text: '',
      count: 1,
      remaining: kind === 'error' ? null : 1000,
      deadline: 0,
      paused: false,
    });
    expect(pickEviction([mk('error', 1), mk('notice', 2), mk('error', 3)])).toBe(1);
    expect(pickEviction([mk('error', 1), mk('error', 2)])).toBe(0);
  });
});

describe('request error classification', () => {
  test('transport failures become NetworkError, HTTP and payload errors do not', async () => {
    const {getMeta, NetworkError, HttpError} = await import('../src/api');
    const realFetch = globalThis.fetch;
    try {
      (globalThis as any).fetch = () => Promise.reject(new TypeError('Failed to fetch'));
      await expect(getMeta()).rejects.toBeInstanceOf(NetworkError);
      (globalThis as any).fetch = async () => ({ok: false, status: 500, json: async () => ({detail: 'boom'})});
      await expect(getMeta()).rejects.toBeInstanceOf(HttpError);
      (globalThis as any).fetch = async () => ({ok: true, status: 200, json: async () => { throw new SyntaxError('bad json'); }});
      const err = await getMeta().catch((e) => e);
      expect(err).not.toBeInstanceOf(NetworkError);
      expect(err).not.toBeInstanceOf(HttpError);
      const ctl = new AbortController();
      ctl.abort();
      (globalThis as any).fetch = () => Promise.reject(Object.assign(new Error('aborted'), {name: 'AbortError'}));
      const aborted = await getMeta().catch((e) => e); // getMeta has no signal, but name check keeps aborts
      expect(aborted.name).toBe('AbortError');
    } finally {
      globalThis.fetch = realFetch;
    }
  });
});
