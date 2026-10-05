import {describe, expect, test} from 'bun:test';
import {SearchBox} from '../src/search';

// --- minimal DOM/fetch fakes ------------------------------------------------

function fakeInput() {
  const listeners: Record<string, ((e?: unknown) => void)[]> = {};
  let value = '';
  return {
    get value() {
      return value;
    },
    set value(v: string) {
      value = v;
    },
    addEventListener(t: string, f: (e?: unknown) => void) {
      (listeners[t] ??= []).push(f);
    },
    closest: () => null,
    getBoundingClientRect: () => ({left: 0, bottom: 20}),
    fire(t: string, e: unknown = {}) {
      (listeners[t] ?? []).forEach((f) => f(e));
    },
  };
}

function fakeResults() {
  const classes = new Set<string>();
  return {
    innerHTML: '',
    classList: {
      add: (c: string) => classes.add(c),
      remove: (c: string) => classes.delete(c),
      contains: (c: string) => classes.has(c),
    },
    classes,
    contains: () => false,
    querySelectorAll: () => [],
    style: {} as Record<string, string>,
    offsetWidth: 300,
  };
}

interface Call {
  resolve: (v: unknown) => void;
}
const tick = () => new Promise((r) => setTimeout(r, 0));

function setup(getOrder: () => string) {
  (globalThis as any).window = globalThis;
  (globalThis as any).document = {addEventListener() {}};
  const calls: Call[] = [];
  (globalThis as any).fetch = () =>
    new Promise((resolve) => {
      calls.push({resolve});
    });
  const input = fakeInput();
  const results = fakeResults();
  const selected: unknown[] = [];
  const box = new SearchBox({
    input: input as unknown as HTMLInputElement,
    results: results as unknown as HTMLElement,
    orderButton: {addEventListener() {}} as unknown as HTMLButtonElement,
    getOrder,
    onSelect: (h) => selected.push(h),
    onOrderByHits() {},
  });
  return {box, input, results, calls, selected};
}

const hitsResponse = (ids: string[]) => ({
  ok: true,
  json: async () => ({hits: ids.map((orthogroup, i) => ({orthogroup, og_idx: i, position: i, field: 'orthogroup', text: orthogroup})), total: ids.length}),
});

describe('SearchBox staleness', () => {
  test('a delayed older response cannot overwrite a newer query', async () => {
    const {box, input, results, calls} = setup(() => 'tree');

    input.value = 'a';
    input.fire('input');
    void box.run(); // request 0
    input.value = 'b';
    input.fire('input');
    void box.run(); // request 1 (query changed -> seq bumped)

    expect(calls.length).toBe(2);
    calls[1].resolve(hitsResponse(['FRESH']));
    await tick();
    calls[0].resolve(hitsResponse(['STALE']));
    await tick();

    expect(results.innerHTML).toContain('FRESH');
    expect(results.innerHTML).not.toContain('STALE');
  });

  test('an order change invalidates an in-flight response', async () => {
    let order = 'tree';
    const {box, input, results, calls} = setup(() => order);

    input.value = 'a';
    input.fire('input');
    void box.run();
    expect(calls.length).toBe(1);

    order = 'frequency';
    box.invalidate(); // setOrder() does this
    calls[0].resolve(hitsResponse(['STALE']));
    await tick();

    expect(results.innerHTML).not.toContain('STALE');
    expect(results.classes.has('visible')).toBe(false);
  });

  test('an empty query bumps the sequence and fires no request', async () => {
    const {box, input, calls} = setup(() => 'tree');
    input.value = 'a';
    input.fire('input');
    void box.run();
    expect(calls.length).toBe(1);

    input.value = '';
    input.fire('input');
    void box.run();
    expect(calls.length).toBe(1); // empty query does not fetch
  });

  test('clear() cancels the pending debounce', async () => {
    const {box, input, calls} = setup(() => 'tree');
    input.value = 'a';
    input.fire('input'); // schedules a 220ms debounce
    box.clear();
    await new Promise((r) => setTimeout(r, 260));
    expect(calls.length).toBe(0); // the debounce never fired
  });
});
