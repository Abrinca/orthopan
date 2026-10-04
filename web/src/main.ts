// OrthoPan viewer: tree, metadata, heatmap, strips, search and detail panel.
import {
  HttpError,
  getAnnotationTrack,
  getGenomes,
  getMatrix,
  getMeta,
  getOrthogroup,
  getOrthogroups,
  getPattern,
  getTracks,
  postOrder,
  search as apiSearch,
  type AnnotationTrackResult,
  type Genome,
  type GenomeMetadataInfo,
  type MatrixResult,
  type Meta,
  type SearchHit,
  type TracksResult,
} from './api';
import {annotationCategory, categoryColor} from './annotations';
import {formatCopiesExact, formatCopiesMean} from './decode';
import {Heatmap, colorForValue, copiesLegend, toImageDataFor} from './matrix';
import {Metadata, META_COL_W} from './metadata';
import {PendingCache} from './pending';
import {SearchBox, hitViewport} from './search';
import {
  binRange,
  bufferAt,
  cellRect,
  clampView,
  colRange,
  computeWindow,
  fitView,
  makeTransform,
  prepareCanvas,
  rowRange,
  type View,
  type Window,
} from './transform';
import {Tracks, fmt3, type Strip} from './tracks';
import {Tree} from './tree';

const $ = <T extends HTMLElement>(id: string): T => document.getElementById(id) as T;

const matrixWrap = $('matrix-wrap');
const treeCanvas = $<HTMLCanvasElement>('tree');
const tracksCanvas = $<HTMLCanvasElement>('tracks-canvas');
const metaCol = $('meta-col');
const metaCanvas = $<HTMLCanvasElement>('meta');
const labelsCol = $('labels-col');
const labelsCanvas = $<HTMLCanvasElement>('labels');
const overlay = $<HTMLCanvasElement>('overlay');
const tooltip = $('tooltip');
const panel = $('panel');
const orderSelect = $<HTMLSelectElement>('order');
const modeSelect = $<HTMLSelectElement>('mode');
const sortBtn = $<HTMLButtonElement>('sortbtn');
const statusEl = $('status');
const noticeEl = $('notice');
const legendEl = $('legend');
const layerList = $('layer-list');
const searchInput = $<HTMLInputElement>('search');
const searchResults = $('search-results');
const searchOrderBtn = $<HTMLButtonElement>('search-order');

let meta: Meta;
let genomes: Genome[] = [];
let tree: Tree;
let heatmap: Heatmap;
let tracks: Tracks;
let metadata: Metadata;
let searchBox: SearchBox;
let order = 'tree';
let mode: 'presence' | 'copies' = 'presence';
let view: View = {cx: 0, cy: 0, zx: 0, zy: 0};
let fitted = true;
let activeTracks: string[] = [];
let activeAnnos: string[] = [];
let metaFields: GenomeMetadataInfo[] = [];
let lastWin: Window | null = null;
let lastMat: MatrixResult | null = null;
let dataOrder: string | null = null;
let lastTracks: TracksResult | null = null;
let lastTracksWorld: {og0: number; og1: number} | null = null;
let lastAnnos: AnnotationTrackResult[] | null = null;
let lastAnnoWorld: {og0: number; og1: number} | null = null;
let hover: {row: number; col: number} | null = null;
let highlight: {pos: number; until: number} | null = null;

// fetch / request lifecycle
let gen = 0;
let inflight: AbortController | null = null;
let fetchTimer = 0;
let hashTimer = 0;
let detailSeq = 0;
let hoverSeq = 0;

let activeOgs: {order: string; og0: number; ids: string[]} | null = null;
const ogCache = new Map<string, string[]>();
const ogAnnoCache = new PendingCache<Record<string, unknown>>(64);

function dims(): {pw: number; ph: number} {
  return {pw: matrixWrap.clientWidth, ph: matrixWrap.clientHeight};
}

const esc = (s: string) =>
  s.replace(/[&<>"]/g, (c) => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c] as string);

const fmtVal = (v: unknown) => (v === null || v === undefined || v === '' ? '—' : String(v));

let noticeTimer = 0;
function notice(msg: string): void {
  noticeEl.textContent = msg;
  window.clearTimeout(noticeTimer);
  noticeTimer = window.setTimeout(() => (noticeEl.textContent = ''), 6000);
}

// --- initialization ---------------------------------------------------------

async function boot(): Promise<void> {
  meta = await getMeta();
  genomes = await getGenomes();
  tree = new Tree(treeCanvas);
  await tree.load();
  tracks = new Tracks(tracksCanvas);
  metadata = new Metadata(metaCanvas);

  metaFields = meta.genome_metadata ?? [];
  const metaW = metaFields.length * META_COL_W;
  metaCol.style.display = metaW ? '' : 'none';
  metaCol.style.flex = `0 0 ${metaW}px`;
  metaCol.style.width = `${metaW}px`;

  activeTracks = pickTracks(meta);
  activeAnnos = (meta.annotations ?? []).filter((a) => a.kind === 'categorical').map((a) => a.name);

  const {pw, ph} = dims();
  view = fitView(meta.n_orthogroups, meta.n_genomes, pw, ph);
  populateOrders(meta);
  populateLayers(meta);
  setupSearch();

  applyHash();

  heatmap = new Heatmap(
    $('deck') as HTMLDivElement,
    (v) => clampView(v, meta.n_orthogroups, meta.n_genomes, dims().pw, dims().ph),
    (v) => {
      view = v;
      afterViewChange();
    },
  );
  heatmap.setView(view);
  syncMode();

  window.addEventListener('resize', onResize);
  new ResizeObserver(onResize).observe(matrixWrap);

  bindInteractions();
  knownDpr = window.devicePixelRatio || 1;
  watchDevicePixelRatio();
  window.clearInterval(dprTimer);
  dprTimer = window.setInterval(checkDevicePixelRatio, 400);
  updateHash();
  scheduleFetch();
}

// A change in devicePixelRatio at the same CSS size does not fire resize or the
// ResizeObserver, so the Canvas2D backing stores would keep the old scale.
// Re-arm a resolution media query and redraw (no refit, no refetch) on change.
// Some environments (notably CDP device-metric emulation) update the query's
// `matches` without dispatching `change`, so a cheap poll is the safety net.
let dprQuery: MediaQueryList | null = null;
let knownDpr = 1;
let dprTimer = 0;
function watchDevicePixelRatio(): void {
  dprQuery?.removeEventListener('change', onDprChange);
  dprQuery = window.matchMedia(`(resolution: ${window.devicePixelRatio || 1}dppx)`);
  dprQuery.addEventListener('change', onDprChange);
}

function onDprChange(): void {
  knownDpr = window.devicePixelRatio || 1;
  watchDevicePixelRatio();
  heatmap?.resyncDpr();
  render();
}

function checkDevicePixelRatio(): void {
  if ((window.devicePixelRatio || 1) !== knownDpr) onDprChange();
}

function pickTracks(m: Meta): string[] {
  const out: string[] = [];
  if (m.tracks.some((t) => t.name === 'copy_number_mean')) out.push('copy_number_mean');
  const hom = m.tracks.find((t) => t.name.includes('homogeneity'));
  if (hom) out.push(hom.name);
  if (!out.length && m.tracks.length) out.push(m.tracks[0].name);
  return out.slice(0, 3);
}

function hasOption(value: string): boolean {
  return [...orderSelect.options].some((o) => o.value === value);
}

function addOption(value: string): void {
  orderSelect.insertAdjacentHTML('beforeend', `<option value="${esc(value)}">${esc(value)}</option>`);
}

function populateOrders(m: Meta): void {
  const opts: string[] = [...m.orders];
  for (const t of m.tracks) opts.push(`track:${t.name}`, `track:${t.name}:desc`);
  orderSelect.innerHTML = opts.map((o) => `<option value="${esc(o)}">${esc(o)}</option>`).join('');
  if (!hasOption(order)) addOption(order);
  orderSelect.value = order;
  orderSelect.onchange = () => setOrder(orderSelect.value);
}

function populateLayers(m: Meta): void {
  const boxes: string[] = [];
  boxes.push('<h4>Tracks</h4>');
  for (const t of m.tracks) {
    boxes.push(
      `<label><input type="checkbox" data-kind="track" value="${esc(t.name)}" ${activeTracks.includes(t.name) ? 'checked' : ''}/> ${esc(t.name)}</label>`,
    );
  }
  const cats = (m.annotations ?? []).filter((a) => a.kind === 'categorical');
  if (cats.length) {
    boxes.push('<h4>Annotations</h4>');
    for (const a of cats) {
      const swatches = (a.categories ?? [])
        .slice(0, 8)
        .map((c, i) => `<span class="swatch" style="background:${categoryColor(i)}" title="${esc(c)}"></span>`)
        .join('');
      boxes.push(
        `<label><input type="checkbox" data-kind="anno" value="${esc(a.name)}" ${activeAnnos.includes(a.name) ? 'checked' : ''}/> ${esc(a.name)} ${swatches}</label>`,
      );
    }
  }
  layerList.innerHTML = boxes.join('');
  layerList.querySelectorAll<HTMLInputElement>('input[type=checkbox]').forEach((cb) => {
    cb.addEventListener('change', () => {
      if (cb.dataset.kind === 'track') {
        activeTracks = cb.checked ? [...activeTracks, cb.value] : activeTracks.filter((n) => n !== cb.value);
      } else {
        activeAnnos = cb.checked ? [...activeAnnos, cb.value] : activeAnnos.filter((n) => n !== cb.value);
      }
      scheduleFetch();
    });
  });
}

function setupSearch(): void {
  searchBox = new SearchBox({
    input: searchInput,
    results: searchResults,
    orderButton: searchOrderBtn,
    getOrder: () => order,
    onSelect: (h) => selectHit(h),
    onOrderByHits: (q) => void orderByHits(q),
    onQueryChange: () => scheduleHash(),
  });
}

function setOrder(o: string): void {
  order = o;
  if (!hasOption(o)) addOption(o);
  orderSelect.value = o;
  searchBox?.invalidate(); // order changed: stale hits can no longer be selected
  activeOgs = null;
  ogCache.clear();
  dataOrder = null;
  lastTracks = null;
  lastTracksWorld = null;
  lastAnnos = null;
  lastAnnoWorld = null;
  hover = null;
  detailSeq++;
  hideTooltip();
  scheduleRender();
  scheduleHash();
  scheduleFetch();
}

function syncMode(): void {
  modeSelect.value = mode;
  setupLegend();
}

// --- resize -----------------------------------------------------------------

function onResize(): void {
  const {pw, ph} = dims();
  if (pw < 2 || ph < 2) return;
  if (fitted) view = fitView(meta.n_orthogroups, meta.n_genomes, pw, ph);
  view = clampView(view, meta.n_orthogroups, meta.n_genomes, pw, ph);
  heatmap.setView(view);
  render();
  scheduleFetch();
  scheduleHash();
}

// --- data fetching ----------------------------------------------------------

function scheduleFetch(): void {
  gen++;
  if (inflight) {
    inflight.abort();
    inflight = null;
  }
  window.clearTimeout(fetchTimer);
  fetchTimer = window.setTimeout(() => void fetchData(), 140);
}

async function fetchData(): Promise<void> {
  if (!meta) return;
  const myGen = gen;
  const snapOrder = order;
  const snapMode = mode;
  const {pw, ph} = dims();
  if (pw < 2 || ph < 2) return;
  const win = computeWindow(view, meta.n_orthogroups, meta.n_genomes, pw, ph);
  const ac = new AbortController();
  inflight = ac;
  statusEl.textContent = 'loading…';
  try {
    const m = await getMatrix({order: snapOrder, value: snapMode, ...win}, ac.signal);
    if (myGen !== gen || ac.signal.aborted) return;
    lastMat = m;
    lastWin = win;
    dataOrder = snapOrder;
    heatmap.setData(toImageDataFor(m.value, m.data, m.width, m.height), [win.og0, win.g1, win.og1, win.g0]);
    statusEl.textContent = `${m.width}×${m.height} cells · ${snapOrder} · ${m.value}`;
    render();

    const [tr, annos] = await Promise.all([
      activeTracks.length
        ? getTracks({order: snapOrder, og0: win.og0, og1: win.og1, w: win.w, names: activeTracks}, ac.signal)
        : Promise.resolve(null),
      Promise.all(
        activeAnnos.map((name) =>
          getAnnotationTrack({order: snapOrder, og0: win.og0, og1: win.og1, w: win.w, name}, ac.signal).catch(() => null),
        ),
      ),
    ]);
    if (myGen !== gen || ac.signal.aborted) return;
    lastTracks = tr;
    lastTracksWorld = tr ? {og0: win.og0, og1: win.og1} : null;
    lastAnnos = annos.filter((a): a is AnnotationTrackResult => a !== null);
    lastAnnoWorld = lastAnnos.length ? {og0: win.og0, og1: win.og1} : null;
    render();
    void ensureOgNames(snapOrder, win, ac.signal, myGen);
  } catch (err) {
    if (ac.signal.aborted || myGen !== gen) return;
    if (err instanceof HttpError && err.status === 404 && snapOrder !== 'tree') {
      notice(`unknown order '${snapOrder}' — using tree`);
      setOrder('tree');
      return;
    }
    statusEl.textContent = `error: ${(err as Error).message}`;
  }
}

async function ensureOgNames(
  snapOrder: string,
  win: Window,
  signal: AbortSignal,
  myGen: number,
): Promise<void> {
  if (win.og1 - win.og0 > 5000) return;
  const key = `${snapOrder}:${win.og0}:${win.og1}`;
  const apply = (ids: string[]) => {
    if (myGen === gen && snapOrder === order) {
      activeOgs = {order: snapOrder, og0: win.og0, ids};
      render();
    }
  };
  const cached = ogCache.get(key);
  if (cached) {
    apply(cached);
    return;
  }
  try {
    const ids = await getOrthogroups(snapOrder, win.og0, win.og1, signal);
    if (signal.aborted) return;
    ogCache.set(key, ids);
    while (ogCache.size > 8) {
      const oldest = ogCache.keys().next().value as string | undefined;
      if (oldest === undefined) break;
      ogCache.delete(oldest);
    }
    apply(ids);
  } catch {
    /* labels are optional */
  }
}

// At most one hover detail request in flight: a request for a different OG aborts
// the previous one, so moving across many uncached OGs never piles up work.
let annoHover: {og: string; ctl: AbortController} | null = null;

async function ensureOgAnnotations(og: string): Promise<Record<string, unknown>> {
  const cached = ogAnnoCache.get(og);
  if (cached !== undefined) return cached;
  if (annoHover && annoHover.og !== og) {
    annoHover.ctl.abort();
    ogAnnoCache.forget(annoHover.og);
  }
  if (!annoHover || annoHover.og !== og) annoHover = {og, ctl: new AbortController()};
  const {signal} = annoHover.ctl;
  return ogAnnoCache.load(og, async () => {
    try {
      const d = await getOrthogroup(og, signal);
      return d.annotations ?? {};
    } catch (err) {
      if (signal.aborted) throw err; // not cached; retried on next hover
      return {};
    } finally {
      if (annoHover?.ctl.signal === signal) annoHover = null;
    }
  });
}

function ogAt(pos: number): string | null {
  if (activeOgs && activeOgs.order === order && pos >= activeOgs.og0 && pos < activeOgs.og0 + activeOgs.ids.length) {
    return activeOgs.ids[pos - activeOgs.og0] ?? null;
  }
  return null;
}

// --- strips -----------------------------------------------------------------

function buildStrips(): Strip[] {
  const strips: Strip[] = [];
  if (lastTracks && lastTracksWorld) {
    for (const name of activeTracks) {
      if (!lastTracks.names.includes(name)) continue;
      strips.push({kind: 'numeric', name, info: meta.tracks.find((t) => t.name === name), result: lastTracks, world: lastTracksWorld});
    }
  }
  if (lastAnnos && lastAnnoWorld) {
    for (const res of lastAnnos) {
      const info = (meta.annotations ?? []).find((a) => a.name === res.name);
      strips.push({kind: 'category', name: res.name, categories: info?.categories ?? [], result: res, world: lastAnnoWorld});
    }
  }
  return strips;
}

// --- rendering --------------------------------------------------------------

function afterViewChange(): void {
  fitted = false;
  scheduleRender();
  scheduleFetch();
  scheduleHash();
}

let renderPending = false;
function scheduleRender(): void {
  if (renderPending) return;
  renderPending = true;
  requestAnimationFrame(() => {
    renderPending = false;
    render();
  });
}

function renderOverlayNow(): void {
  if (!meta) return;
  const {pw, ph} = dims();
  if (pw < 2 || ph < 2) return;
  drawOverlay(makeTransform(view, pw, ph), pw, ph);
}

function hideTooltip(): void {
  tooltip.style.display = 'none';
}

function render(): void {
  if (!meta) return;
  syncLabelColumn();
  const {pw, ph} = dims();
  if (pw < 2 || ph < 2) return;
  const tr = makeTransform(view, pw, ph);
  tree.draw(tr, treeCanvas.clientWidth, treeCanvas.clientHeight);
  stripsForHit = buildStrips();
  tracks.draw(tr, stripsForHit, tracksCanvas.clientWidth, tracksCanvas.clientHeight);
  if (metaFields.length) metadata.draw(tr, genomes, metaFields, metaCanvas.clientWidth, metaCanvas.clientHeight);
  drawLabels(tr, labelsCanvas.clientWidth, labelsCanvas.clientHeight);
  drawOverlay(tr, pw, ph);
  updateSortButton();
}

let stripsForHit: Strip[] = [];

function syncLabelColumn(): void {
  const collapsed = 2 ** view.zy < 11;
  labelsCol.classList.toggle('collapsed', collapsed);
  panel.style.right = `${labelsCol.offsetWidth}px`;
}

function drawCladeTint(ctx: CanvasRenderingContext2D, tr: ReturnType<typeof makeTransform>, pw: number): void {
  const range = tree.selectedRange;
  if (!range) return;
  ctx.fillStyle = 'rgba(37,99,235,0.12)';
  const y0 = tr.yOf(range[0]);
  const y1 = tr.yOf(range[1] + 1);
  ctx.fillRect(0, y0, pw, y1 - y0);
}

function drawLabels(tr: ReturnType<typeof makeTransform>, pw: number, ph: number): void {
  const ctx = prepareCanvas(labelsCanvas, pw, ph);
  if (!ctx) return;
  drawCladeTint(ctx, tr, pw);
  if (tr.sy < 11) return;
  const r0 = Math.max(0, Math.floor(tr.rowAt(0)) - 1);
  const r1 = Math.min(meta.n_genomes - 1, Math.ceil(tr.rowAt(ph)) + 1);
  ctx.font = '11px system-ui, sans-serif';
  ctx.textBaseline = 'middle';
  ctx.fillStyle = '#344054';
  for (let r = r0; r <= r1; r++) {
    const y = tr.yOf(r + 0.5);
    if (y < -8 || y > ph + 8) continue;
    ctx.fillText(genomes[r]?.name ?? String(r), 6, y);
  }
}

function drawOverlay(tr: ReturnType<typeof makeTransform>, pw: number, ph: number): void {
  const ctx = prepareCanvas(overlay, pw, ph);
  if (!ctx) return;

  drawCladeTint(ctx, tr, pw);

  // search-hit highlight
  if (highlight && performance.now() < highlight.until) {
    const x0 = tr.xOf(highlight.pos);
    const x1 = tr.xOf(highlight.pos + 1);
    ctx.fillStyle = 'rgba(37,99,235,0.20)';
    ctx.fillRect(x0, 0, Math.max(1, x1 - x0), ph);
    ctx.strokeStyle = 'rgba(37,99,235,0.9)';
    ctx.strokeRect(x0 + 0.5, 0.5, Math.max(1, x1 - x0), ph - 1);
  }

  if (activeOgs && activeOgs.order === order && tr.sx >= 24) {
    ctx.font = '10px system-ui, sans-serif';
    ctx.fillStyle = '#475467';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    for (let i = 0; i < activeOgs.ids.length; i++) {
      const x = tr.xOf(activeOgs.og0 + i + 0.5);
      if (x < 0 || x > pw) continue;
      ctx.save();
      ctx.translate(x, 2);
      ctx.rotate(Math.PI / 2);
      ctx.fillText(activeOgs.ids[i], 0, 0);
      ctx.restore();
    }
    ctx.textAlign = 'left';
  }

  if (hover && lastWin) {
    const [x, y, w, h] = cellRect(lastWin, tr, hover.col, hover.row);
    ctx.strokeStyle = 'rgba(17,24,39,0.7)';
    ctx.lineWidth = 1;
    ctx.strokeRect(x + 0.5, y + 0.5, w, h);
  }
}

function setupLegend(): void {
  if (mode === 'copies') {
    legendEl.innerHTML = copiesLegend()
      .map((b) => `<span style="background:${b.color}" title="${b.label}">${b.label === 'absent' ? '' : b.label}</span>`)
      .join('');
    legendEl.title = 'Copy number: absent, 1, 2, 3, 4, 5+';
    return;
  }
  const stops: string[] = [];
  const n = 16;
  for (let i = 0; i <= n; i++) {
    const [r, g, b] = colorForValue(Math.round((i / n) * 255));
    stops.push(`rgb(${r},${g},${b}) ${((i / n) * 100).toFixed(0)}%`);
  }
  legendEl.innerHTML = `<span style="width:110px;background:linear-gradient(90deg, ${stops.join(',')})"></span>`;
  legendEl.title = 'Presence fraction: left = absent, right = present';
}

// --- interaction ------------------------------------------------------------

let pointerDown: {x: number; y: number} | null = null;

function bindInteractions(): void {
  matrixWrap.addEventListener(
    'wheel',
    (e) => {
      e.preventDefault();
      const {pw, ph} = dims();
      const rect = matrixWrap.getBoundingClientRect();
      const px = e.clientX - rect.left;
      const py = e.clientY - rect.top;
      const tr = makeTransform(view, pw, ph);
      const wx = tr.ogAt(px);
      const wy = tr.rowAt(py);
      const factor = 2 ** (-e.deltaY * 0.0015);
      const nv: View = {...view};
      if (e.shiftKey) {
        nv.zy = view.zy + Math.log2(factor);
        nv.cy = wy - (py - ph / 2) / 2 ** nv.zy;
      } else {
        nv.zx = view.zx + Math.log2(factor);
        nv.cx = wx - (px - pw / 2) / 2 ** nv.zx;
      }
      view = clampView(nv, meta.n_orthogroups, meta.n_genomes, pw, ph);
      heatmap.setView(view);
      afterViewChange();
    },
    {passive: false},
  );

  matrixWrap.addEventListener('pointerdown', (e) => {
    pointerDown = {x: e.clientX, y: e.clientY};
  });
  matrixWrap.addEventListener('mousemove', (e) => onHover(e));
  matrixWrap.addEventListener('mouseleave', () => {
    hover = null;
    hideTooltip();
    renderOverlayNow();
  });
  matrixWrap.addEventListener('click', (e) => {
    const down = pointerDown;
    pointerDown = null;
    if (!down) return;
    if (Math.hypot(e.clientX - down.x, e.clientY - down.y) > 4) return;
    void onClick(e);
  });

  tracksCanvas.addEventListener('mousemove', (e) => onStripHover(e));
  tracksCanvas.addEventListener('mouseleave', hideTooltip);
  metaCanvas.addEventListener('mousemove', (e) => onMetaHover(e));
  metaCanvas.addEventListener('mouseleave', hideTooltip);

  treeCanvas.addEventListener('click', (e) => {
    const rect = treeCanvas.getBoundingClientRect();
    const idx = tree.hitTest(e.clientX - rect.left, e.clientY - rect.top);
    if (idx !== null) {
      tree.selected = tree.selected === idx ? null : idx;
      render();
    }
  });

  modeSelect.onchange = () => {
    mode = modeSelect.value as 'presence' | 'copies';
    syncMode();
    scheduleHash();
    scheduleFetch();
  };

  sortBtn.onclick = () => void sortByContrast();
}

function placeTooltip(clientX: number, clientY: number): void {
  tooltip.style.left = `${clientX + 14}px`;
  tooltip.style.top = `${clientY + 14}px`;
}

function ogRangeText(a: string, b: string): string {
  return a === b ? a : `${a} – ${b}`;
}

function onStripHover(e: MouseEvent): void {
  if (stripsForHit.length === 0 || !lastWin) return;
  const rect = tracksCanvas.getBoundingClientRect();
  const px = e.clientX - rect.left;
  const py = e.clientY - rect.top;
  const tr = makeTransform(view, dims().pw, dims().ph);
  const hit = tracks.hit(tr, stripsForHit, px, py, tracksCanvas.clientHeight);
  if (!hit) {
    hideTooltip();
    return;
  }
  const {og0, og1} = hit.strip.world;
  const [c0, c1] = binRange(og0, og1, hit.strip.result.bins, hit.index);
  const ogA = ogAt(c0) ?? `#${c0}`;
  const ogB = ogAt(Math.max(c0, c1 - 1)) ?? `#${c1 - 1}`;
  const og = c1 - c0 > 1 ? ogRangeText(ogA, ogB) : ogA;
  let body: string;
  if (hit.strip.kind === 'category') {
    body = `${esc(hit.strip.name)}: <b>${esc(annotationCategory(hit.value, hit.strip.categories))}</b>`;
  } else {
    body = `${esc(hit.strip.name)}: <b>${fmt3(hit.value)}</b>`;
  }
  tooltip.innerHTML = `${body}<div class="tip-dim">${esc(og)}</div>`;
  tooltip.style.display = 'block';
  placeTooltip(e.clientX, e.clientY);
}

function onMetaHover(e: MouseEvent): void {
  if (!metaFields.length || !meta) return;
  const rect = metaCanvas.getBoundingClientRect();
  const px = e.clientX - rect.left;
  const py = e.clientY - rect.top;
  const tr = makeTransform(view, dims().pw, dims().ph);
  const hit = metadata.hit(px, py, tr, meta.n_genomes);
  if (!hit || hit.field >= metaFields.length) {
    hideTooltip();
    return;
  }
  const field = metaFields[hit.field];
  const value = genomes[hit.row]?.metadata?.[field.name] ?? null;
  tooltip.innerHTML = `<b>${esc(genomes[hit.row]?.name ?? '')}</b><br>${esc(field.name)}: <b>${esc(fmtVal(value))}</b>`;
  tooltip.style.display = 'block';
  placeTooltip(e.clientX, e.clientY);
}

function setTooltip(baseHtml: string, annos?: Record<string, unknown>): void {
  let html = baseHtml;
  if (annos) {
    const lines = Object.entries(annos)
      .filter(([, v]) => v !== null && v !== undefined && v !== '')
      .map(([k, v]) => `<div><span class="tip-dim">${esc(k)}</span> ${esc(String(v))}</div>`)
      .join('');
    if (lines) html += lines;
  }
  tooltip.innerHTML = html;
}

function onHover(e: MouseEvent): void {
  if (!lastMat || !lastWin || !meta || dataOrder !== order || lastMat.value !== mode) {
    hover = null;
    hideTooltip();
    renderOverlayNow();
    return;
  }
  const {pw, ph} = dims();
  const rect = matrixWrap.getBoundingClientRect();
  const px = e.clientX - rect.left;
  const py = e.clientY - rect.top;
  const tr = makeTransform(view, pw, ph);
  const hit = bufferAt(lastWin, tr, px, py, meta.n_genomes, meta.n_orthogroups);
  if (!hit) {
    hover = null;
    hideTooltip();
    drawOverlay(tr, pw, ph);
    return;
  }
  hover = {row: hit.row, col: hit.col};
  const value = lastMat.data[hit.row * lastMat.width + hit.col] ?? 0;

  let nameText: string;
  let detail: string;
  if (hit.binnedX || hit.binnedY) {
    const [r0, r1] = rowRange(lastWin, hit.row);
    const [c0, c1] = colRange(lastWin, hit.col);
    const gA = genomes[r0]?.name ?? `row ${r0}`;
    const gB = genomes[Math.max(r0, r1 - 1)]?.name ?? `row ${r1 - 1}`;
    nameText = r1 - r0 > 1 ? `${gA} – ${gB}` : gA;
    const ogA = ogAt(c0) ?? `#${c0}`;
    const ogB = ogAt(Math.max(c0, c1 - 1)) ?? `#${c1 - 1}`;
    detail =
      mode === 'copies'
        ? `${ogRangeText(ogA, ogB)} · mean ${formatCopiesMean(value)} copies`
        : `${ogRangeText(ogA, ogB)} · present ${Math.round((value / 255) * 100)}%`;
  } else {
    nameText = genomes[hit.genome]?.name ?? `row ${hit.genome}`;
    const og = ogAt(hit.og) ?? `#${hit.og}`;
    detail =
      mode === 'copies'
        ? `${og} · ${value > 0 ? `${formatCopiesExact(value)} copies` : 'absent'}`
        : `${og} · ${value > 0 ? 'present' : 'absent'}`;
  }
  const baseHtml = `<b>${esc(nameText)}</b><br>${esc(detail)}`;

  // annotations for the hovered OG at full detail (deduped per OG)
  const og = !hit.binnedX && !hit.binnedY ? ogAt(hit.og) : null;
  const my = ++hoverSeq;
  if (og && (meta.annotations ?? []).length) {
    const cached = ogAnnoCache.get(og);
    if (cached) {
      setTooltip(baseHtml, cached);
    } else {
      setTooltip(baseHtml);
      void ensureOgAnnotations(og).then((a) => {
        if (my === hoverSeq && hover && ogAt(hit.og) === og) setTooltip(baseHtml, a);
      }, () => { /* aborted by a newer hover */ });
    }
  } else {
    setTooltip(baseHtml);
  }
  tooltip.style.display = 'block';
  placeTooltip(e.clientX, e.clientY);
  drawOverlay(tr, pw, ph);
}

async function onClick(e: MouseEvent): Promise<void> {
  if (!lastMat || !lastWin || !meta || dataOrder !== order || lastMat.value !== mode) return;
  const {pw, ph} = dims();
  const rect = matrixWrap.getBoundingClientRect();
  const px = e.clientX - rect.left;
  const py = e.clientY - rect.top;
  const tr = makeTransform(view, pw, ph);
  const hit = bufferAt(lastWin, tr, px, py, meta.n_genomes, meta.n_orthogroups);
  if (!hit) return;
  if (hit.binnedX || hit.binnedY) {
    notice('zoom in until the cell covers a single orthogroup to open it');
    return;
  }
  const my = ++detailSeq;
  const pos = hit.og;
  let og = ogAt(pos);
  if (!og) {
    try {
      const ids = await getOrthogroups(order, pos, pos + 1);
      if (my !== detailSeq) return;
      og = ids[0] ?? null;
    } catch {
      og = null;
    }
  }
  if (my !== detailSeq || !og) return;
  await openPanel(og, my);
}

// --- search / jumps ---------------------------------------------------------

function selectHit(hit: SearchHit): void {
  jumpToPosition(hit.position);
  setHighlight(hit.position);
}

function jumpToPosition(pos: number): void {
  const {pw, ph} = dims();
  let v = hitViewport(pos, view, 40);
  v = clampView(v, meta.n_orthogroups, meta.n_genomes, pw, ph);
  view = v;
  fitted = false;
  heatmap.setView(view);
  afterViewChange();
}

function setHighlight(pos: number): void {
  highlight = {pos, until: performance.now() + 2000};
  scheduleRender();
  window.setTimeout(() => {
    if (highlight && performance.now() >= highlight.until) {
      highlight = null;
      renderOverlayNow();
    }
  }, 2100);
}

async function jumpToOg(id: string): Promise<void> {
  try {
    const res = await apiSearch(id, order, 5);
    const hit = res.hits.find((h) => h.orthogroup === id) ?? res.hits[0];
    if (hit) selectHit(hit);
  } catch {
    /* ignore */
  }
}

async function orderByHits(q: string): Promise<void> {
  try {
    const {order: id} = await postOrder({by: 'search', q});
    setOrder(id);
  } catch (err) {
    notice(`search order failed: ${(err as Error).message}`);
  }
}

// --- panel ------------------------------------------------------------------

async function openPanel(og: string, my: number): Promise<void> {
  panel.classList.add('visible');
  panel.innerHTML = '<button class="close" aria-label="close">×</button><h2>loading…</h2>';
  panel.querySelector('.close')!.addEventListener('click', closePanel);
  try {
    const d = await getOrthogroup(og);
    if (my !== detailSeq) return;
    const trackRows = Object.entries(d.tracks ?? {})
      .map(([k, v]) => `<tr><td>${esc(k)}</td><td>${v === null || v === undefined ? '—' : fmt3(v)}</td></tr>`)
      .join('');
    const annoRows = Object.entries(d.annotations ?? {})
      .map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(fmtVal(v))}</td></tr>`)
      .join('');
    const members = (d.members ?? [])
      .map((m) => `<li><code title="${esc(m.genes.map((g) => g.raw).join(', '))}">${esc(m.genome)}</code> — ${m.genes.map((g) => esc(g.id)).join(', ')}</li>`)
      .join('');
    panel.innerHTML = `
      <button class="close" aria-label="close">×</button>
      <h2>${esc(d.orthogroup)}</h2>
      <table>
        <tr><td>og_idx</td><td>${d.og_idx}</td></tr>
        <tr><td>genomes</td><td>${d.n_genomes}</td></tr>
        <tr><td>genes</td><td>${d.n_genes}</td></tr>
        <tr><td>MRCA</td><td>${esc(d.mrca_label ?? '—')}</td></tr>
        <tr><td>purity</td><td>${d.purity === null ? '—' : fmt3(d.purity)}</td></tr>
        ${d.pattern_id !== undefined ? `<tr><td>pattern</td><td>${d.pattern_id} (${d.n_same_pattern ?? 0} same)</td></tr>` : ''}
        ${annoRows}
        ${trackRows}
      </table>
      <h2>members (${d.members?.length ?? 0})</h2>
      <ul>${members}</ul>
      <div id="pattern-section"></div>`;
    panel.querySelector('.close')!.addEventListener('click', closePanel);
    void loadPattern(d.orthogroup, my);
  } catch (err) {
    if (my !== detailSeq) return;
    panel.innerHTML = `<button class="close">×</button><h2>error</h2><p>${esc((err as Error).message)}</p>`;
    panel.querySelector('.close')!.addEventListener('click', closePanel);
  }
}

async function loadPattern(og: string, my: number): Promise<void> {
  const host = document.getElementById('pattern-section');
  if (!host) return;
  try {
    const pat = await getPattern(og);
    if (my !== detailSeq) return;
    const items = pat.orthogroups
      .slice(0, 200)
      .map((p) => {
        const desc = Object.entries(p.annotations ?? {})
          .map(([k, v]) => `${esc(k)}=${esc(fmtVal(v))}`)
          .join(' ');
        return `<li><span class="linkish" data-og="${esc(p.orthogroup)}">${esc(p.orthogroup)}</span>${desc ? ` — ${desc}` : ''}</li>`;
      })
      .join('');
    host.innerHTML = `
      <h2>same pattern (${pat.total})</h2>
      <button id="order-pattern">Order by pattern</button>
      <ul>${items}</ul>`;
    host.querySelectorAll('.linkish').forEach((el) => {
      el.addEventListener('click', () => void jumpToOg((el as HTMLElement).dataset.og!));
    });
    host.querySelector('#order-pattern')!.addEventListener('click', () => void orderByPattern(og));
  } catch {
    host.innerHTML = '';
  }
}

async function orderByPattern(og: string): Promise<void> {
  try {
    const {order: id} = await postOrder({by: 'pattern', og});
    setOrder(id);
  } catch (err) {
    notice(`pattern order failed: ${(err as Error).message}`);
  }
}

function closePanel(): void {
  detailSeq++;
  panel.classList.remove('visible');
}

function updateSortButton(): void {
  const range = tree.selectedRange;
  sortBtn.classList.toggle('visible', !!range && range[1] > range[0]);
}

async function sortByContrast(): Promise<void> {
  const range = tree.selectedRange;
  if (!range) return;
  const names: string[] = [];
  for (let r = range[0]; r <= range[1]; r++) if (genomes[r]) names.push(genomes[r].name);
  if (!names.length) return;
  sortBtn.disabled = true;
  try {
    const {order: id} = await postOrder({by: 'contrast', genomes: names});
    setOrder(id);
  } catch (err) {
    notice(`contrast sort failed: ${(err as Error).message}`);
  } finally {
    sortBtn.disabled = false;
  }
}

// --- URL hash (order, viewport, mode, search) -------------------------------

function scheduleHash(): void {
  window.clearTimeout(hashTimer);
  hashTimer = window.setTimeout(updateHash, 200);
}

function updateHash(): void {
  const q = searchBox?.query ?? '';
  const h =
    `order=${encodeURIComponent(order)}&cx=${view.cx.toFixed(2)}&cy=${view.cy.toFixed(2)}` +
    `&zx=${view.zx.toFixed(3)}&zy=${view.zy.toFixed(3)}&mode=${mode}` +
    (q ? `&q=${encodeURIComponent(q)}` : '');
  history.replaceState(null, '', `#${h}`);
}

function applyHash(): void {
  const p = new URLSearchParams(location.hash.slice(1));
  const o = p.get('order');
  if (o) {
    order = o;
    if (!hasOption(o)) addOption(o);
  }
  const m = p.get('mode');
  if (m === 'presence' || m === 'copies') mode = m;
  const cx = Number(p.get('cx'));
  const cy = Number(p.get('cy'));
  if (p.has('cx') && p.has('cy') && Number.isFinite(cx) && Number.isFinite(cy)) {
    view = {...view, cx, cy};
    fitted = false;
  }
  const zx = Number(p.get('zx'));
  const zy = Number(p.get('zy'));
  if (p.has('zx') && p.has('zy') && Number.isFinite(zx) && Number.isFinite(zy)) {
    view = {...view, zx, zy};
    fitted = false;
  }
  orderSelect.value = order;
  const q = p.get('q');
  if (q) searchBox.setQuery(q);
}

boot().catch((err) => {
  statusEl.textContent = `fatal: ${(err as Error).message}`;
});
