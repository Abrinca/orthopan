// OrthoPan viewer: tree, metadata, heatmap, strips, search and detail panel.
import {
  HttpError,
  getAnnotationTrack,
  getGenomes,
  getMatrix,
  getMeta,
  getMsa,
  getGeneTree,
  getOrthogroup,
  getOrthogroups,
  getPattern,
  getTracks,
  fastaUrl,
  geneTreeNwkUrl,
  postOrder,
  search as apiSearch,
  type AnnotationTrackResult,
  type Genome,
  type GenomeMetadataInfo,
  type GeneTreeArrays,
  type MatrixResult,
  type Meta,
  type SearchHit,
  type TracksResult,
} from './api';
import {annotationCategory, categoryColor} from './annotations';
import {formatCopiesExact, formatCopiesMean} from './decode';
import {GeneTreeView, MISSING_COLOR, type GeneTreeHover} from './genetree';
import {Generation} from './generation';
import {Heatmap, colorForValue, copiesLegend, toImageDataFor} from './matrix';
import {Metadata, META_COL_W, numericColor} from './metadata';
import {MsaView, type MsaHover} from './msa';
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
const drawerTitle = $('drawer-title');
const drawerClose = $<HTMLButtonElement>('drawer-close');
const drawerExpand = $<HTMLButtonElement>('drawer-expand');
const drawerResize = $('drawer-resize');
const idModeSelect = $<HTMLSelectElement>('id-mode');
const dlFasta = $<HTMLAnchorElement>('dl-fasta');
const dlAligned = $<HTMLAnchorElement>('dl-aligned');
const dlNwk = $<HTMLAnchorElement>('dl-nwk');
const tabMembers = $('tab-members');
const msaCanvas = $<HTMLCanvasElement>('msa-canvas');
const msaStatus = $('msa-status');
const msaBadges = $('msa-badges');
const msaRetry = $<HTMLButtonElement>('msa-retry');
const geneTreeCanvas = $<HTMLCanvasElement>('tree-canvas');
const treeFieldSelect = $<HTMLSelectElement>('tree-field');
const treeModeSelect = $<HTMLSelectElement>('tree-mode');
const treeLegend = $('tree-legend');
const treeStatus = $('tree-status');
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
let hoverSeq = 0;
let clickSeq = 0;

// drawer state
let currentOg: string | null = null;
let currentTab: 'members' | 'msa' | 'gene-tree' = 'members';
let highlightRow: number | null = null;
let msaView: MsaView;
let geneTreeView: GeneTreeView;
let treeField = '';

let activeOgs: {order: string; og0: number; ids: string[]} | null = null;
const ogCache = new Map<string, string[]>();
const ogAnnoCache = new PendingCache<Record<string, unknown>>(64);

// drawer request lifecycle: one generation guards every pane; changing OG or
// closing aborts all in-flight drawer requests and discards late responses.
const drawerGen = new Generation();
let drawerSignal: AbortSignal = new AbortController().signal;
const drawerLoaded = {members: null as string | null, msa: null as string | null, tree: null as string | null};
const treeData = new Map<string, GeneTreeArrays>();
const treeState = new Map<string, 'ok' | 'missing'>();
const treeError = new Map<string, string>();
const treeInflight = new Map<string, Promise<GeneTreeArrays>>();

function beginDrawer(): void {
  drawerSignal = drawerGen.begin().signal;
  treeInflight.clear(); // previous requests were aborted; never reuse them
}

const TREE_CACHE_MAX = 16;

/** Remember a tree result; keeps at most TREE_CACHE_MAX OGs (oldest evicted from all maps). */
function rememberTree(og: string, data: GeneTreeArrays | null, error?: string): void {
  treeState.delete(og); // re-insert to mark as most recent
  treeState.set(og, data ? 'ok' : 'missing');
  if (data) treeData.set(og, data);
  if (error !== undefined) treeError.set(og, error);
  while (treeState.size > TREE_CACHE_MAX) {
    const oldest = treeState.keys().next().value as string;
    treeState.delete(oldest);
    treeData.delete(oldest);
    treeError.delete(oldest);
  }
}

/** Fetch a gene tree once per OG, sharing the in-flight promise. Only a 404 is
 *  cached as "missing"; aborts and transient errors stay retryable. */
function fetchTreeData(og: string): Promise<GeneTreeArrays> {
  let p = treeInflight.get(og);
  if (!p) {
    const req: Promise<GeneTreeArrays> = getGeneTree(og, drawerSignal);
    p = req;
    treeInflight.set(og, req);
    req.then(
      (d) => rememberTree(og, d),
      (err) => {
        if (err instanceof HttpError && err.status === 404) rememberTree(og, null, (err as Error).message);
      },
    ).finally(() => {
      if (treeInflight.get(og) === req) treeInflight.delete(og);
    });
  }
  return p;
}

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
  msaView = new MsaView(msaCanvas, onMsaHover);
  geneTreeView = new GeneTreeView(geneTreeCanvas, onTreeHover, onLeafHover);
  setupDrawer();

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
  msaView?.draw();
  geneTreeView?.draw();
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

  // gene-tree leaf hover highlights that genome's matrix row
  if (highlightRow !== null) {
    ctx.fillStyle = 'rgba(220,38,38,0.14)';
    const y0 = tr.yOf(highlightRow);
    const y1 = tr.yOf(highlightRow + 1);
    ctx.fillRect(0, y0, pw, y1 - y0);
  }

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
  const my = ++clickSeq;
  const pos = hit.og;
  let og = ogAt(pos);
  if (!og) {
    try {
      const ids = await getOrthogroups(order, pos, pos + 1);
      if (my !== clickSeq) return;
      og = ids[0] ?? null;
    } catch {
      og = null;
    }
  }
  if (my !== clickSeq || !og) return;
  openDrawer(og);
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

// --- detail drawer ----------------------------------------------------------

function setupDrawer(): void {
  document.querySelectorAll<HTMLButtonElement>('#panel .tab').forEach((b) => {
    b.addEventListener('click', () => setDrawerTab(b.dataset.tab as DrawerTab));
  });
  drawerClose.addEventListener('click', closeDrawer);
  drawerExpand.addEventListener('click', () => {
    panel.classList.toggle('expanded');
    panel.style.width = '';
    saveDrawerState();
    msaView.draw();
    geneTreeView.draw();
  });
  idModeSelect.addEventListener('change', () => {
    updateDownloads();
    updateNwkLink();
  });
  const fields = meta.genome_metadata ?? [];
  treeFieldSelect.innerHTML =
    '<option value="">genome identity</option>' +
    fields.map((f) => `<option value="${esc(f.name)}">${esc(f.name)} (${esc(f.kind)})</option>`).join('');
  const firstCat = fields.find((f) => f.kind === 'categorical');
  treeField = firstCat ? firstCat.name : '';
  treeFieldSelect.value = treeField;
  renderTreeLegend();
  treeFieldSelect.addEventListener('change', () => {
    treeField = treeFieldSelect.value;
    renderTreeLegend();
    geneTreeView.setSpec(metaFields.find((f) => f.name === treeField));
  });
  treeModeSelect.addEventListener('change', () => {
    geneTreeView.setMode(treeModeSelect.value as 'phylogram' | 'cladogram');
  });
  $('msa-zoom-in').addEventListener('click', () => msaView.zoomBy(1.5));
  $('msa-zoom-out').addEventListener('click', () => msaView.zoomBy(1 / 1.5));
  $('tree-zoom-in').addEventListener('click', () => geneTreeView.zoomBy(1.25));
  $('tree-zoom-out').addEventListener('click', () => geneTreeView.zoomBy(1 / 1.25));
  msaRetry.addEventListener('click', () => {
    if (currentOg) void renderMsa(currentOg);
  });
  // restore persisted drawer width / expanded state
  try {
    const w = Number(localStorage.getItem('orthopan.drawerWidth'));
    if (Number.isFinite(w) && w >= 300) panel.style.width = `${Math.min(window.innerWidth * 0.9, w)}px`;
    if (localStorage.getItem('orthopan.drawerExpanded') === '1') panel.classList.add('expanded');
  } catch {
    /* storage unavailable */
  }
  // cheap horizontal resize of the drawer
  drawerResize.addEventListener('pointerdown', (e) => {
    e.preventDefault();
    const startX = e.clientX;
    const startW = panel.offsetWidth;
    const move = (ev: PointerEvent) => {
      panel.classList.remove('expanded');
      panel.style.width = `${Math.max(300, Math.min(window.innerWidth * 0.9, startW + (startX - ev.clientX)))}px`;
      msaView.draw();
      geneTreeView.draw();
    };
    const up = () => {
      saveDrawerState();
      window.removeEventListener('pointermove', move);
      window.removeEventListener('pointerup', up);
    };
    window.addEventListener('pointermove', move);
    window.addEventListener('pointerup', up);
  });
  // Esc closes the drawer unless an input/select/textarea has focus
  document.addEventListener('keydown', (e) => {
    if (e.key !== 'Escape') return;
    const tag = (document.activeElement as HTMLElement | null)?.tagName;
    if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA') return;
    if (panel.classList.contains('visible')) {
      e.preventDefault();
      closeDrawer();
    }
  });
  // back button / hand-edited hash opens the OG and tab
  window.addEventListener('hashchange', onHashChange);
}

function saveDrawerState(): void {
  try {
    localStorage.setItem('orthopan.drawerWidth', String(panel.offsetWidth));
    localStorage.setItem('orthopan.drawerExpanded', panel.classList.contains('expanded') ? '1' : '0');
  } catch {
    /* storage unavailable */
  }
}

function renderTreeLegend(): void {
  const spec = (meta?.genome_metadata ?? []).find((f) => f.name === treeField);
  let html = '';
  if (!spec || spec.kind === 'categorical') {
    const cats = (spec?.categories ?? []).slice(0, 16);
    html = cats.map((c, i) => `<span class="lg"><i style="background:${categoryColor(i)}"></i>${esc(c)}</span>`).join('');
    if (spec) html += `<span class="lg"><i style="background:${MISSING_COLOR}"></i>missing</span>`;
    if (!html) html = `<span class="lg">genome identity</span>`;
  } else {
    const stops = [0, 0.25, 0.5, 0.75, 1].map((f) => `<i style="background:${numericColor(f)}"></i>`).join('');
    html = `<span class="lg">${fmt3(spec.min ?? 0)}${stops}${fmt3(spec.max ?? 1)}</span><span class="lg"><i style="background:${MISSING_COLOR}"></i>missing</span>`;
  }
  html += `<span class="lg"><i class="dup"></i>duplication</span>`;
  treeLegend.innerHTML = html;
}

type DrawerTab = 'members' | 'msa' | 'gene-tree';

function openDrawer(og: string, tab?: DrawerTab): void {
  const changed = og !== currentOg;
  if (changed) {
    beginDrawer();
    drawerLoaded.members = drawerLoaded.msa = drawerLoaded.tree = null;
    clearDrawerPanes();
    currentOg = og;
    drawerTitle.textContent = og;
    updateDownloads();
    updateNwkLink();
    probeTree(og);
  }
  panel.classList.add('visible');
  setDrawerTab(tab ?? currentTab);
}

function clearDrawerPanes(): void {
  tabMembers.innerHTML = '';
  msaView.setData(null);
  msaStatus.textContent = '';
  msaBadges.innerHTML = '';
  msaRetry.hidden = true;
  geneTreeView.setData(null, genomes, undefined);
  treeStatus.textContent = '';
  highlightRow = null;
  hideTooltip();
  renderOverlayNow();
}

function setDrawerTab(tab: DrawerTab): void {
  currentTab = tab;
  document.querySelectorAll<HTMLButtonElement>('#panel .tab').forEach((b) => b.classList.toggle('active', b.dataset.tab === tab));
  document.querySelectorAll<HTMLElement>('#panel .tab-pane').forEach((p) => p.classList.toggle('active', p.id === `tab-${tab}`));
  scheduleHash();
  const og = currentOg;
  if (!og) return;
  if (tab === 'members') {
    if (drawerLoaded.members !== og) void renderMembers(og);
  } else if (tab === 'msa') {
    if (drawerLoaded.msa !== og) void renderMsa(og);
    else msaView.draw();
  } else {
    if (drawerLoaded.tree !== og) void renderGeneTree(og);
    else geneTreeView.draw();
  }
}

function updateDownloads(): void {
  if (!currentOg) return;
  const ids = idModeSelect.value === 'raw' ? 'raw' : 'display';
  dlFasta.href = fastaUrl(currentOg, false, ids);
  dlAligned.href = fastaUrl(currentOg, true, ids);
  dlNwk.href = geneTreeNwkUrl(currentOg, ids);
}

/** Enable/disable the Newick link from what we already know about the OG. */
function updateNwkLink(): void {
  if (!currentOg) return;
  updateDownloads();
  const enabled = treeState.get(currentOg) !== 'missing';
  dlNwk.classList.toggle('disabled', !enabled);
  dlNwk.setAttribute('aria-disabled', enabled ? 'false' : 'true');
  if (!enabled) dlNwk.removeAttribute('href');
}

/** In the background, find out whether the OG has a gene tree so the Newick
 *  link can be disabled before the user clicks it. */
function probeTree(og: string): void {
  if (treeState.has(og)) {
    updateNwkLink();
    return;
  }
  fetchTreeData(og).then(
    () => {
      if (currentOg === og) updateNwkLink();
    },
    () => {
      if (currentOg === og) updateNwkLink();
    },
  );
}

async function renderMembers(og: string): Promise<void> {
  const my = drawerGen.current;
  tabMembers.innerHTML = '<h2>loading…</h2>';
  try {
    const d = await getOrthogroup(og, drawerSignal);
    if (!drawerGen.isCurrent(my) || og !== currentOg) return;
    const trackRows = Object.entries(d.tracks ?? {})
      .map(([k, v]) => `<tr><td>${esc(k)}</td><td>${v === null || v === undefined ? '—' : fmt3(v)}</td></tr>`)
      .join('');
    const annoRows = Object.entries(d.annotations ?? {})
      .map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(fmtVal(v))}</td></tr>`)
      .join('');
    const members = (d.members ?? [])
      .map(
        (m) =>
          `<li class="genome-group"><b>${esc(m.genome)}</b><ul>` +
          m.genes
            .map(
              (g) =>
                `<li><code title="${esc(g.raw)}">${esc(g.id)}</code> ` +
                `<button class="copyid" data-id="${esc(g.raw)}" title="copy raw id">copy</button></li>`,
            )
            .join('') +
          `</ul></li>`,
      )
      .join('');
    tabMembers.innerHTML = `
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
      <h2>member genomes (${d.members?.length ?? 0})</h2>
      <ul class="genome-list">${members}</ul>
      <div id="pattern-section"></div>`;
    tabMembers.querySelectorAll<HTMLButtonElement>('.copyid').forEach((b) => {
      b.addEventListener('click', () => {
        const id = b.dataset.id ?? '';
        navigator.clipboard?.writeText(id).then(
          () => notice(`copied ${id}`),
          () => notice('copy failed'),
        );
      });
    });
    drawerLoaded.members = og;
    void loadPattern(d.orthogroup, my);
  } catch (err) {
    if (!drawerGen.isCurrent(my) || og !== currentOg) return;
    tabMembers.innerHTML = `<h2>error</h2><p>${esc((err as Error).message)}</p>`;
  }
}

async function renderMsa(og: string): Promise<void> {
  const my = drawerGen.current;
  msaStatus.textContent = 'loading…';
  msaBadges.innerHTML = '';
  msaRetry.hidden = true;
  try {
    const d = await getMsa(og, drawerSignal);
    if (!drawerGen.isCurrent(my) || og !== currentOg) return;
    msaView.setData(d);
    drawerLoaded.msa = og;
    msaStatus.textContent = `${d.rows.length} rows × ${d.length} cols`;
    msaBadges.innerHTML =
      `<span class="badge src">${esc(d.source)}</span>` + (d.trimmed === true ? ' <span class="badge">trimmed</span>' : '');
  } catch (err) {
    if (!drawerGen.isCurrent(my) || og !== currentOg) return;
    msaView.setData(null);
    msaStatus.textContent = (err as Error).message;
    if (err instanceof HttpError && err.status === 503) msaRetry.hidden = false;
  }
}

async function renderGeneTree(og: string): Promise<void> {
  const my = drawerGen.current;
  const spec = metaFields.find((f) => f.name === treeField);
  if (treeState.get(og) === 'missing') {
    geneTreeView.setData(null, genomes, spec);
    treeStatus.textContent = treeError.get(og) ?? 'no gene tree for this orthogroup';
    drawerLoaded.tree = og;
    return;
  }
  treeStatus.textContent = 'loading…';
  try {
    const d = treeData.get(og) ?? (await fetchTreeData(og));
    if (!drawerGen.isCurrent(my) || og !== currentOg) return;
    geneTreeView.setData(d, genomes, spec);
    treeStatus.textContent = `${d.parent.length} nodes`;
    drawerLoaded.tree = og;
    updateNwkLink();
  } catch (err) {
    if (!drawerGen.isCurrent(my) || og !== currentOg) return;
    geneTreeView.setData(null, genomes, spec);
    treeStatus.textContent = (err as Error).message;
    if (treeState.get(og) === 'missing') updateNwkLink();
  }
}

function pct(x: number): string {
  return `${Math.round(x * 100)}%`;
}

function onMsaHover(h: MsaHover | null, clientX: number, clientY: number): void {
  if (!h) {
    hideTooltip();
    return;
  }
  if (h.bar) {
    tooltip.innerHTML = `col ${h.col + 1} · conservation ${pct(h.conservation)} · gap ${pct(h.gap)}`;
    tooltip.style.display = 'block';
    placeTooltip(clientX, clientY);
    return;
  }
  tooltip.innerHTML =
    `<b>${esc(h.gene)}</b> · <span class="tip-dim">${esc(h.genome)}</span><br>` +
    `col ${h.col + 1} · ${esc(h.residue)}${h.ungapped > 0 ? ` · ungapped ${h.ungapped}` : ''}<br>` +
    `conservation ${pct(h.conservation)} · gap ${pct(h.gap)}`;
  tooltip.style.display = 'block';
  placeTooltip(clientX, clientY);
}

function onTreeHover(h: GeneTreeHover | null, clientX: number, clientY: number): void {
  if (!h) {
    hideTooltip();
    return;
  }
  const gene = h.gene ?? 'internal node';
  const support = h.support === null || h.support === undefined ? '' : ` · support ${h.support}`;
  tooltip.innerHTML = `<b>${esc(gene)}</b>${h.genome ? ` · <span class="tip-dim">${esc(h.genome)}</span>` : ''}${h.duplication ? ' · <b>duplication</b>' : ''}${support}`;
  tooltip.style.display = 'block';
  placeTooltip(clientX, clientY);
}

function onLeafHover(genomeRow: number | null): void {
  highlightRow = genomeRow;
  renderOverlayNow();
}

function onHashChange(): void {
  const p = new URLSearchParams(location.hash.slice(1));
  const og = p.get('og');
  const tab = p.get('tab');
  if (og && (tab === 'members' || tab === 'msa' || tab === 'gene-tree')) openDrawer(og, tab);
  else if (!og && panel.classList.contains('visible')) closeDrawer();
}

function closeDrawer(): void {
  beginDrawer();
  panel.classList.remove('visible');
  currentOg = null;
  drawerLoaded.members = drawerLoaded.msa = drawerLoaded.tree = null;
  clearDrawerPanes();
  scheduleHash();
}

async function loadPattern(og: string, my: number): Promise<void> {
  const host = document.getElementById('pattern-section');
  if (!host) return;
  try {
    const pat = await getPattern(og, 200, drawerSignal);
    if (!drawerGen.isCurrent(my) || og !== currentOg) return;
    const others = pat.orthogroups.filter((p) => p.orthogroup !== og);
    const items = others
      .slice(0, 200)
      .map((p) => {
        const desc = Object.entries(p.annotations ?? {})
          .map(([k, v]) => `${esc(k)}=${esc(fmtVal(v))}`)
          .join(' ');
        return `<li><span class="linkish" data-og="${esc(p.orthogroup)}">${esc(p.orthogroup)}</span>${desc ? ` — ${desc}` : ''}</li>`;
      })
      .join('');
    host.innerHTML = others.length
      ? `<h2>same pattern (${others.length})</h2><button id="order-pattern">Order by pattern</button><ul>${items}</ul>`
      : `<h2>same pattern</h2><p class="muted">no other orthogroup has this pattern</p><button id="order-pattern">Order by pattern</button>`;
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
  const drawer = currentOg ? `&og=${encodeURIComponent(currentOg)}&tab=${currentTab}` : '';
  const h =
    `order=${encodeURIComponent(order)}&cx=${view.cx.toFixed(2)}&cy=${view.cy.toFixed(2)}` +
    `&zx=${view.zx.toFixed(3)}&zy=${view.zy.toFixed(3)}&mode=${mode}` +
    drawer +
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
  const og = p.get('og');
  const tab = p.get('tab');
  if (og && (tab === 'members' || tab === 'msa' || tab === 'gene-tree')) openDrawer(og, tab);
}

boot().catch((err) => {
  statusEl.textContent = `fatal: ${(err as Error).message}`;
});
