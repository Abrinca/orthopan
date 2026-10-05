// OrthoPan viewer: tree, metadata, heatmap, strips, search and detail panel.
import {
  HttpError,
  NetworkError,
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
  type MsaResult,
  type SearchHit,
  type TracksResult,
} from './api';
import {annotationCategory, categoryColor} from './annotations';
import {AlignmentView, orderAlignmentRows, locateOg, stepPosition, type AlnContent, type AlnHover, type AlnBarHover} from './alignment';
import {formatCopiesExact, formatCopiesMean} from './decode';
import {Generation, Sequence} from './generation';
import {Heatmap, colorForValue, copiesLegend, toImageDataFor, FLAT_TILE, tileParams as tileParamsFor, type TileParams} from './matrix';
import {Metadata, META_COL_W, numericColor} from './metadata';
import {toasts, errorIdentity} from './toast';
import {computeConservation} from './msa';
import {memberSlice, type MemberGroup} from './members';
import {PendingCache} from './pending';
import {needsFetch, requestKeys} from './requestkey';
import {decodeSelection, encodeSelection, rangeSelection, selectionNames, singleSelection, summarizeSelection, type GenomeSelection} from './selection';
import {currentTheme, getThemeSetting, initTheme, onThemeChange, palette, setThemeSetting, type ThemeSetting} from './theme';
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
  zoomedView,
  type View,
  type Window,
} from './transform';
import {Tracks, fmt3, type Strip} from './tracks';
import {Tree} from './tree';
import {placeDropdown, trackDropdown} from './dropdown';

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
const drawerPrev = $<HTMLButtonElement>('drawer-prev');
const drawerNext = $<HTMLButtonElement>('drawer-next');
const openAln = $<HTMLButtonElement>('open-aln');
const idModeSelect = $<HTMLSelectElement>('id-mode');
const dlFasta = $<HTMLAnchorElement>('dl-fasta');
const dlAligned = $<HTMLAnchorElement>('dl-aligned');
const dlNwk = $<HTMLAnchorElement>('dl-nwk');
const tabMembers = $('tab-members');
// alignment view
const alnView = $('aln-view');
const alnScroll = $('aln-scroll');
const alnRows = $('aln-rows');
const alnCanvas = $<HTMLCanvasElement>('aln-canvas');
const alnBar = $<HTMLCanvasElement>('aln-bar');
const alnTitle = $('aln-title');
const alnBadges = $('aln-badges');
const alnStatusEl = $('aln-status');
const alnLegend = $('aln-legend');
const alnPrev = $<HTMLButtonElement>('aln-prev');
const alnNext = $<HTMLButtonElement>('aln-next');
const alnClose = $<HTMLButtonElement>('aln-close');
const alnIdMode = $<HTMLSelectElement>('aln-id-mode');
const alnField = $<HTMLSelectElement>('aln-field');
const alnMode = $<HTMLSelectElement>('aln-mode');
const alnRetry = $<HTMLButtonElement>('aln-retry');
const alnDlFasta = $<HTMLAnchorElement>('aln-dl-fasta');
const alnDlAligned = $<HTMLAnchorElement>('aln-dl-aligned');
const alnDlNwk = $<HTMLAnchorElement>('aln-dl-nwk');
const orderSelect = $<HTMLSelectElement>('order');
const modeSelect = $<HTMLSelectElement>('mode');
const themeSelect = $<HTMLSelectElement>('theme');
const treeShapeSelect = $<HTMLSelectElement>('tree-shape');
const sortBtn = $<HTMLButtonElement>('sortbtn');
const statusEl = $('status');
const legendEl = $('legend');
const layersEl = $('layers');
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
// request keys so an unchanged window/order/value is not fetched again
let lastMatrixKey: string | null = null;
let lastTracksKey: string | null = null;
let hover: {row: number; col: number} | null = null;
let highlight: {pos: number; until: number} | null = null;

// fetch / request lifecycle
let gen = 0;
let inflight: AbortController | null = null;
let fetchTimer = 0;
let hashTimer = 0;
let hoverSeq = 0;
let clickSeq = 0;

// detail drawer + alignment view state
let currentOg: string | null = null;
let alignmentView: AlignmentView;
let alnOpen = false;
let alnOg: string | null = null;
let alnStartOg: string | null = null;
let alnFromDrawer = false;
// the inline FAMSA 503 retry also shows a notice; remember it so it can be
// dismissed when the alignment view closes or moves to another OG
let alnRetryToast: {id: number; og: string} | null = null;
let treeField = '';

// members pane: page in genes so huge OGs stay responsive
const MEMBERS_PAGE = 400;
let memberGroups: MemberGroup[] = [];
let membersLimit = MEMBERS_PAGE;
// genome selected by the matrix click (null when the drawer was opened without one)
let selectedGenome: GenomeSelection | null = null;

let activeOgs: {order: string; og0: number; ids: string[]} | null = null;
const ogCache = new Map<string, string[]>();
const ogAnnoCache = new PendingCache<Record<string, unknown>>(64);

// drawer request lifecycle: one generation guards the members pane; the
// alignment view has its own generation because it fetches MSA + tree.
const drawerGen = new Generation();
let drawerSignal: AbortSignal = new AbortController().signal;
const alnGen = new Generation();
let alnSignal: AbortSignal = new AbortController().signal;
const msaCache = new Map<string, MsaResult>();
// guards async prev/next lookups: invalidated on close, route and order changes
const navSeq = new Sequence();
let navChain: Promise<void> = Promise.resolve();
const treeData = new Map<string, GeneTreeArrays>();
const treeState = new Map<string, 'ok' | 'missing'>();
const treeError = new Map<string, string>();
const treeInflight = new Map<string, Promise<GeneTreeArrays>>();

function beginDrawer(): void {
  drawerSignal = drawerGen.begin().signal;
  treeInflight.clear(); // previous requests were aborted; never reuse them
}

function beginAlignment(): number {
  alnSignal = alnGen.begin().signal;
  treeInflight.clear();
  return alnGen.current;
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

/** Alignment-view tree load: uses the alignment generation's signal and the
 *  same caches as the drawer probe. A cached 404 throws immediately. */
async function loadTree(og: string, signal: AbortSignal): Promise<GeneTreeArrays> {
  const cached = treeData.get(og);
  if (cached) return cached;
  if (treeState.get(og) === 'missing') throw new Error(treeError.get(og) ?? 'no gene tree');
  try {
    const d = await getGeneTree(og, signal);
    rememberTree(og, d);
    return d;
  } catch (err) {
    if (err instanceof HttpError && err.status === 404) rememberTree(og, null, (err as Error).message);
    throw err;
  }
}

function dims(): {pw: number; ph: number} {
  return {pw: matrixWrap.clientWidth, ph: matrixWrap.clientHeight};
}

/** Rounded-tile params for the live viewport; read every deck.gl frame. */
function heatmapTileParams(): TileParams {
  if (!lastWin || !lastMat) return FLAT_TILE;
  const binnedX = lastWin.w < lastWin.og1 - lastWin.og0;
  const binnedY = lastWin.h < lastWin.g1 - lastWin.g0;
  const dpr = window.devicePixelRatio || 1;
  return tileParamsFor(2 ** view.zx, 2 ** view.zy, [lastWin.w, lastWin.h], !binnedX && !binnedY, dpr, palette().absentDarken);
}

/** Re-rasterise the current matrix window for the active theme (no refetch). */
function refreshHeatmapRaster(): void {
  if (!lastMat || !lastWin || !heatmap) return;
  heatmap.setData(
    toImageDataFor(lastMat.value, lastMat.data, lastMat.width, lastMat.height, currentTheme()),
    [lastWin.og0, lastWin.g1, lastWin.og1, lastWin.g0],
    heatmapTileParams,
  );
}

/** Redraw every themed surface after a light/dark change. */
function applyTheme(): void {
  setupLegend();
  if (meta) {
    renderTreeLegend(alnLegend, treeField);
    populateLayers(meta);
  }
  refreshHeatmapRaster();
  if (meta) render();
  alignmentView?.draw();
}

const esc = (s: string) =>
  s.replace(/[&<>"]/g, (c) => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c] as string);

const fmtVal = (v: unknown) => (v === null || v === undefined || v === '' ? '—' : String(v));

/** Show a sticky error toast for a failed request (aborts are ignored). */
function toastError(err: unknown, context?: string): void {
  const e = err as Error | undefined;
  if (e?.name === 'AbortError') return;
  if (err instanceof NetworkError) {
    // transport failure from any endpoint collapses to one toast
    toasts().error('server unreachable', {key: 'network', details: {url: err.url, message: err.message}});
    return;
  }
  const http = err instanceof HttpError ? err : null;
  if (!http) {
    // client-side failure (bad payload, parse error, ...): keep it context-specific
    const message = e?.message ?? String(err);
    const text = context ? `${context} failed: ${message}` : message;
    toasts().error(text, {key: `client|${context ?? ''}|${message}`, details: {message}});
    return;
  }
  const message = http.detail ?? http.message;
  const text = context ? `${context} failed: ${message}` : message;
  toasts().error(text, {
    key: errorIdentity(context, http.status, message),
    details: {status: http.status, url: http.url, message},
  });
}

// --- initialization ---------------------------------------------------------

async function boot(): Promise<void> {
  initTheme();
  meta = await getMeta();
  genomes = await getGenomes();
  tree = new Tree(treeCanvas);
  await tree.load();
  tracks = new Tracks(tracksCanvas);
  metadata = new Metadata(metaCanvas);
  alignmentView = new AlignmentView(alnScroll, alnRows, alnCanvas, alnBar, onAlnHover, onAlnBarHover);
  setupDrawer();
  setupAlignment();

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
  alignmentView?.draw();
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
        .map((c, i) => `<span class="swatch" style="background:${categoryColor(i, currentTheme())}" title="${esc(c)}"></span>`)
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
    onError: (err) => toastError(err, 'search'),
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
  // the displayed state was cleared: never dedupe the next request against it
  lastMatrixKey = null;
  lastTracksKey = null;
  hover = null;
  hideTooltip();
  navSeq.invalidate(); // in-flight prev/next must not apply after a reorder
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
  const keys = requestKeys(snapOrder, snapMode, win, activeTracks, activeAnnos);
  const hasMatrix = lastMat !== null && dataOrder === snapOrder;
  const {matrix: needMatrix, tracks: needTracks} = needsFetch(keys, lastMatrixKey, lastTracksKey, hasMatrix);
  const setStatus = () => {
    if (lastMat) statusEl.textContent = `${lastMat.width}×${lastMat.height} cells · ${order} · ${lastMat.value}`;
  };
  if (!needMatrix && !needTracks) {
    // an unchanged window/order/value/layers set: nothing to fetch
    setStatus();
    render();
    return;
  }
  const ac = new AbortController();
  inflight = ac;
  const aborted = () => myGen !== gen || ac.signal.aborted;
  const tasks: Promise<void>[] = [];

  if (needMatrix) {
    statusEl.textContent = 'loading…';
    tasks.push(
      (async () => {
        try {
          const m = await getMatrix({order: snapOrder, value: snapMode, ...win}, ac.signal);
          if (aborted()) return;
          lastMat = m;
          lastWin = win;
          dataOrder = snapOrder;
          lastMatrixKey = keys.matrix;
          heatmap.setData(toImageDataFor(m.value, m.data, m.width, m.height, currentTheme()), [win.og0, win.g1, win.og1, win.g0], heatmapTileParams);
          setStatus();
          render();
        } catch (err) {
          if (ac.signal.aborted || myGen !== gen) return;
          if (err instanceof HttpError && err.status === 404 && snapOrder !== 'tree') {
            toasts().notice(`unknown order '${snapOrder}' — using tree`);
            setOrder('tree');
            return;
          }
          toastError(err, 'matrix');
          statusEl.textContent = `error: ${(err as Error).message}`;
        }
      })(),
    );
  }

  if (needTracks) {
    tasks.push(
      (async () => {
        try {
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
          if (aborted()) return;
          lastTracks = tr;
          lastTracksWorld = tr ? {og0: win.og0, og1: win.og1} : null;
          lastAnnos = annos.filter((a): a is AnnotationTrackResult => a !== null);
          lastAnnoWorld = lastAnnos.length ? {og0: win.og0, og1: win.og1} : null;
          lastTracksKey = keys.tracks;
          render();
        } catch (err) {
          if (aborted()) return;
          toastError(err, 'tracks');
        }
      })(),
    );
  }

  // labels/OG names do not depend on the matrix fetch or the GL frame
  tasks.push(ensureOgNames(snapOrder, win, ac.signal, myGen));
  await Promise.all(tasks);
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
  } catch (err) {
    if (signal.aborted) return;
    toastError(err, 'orthogroup labels');
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
      toastError(err, 'orthogroup details');
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
  const theme = currentTheme();
  tree.draw(tr, treeCanvas.clientWidth, treeCanvas.clientHeight, theme);
  stripsForHit = buildStrips();
  tracks.draw(tr, stripsForHit, tracksCanvas.clientWidth, tracksCanvas.clientHeight, theme);
  if (metaFields.length) metadata.draw(tr, genomes, metaFields, metaCanvas.clientWidth, metaCanvas.clientHeight, theme);
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
  ctx.fillStyle = palette().cladeTint;
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
  ctx.fillStyle = palette().genomeLabel;
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

  const pal = palette();
  // search-hit highlight
  if (highlight && performance.now() < highlight.until) {
    const x0 = tr.xOf(highlight.pos);
    const x1 = tr.xOf(highlight.pos + 1);
    ctx.fillStyle = pal.searchFill;
    ctx.fillRect(x0, 0, Math.max(1, x1 - x0), ph);
    ctx.strokeStyle = pal.searchStroke;
    ctx.strokeRect(x0 + 0.5, 0.5, Math.max(1, x1 - x0), ph - 1);
  }

  if (activeOgs && activeOgs.order === order && tr.sx >= 24) {
    ctx.font = '10px system-ui, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.lineJoin = 'round';
    ctx.miterLimit = 2;
    for (let i = 0; i < activeOgs.ids.length; i++) {
      const x = tr.xOf(activeOgs.og0 + i + 0.5);
      if (x < 0 || x > pw) continue;
      ctx.save();
      ctx.translate(x, 2);
      ctx.rotate(Math.PI / 2);
      // bg-coloured halo keeps the rotated label legible over any cell colour
      ctx.lineWidth = 3;
      ctx.strokeStyle = pal.bg;
      ctx.strokeText(activeOgs.ids[i], 0, 0);
      ctx.fillStyle = pal.text;
      ctx.fillText(activeOgs.ids[i], 0, 0);
      ctx.restore();
    }
    ctx.textAlign = 'left';
  }

  if (hover && lastWin) {
    const [x, y, w, h] = cellRect(lastWin, tr, hover.col, hover.row);
    const tp = heatmapTileParams();
    ctx.strokeStyle = palette().hoverOutline;
    ctx.lineWidth = 1;
    if (tp.mix > 0.5) {
      // match the rounded tile
      const dpr = window.devicePixelRatio || 1;
      const ix = tp.gap[0] / 2 / dpr;
      const iy = tp.gap[1] / 2 / dpr;
      const r = Math.min(tp.radius / dpr, Math.max(0, Math.min(w - 2 * ix, h - 2 * iy) / 2));
      roundRectPath(ctx, x + ix, y + iy, w - 2 * ix, h - 2 * iy, r);
      ctx.stroke();
    } else {
      ctx.strokeRect(x + 0.5, y + 0.5, w, h);
    }
  }
}

function roundRectPath(ctx: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, r: number): void {
  const rr = Math.max(0, Math.min(r, Math.min(w, h) / 2));
  ctx.beginPath();
  ctx.moveTo(x + rr, y);
  ctx.lineTo(x + w - rr, y);
  ctx.arcTo(x + w, y, x + w, y + rr, rr);
  ctx.lineTo(x + w, y + h - rr);
  ctx.arcTo(x + w, y + h, x + w - rr, y + h, rr);
  ctx.lineTo(x + rr, y + h);
  ctx.arcTo(x, y + h, x, y + h - rr, rr);
  ctx.lineTo(x, y + rr);
  ctx.arcTo(x, y, x + rr, y, rr);
  ctx.closePath();
}

function setupLegend(): void {
  const theme = currentTheme();
  if (mode === 'copies') {
    legendEl.innerHTML = copiesLegend(theme)
      .map((b) => `<span style="background:${b.color}" title="${b.label}">${b.label === 'absent' ? '' : b.label}</span>`)
      .join('');
    legendEl.title = 'Copy number: absent, 1, 2, 3, 4, 5+';
    return;
  }
  const stops: string[] = [];
  const n = 16;
  for (let i = 0; i <= n; i++) {
    const [r, g, b] = colorForValue(Math.round((i / n) * 255), theme);
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
      const factor = 2 ** (-e.deltaY * 0.0015);
      view = e.shiftKey
        ? zoomedView(view, meta.n_orthogroups, meta.n_genomes, pw, ph, {x: false, w: tr.rowAt(py), p: py}, factor)
        : zoomedView(view, meta.n_orthogroups, meta.n_genomes, pw, ph, {x: true, w: tr.ogAt(px), p: px}, factor);
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

  treeShapeSelect.onchange = () => {
    tree.setMode(treeShapeSelect.value === 'cladogram' ? 'cladogram' : 'phylogram');
    render();
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
  if (hit.binnedX) {
    toasts().hint('zoom in until the cell covers a single orthogroup to open it');
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
    } catch (err) {
      og = null;
      toastError(err, 'orthogroup id');
    }
  }
  if (my !== clickSeq || !og) return;
  // The clicked row: a single genome, or a range when the row is binned.
  let sel: GenomeSelection | null = null;
  if (hit.binnedY) {
    const [r0, r1] = rowRange(lastWin, hit.row);
    const names: string[] = [];
    for (let r = r0; r < r1 && r < meta.n_genomes; r++) if (genomes[r]) names.push(genomes[r].name);
    if (names.length) sel = rangeSelection(names, [r0, r1]);
  } else {
    const g = genomes[hit.genome];
    if (g) sel = singleSelection(g.name);
  }
  openDrawer(og, sel);
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
  } catch (err) {
    toastError(err, 'search');
  }
}

async function orderByHits(q: string): Promise<void> {
  try {
    const {order: id} = await postOrder({by: 'search', q});
    setOrder(id);
  } catch (err) {
    toastError(err, 'search order');
  }
}

// --- detail drawer ----------------------------------------------------------

function setupDrawer(): void {
  drawerClose.addEventListener('click', closeDrawer);
  drawerPrev.addEventListener('click', () => void stepOg(-1));
  drawerNext.addEventListener('click', () => void stepOg(1));
  openAln.addEventListener('click', () => {
    if (currentOg) openAlignment(currentOg);
  });
  idModeSelect.addEventListener('change', () => {
    updateDownloads();
    updateNwkLink();
  });
  // theme setting: auto / light / dark
  themeSelect.value = getThemeSetting();
  themeSelect.addEventListener('change', () => setThemeSetting(themeSelect.value as ThemeSetting));
  onThemeChange(() => applyTheme());
  // the top bar no longer wraps; pin the Layers dropdown to the viewport
  const layersSummary = layersEl.querySelector('summary')!;
  const layersOpen = () => (layersEl as HTMLDetailsElement).open;
  layersEl.addEventListener('toggle', () => {
    if (layersOpen()) placeDropdown(layersSummary, layerList);
  });
  trackDropdown(layersSummary, layerList, $('topbar'), layersOpen);
  document.addEventListener('keydown', onKeyDown);
  window.addEventListener('popstate', onLocationChange);
  window.addEventListener('hashchange', onLocationChange);
  // one delegated handler for every copy button in the panel
  tabMembers.addEventListener('click', (e) => {
    const b = (e.target as HTMLElement | null)?.closest<HTMLButtonElement>('.copyid');
    if (!b) return;
    const id = b.dataset.id ?? '';
    navigator.clipboard?.writeText(id).then(
      () => toasts().hint(`copied ${id}`),
      () => toasts().error('copy failed'),
    );
  });
}

function onKeyDown(e: KeyboardEvent): void {
  const tag = (document.activeElement as HTMLElement | null)?.tagName;
  const typing = tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA';
  if (e.key === 'Escape') {
    if (typing) return;
    if (alnOpen) {
      e.preventDefault();
      closeAlignment();
    } else if (panel.classList.contains('visible')) {
      e.preventDefault();
      closeDrawer();
    }
    return;
  }
  const prev = e.key === 'ArrowLeft' || e.key === '[';
  const next = e.key === 'ArrowRight' || e.key === ']';
  if ((prev || next) && !typing && (alnOpen || panel.classList.contains('visible'))) {
    e.preventDefault();
    void stepOg(prev ? -1 : 1);
  }
}

function renderTreeLegend(el: HTMLElement, field: string): void {
  const theme = currentTheme();
  const missing = palette().missing;
  const spec = (meta?.genome_metadata ?? []).find((f) => f.name === field);
  let html = '';
  if (!spec || spec.kind === 'categorical') {
    const cats = (spec?.categories ?? []).slice(0, 16);
    html = cats.map((c, i) => `<span class="lg"><i style="background:${categoryColor(i, theme)}"></i>${esc(c)}</span>`).join('');
    if (spec) html += `<span class="lg"><i style="background:${missing}"></i>missing</span>`;
    if (!html) html = `<span class="lg">genome identity</span>`;
  } else {
    const stops = [0, 0.25, 0.5, 0.75, 1].map((f) => `<i style="background:${numericColor(f, theme)}"></i>`).join('');
    html = `<span class="lg">${fmt3(spec.min ?? 0)}${stops}${fmt3(spec.max ?? 1)}</span><span class="lg"><i style="background:${missing}"></i>missing</span>`;
  }
  html += `<span class="lg"><i class="dup"></i>duplication</span>`;
  el.innerHTML = html;
}

function openDrawer(og: string, genome: GenomeSelection | null = null): void {
  selectedGenome = genome;
  const changed = og !== currentOg;
  if (changed) {
    beginDrawer();
    clearDrawerPanes();
    currentOg = og;
    drawerTitle.textContent = og;
    updateDownloads();
    updateNwkLink();
    probeTree(og);
    void renderMembers(og);
  } else {
    // same OG but a new genome selection: refresh just the members view
    renderSelectionBlock();
    renderMemberList();
  }
  panel.classList.add('visible');
  updateHash();
}

function clearDrawerPanes(): void {
  tabMembers.innerHTML = '';
  hideTooltip();
  renderOverlayNow();
}

function updateDownloads(): void {
  if (!currentOg) return;
  const ids = idModeSelect.value === 'raw' ? 'raw' : 'display';
  dlFasta.href = fastaUrl(currentOg, false, ids);
  dlAligned.href = fastaUrl(currentOg, true, ids);
  dlNwk.href = geneTreeNwkUrl(currentOg, ids);
  const aids = alnIdMode.value === 'raw' ? 'raw' : 'display';
  alnDlFasta.href = fastaUrl(currentOg, false, aids);
  alnDlAligned.href = fastaUrl(currentOg, true, aids);
  alnDlNwk.href = geneTreeNwkUrl(currentOg, aids);
}

/** Enable/disable the Newick links from what we already know about the OG. */
function updateNwkLink(): void {
  if (!currentOg) return;
  updateDownloads();
  const enabled = treeState.get(currentOg) !== 'missing';
  for (const el of [dlNwk, alnDlNwk]) {
    el.classList.toggle('disabled', !enabled);
    el.setAttribute('aria-disabled', enabled ? 'false' : 'true');
    if (!enabled) el.removeAttribute('href');
  }
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
    memberGroups = d.members ?? [];
    membersLimit = MEMBERS_PAGE;
    tabMembers.innerHTML = `
      <div id="sel-genome"></div>
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
      <h2>member genomes (${memberGroups.length})</h2>
      <div id="members-list"></div>
      <div id="pattern-section"></div>`;
    renderSelectionBlock();
    renderMemberList();
    void loadPattern(d.orthogroup, my);
  } catch (err) {
    if (!drawerGen.isCurrent(my) || og !== currentOg) return;
    tabMembers.innerHTML = `<h2>error</h2><p>${esc((err as Error).message)}</p>`;
    toastError(err, 'orthogroup');
  }
}

function geneListItem(gene: {id: string; raw: string}): string {
  return (
    `<li><code title="${esc(gene.raw)}">${esc(gene.id)}</code> ` +
    `<button class="copyid" data-id="${esc(gene.raw)}" title="copy raw id">copy</button></li>`
  );
}

/** The selected-genome block shown above the member list. */
function renderSelectionBlock(): void {
  const host = document.getElementById('sel-genome');
  if (!host) return;
  if (!selectedGenome) {
    host.innerHTML = '';
    return;
  }
  const sum = summarizeSelection(selectedGenome, memberGroups);
  if (sum.absent) {
    host.innerHTML = `<div class="sel-genome absent">absent in <b>${esc(sum.label)}</b></div>`;
    return;
  }
  if (sum.binned) {
    let html = `<div class="sel-genome"><b>${esc(sum.label)}</b> — ${sum.withGenes} of ${sum.selected} genomes have this OG`;
    if (sum.tooMany) {
      html += ` <span class="muted">(${sum.genes} genes · zoom in to see genes)</span></div>`;
    } else {
      html += ` <span class="muted">(${sum.genes} genes)</span>`;
      html += `<ul class="genome-list">${sum.groups.map((g) => `<li class="genome-group"><b>${esc(g.genome)}</b><ul>${g.genes.map(geneListItem).join('')}</ul></li>`).join('')}</ul></div>`;
    }
    host.innerHTML = html;
    return;
  }
  const group = sum.groups[0];
  host.innerHTML =
    `<div class="sel-genome"><b>In genome ${esc(sum.label)}</b>: ${sum.genes} ${sum.genes === 1 ? 'gene' : 'genes'}` +
    (group ? `<ul>${group.genes.map(geneListItem).join('')}</ul>` : '') +
    `</div>`;
}

/** Render only the first `membersLimit` genes so a 10k-member OG opens fast;
 *  "show more" pages in the rest, keeping the genome grouping. */
function renderMemberList(): void {
  const host = document.getElementById('members-list');
  if (!host) return;
  const pinned = selectionNames(selectedGenome);
  const {groups, shown, total} = memberSlice(memberGroups, membersLimit, pinned);
  const parts = groups
    .map(
      (g) =>
        `<li class="genome-group${pinned.has(g.genome) ? ' selected' : ''}"><b>${esc(g.genome)}</b>${g.hidden > 0 ? ` <span class="muted">(+${g.hidden})</span>` : ''}<ul>` +
        g.genes.map(geneListItem).join('') +
        `</ul></li>`,
    )
    .join('');
  host.innerHTML =
    `<ul class="genome-list">${parts}</ul>` +
    (shown < total ? `<button id="members-more">show more (${total - shown} genes)</button>` : '');
  const more = host.querySelector<HTMLButtonElement>('#members-more');
  if (more) {
    more.addEventListener('click', () => {
      membersLimit += MEMBERS_PAGE;
      renderMemberList();
    });
  }
}

// --- alignment view ---------------------------------------------------------

const ALN_MSA_CACHE_MAX = 8;

function setupAlignment(): void {
  const fields = meta.genome_metadata ?? [];
  alnField.innerHTML =
    '<option value="">genome identity</option>' +
    fields.map((f) => `<option value="${esc(f.name)}">${esc(f.name)} (${esc(f.kind)})</option>`).join('');
  const firstCat = fields.find((f) => f.kind === 'categorical');
  treeField = firstCat ? firstCat.name : '';
  alnField.value = treeField;
  renderTreeLegend(alnLegend, treeField);
  alnField.addEventListener('change', () => {
    treeField = alnField.value;
    renderTreeLegend(alnLegend, treeField);
    alignmentView.setSpec(specField());
  });
  alnMode.addEventListener('change', () => alignmentView.setMode(alnMode.value as 'phylogram' | 'cladogram'));
  alnIdMode.addEventListener('change', () => {
    alignmentView.setIdMode(alnIdMode.value === 'raw' ? 'raw' : 'display');
    updateDownloads();
    updateNwkLink();
  });
  alnPrev.addEventListener('click', () => void stepOg(-1));
  alnNext.addEventListener('click', () => void stepOg(1));
  alnClose.addEventListener('click', closeAlignment);
  alnRetry.addEventListener('click', () => {
    if (alnOg) void loadAlignment(alnOg);
  });
  const zin = document.getElementById('aln-zoom-in');
  const zout = document.getElementById('aln-zoom-out');
  zin?.addEventListener('click', () => alignmentView.zoomBy(1.5));
  zout?.addEventListener('click', () => alignmentView.zoomBy(1 / 1.5));
}

function specField(): GenomeMetadataInfo | undefined {
  return metaFields.find((f) => f.name === treeField);
}

function alnHash(og: string): string {
  return `#${hashParams()}&og=${encodeURIComponent(og)}&view=aln`;
}

function showAlignment(og: string, fromDrawer: boolean, genome?: GenomeSelection | null): void {
  if (genome !== undefined) selectedGenome = genome;
  alnOpen = true;
  alnOg = og;
  alnStartOg = og;
  alnFromDrawer = fromDrawer;
  alnView.hidden = false;
  alnView.classList.add('open');
  alnTitle.textContent = og;
  currentOg = og;
  updateDownloads();
  updateNwkLink();
  void loadAlignment(og);
}

function openAlignment(og: string): void {
  showAlignment(og, panel.classList.contains('visible'));
  history.pushState(null, '', alnHash(og));
}

function hideAlignmentUI(): void {
  alnOpen = false;
  alnOg = null;
  alnGen.abort();
  navSeq.invalidate(); // a step in flight must not reopen the drawer
  alnView.classList.remove('open');
  alnView.hidden = true;
  hideTooltip();
  if (alnRetryToast) {
    toasts().dismiss(alnRetryToast.id);
    alnRetryToast = null;
  }
}

function closeAlignment(): void {
  if (!alnOpen) return;
  const og = alnOg;
  const stepped = og !== alnStartOg;
  const fromDrawer = alnFromDrawer;
  hideAlignmentUI();
  restoreMatrix(og, stepped, fromDrawer);
}

/** After leaving the alignment view: rebuild the drawer state and centre the
 *  matrix on the last viewed OG when it changed. */
function restoreMatrix(og: string | null, stepped: boolean, fromDrawer: boolean): void {
  if (fromDrawer && og) {
    currentOg = og;
    drawerTitle.textContent = og;
    updateDownloads();
    updateNwkLink();
    void renderMembers(og);
  } else {
    currentOg = null;
  }
  alnStartOg = null;
  alnFromDrawer = false;
  updateHash();
  if (stepped && og) void jumpToOg(og);
}

function onLocationChange(): void {
  navSeq.invalidate(); // a route change discards in-flight prev/next
  const p = new URLSearchParams(location.hash.slice(1));
  const og = p.get('og');
  const view = p.get('view');
  const sel = decodeSelection(p.get('genome'), p.get('rows'), genomes.map((g) => g.name));
  if (view === 'aln' && og) {
    const selChanged = encodeSelection(sel) !== encodeSelection(selectedGenome);
    if (!alnOpen || alnOg !== og || selChanged) showAlignment(og, panel.classList.contains('visible'), sel);
    return;
  }
  if (alnOpen) {
    const cur = alnOg;
    // honour an explicit selection in the target hash (e.g. a hand-edited
    // genome/rows); a genome-less Back keeps the current stepped selection
    if (sel) selectedGenome = sel;
    // On Back the hash still names the OG the view was opened from; when the
    // hash was edited to another OG, honour that one instead.
    const target = og && og !== alnStartOg ? og : cur;
    hideAlignmentUI();
    restoreMatrix(target, target !== alnStartOg, alnFromDrawer);
    return;
  }
  if (og) {
    // selection-only changes (same OG, different genome/rows) must apply too
    if (og !== currentOg || encodeSelection(sel) !== encodeSelection(selectedGenome)) openDrawer(og, sel);
  } else if (panel.classList.contains('visible')) closeDrawer();
}

async function positionOfOg(og: string): Promise<number | null> {
  const pos = activeOgs && activeOgs.order === order ? locateOg(activeOgs.ids, activeOgs.og0, og) : null;
  if (pos !== null) return pos;
  try {
    const res = await apiSearch(og, order, 5);
    const hit = res.hits.find((h) => h.orthogroup === og) ?? res.hits[0];
    return hit?.position ?? null;
  } catch {
    return null;
  }
}

/** Enqueue a prev/next step so lookups never overlap. */
function stepOg(delta: number): void {
  const token = navSeq.current(); // captured at enqueue: a later close/route/order change cancels it
  navChain = navChain.then(() => stepOgNow(delta, token)).catch(() => {
    /* ignore */
  });
}

async function stepOgNow(delta: number, token: number): Promise<void> {
  if (!navSeq.isCurrent(token)) return;
  const og = alnOpen ? alnOg : currentOg;
  if (!og) return;
  const pos = await positionOfOg(og);
  if (!navSeq.isCurrent(token)) return; // closed / route or order changed
  if (pos === null) return;
  const next = stepPosition(pos, delta, meta.n_orthogroups);
  if (next === null) return;
  try {
    const ids = await getOrthogroups(order, next, next + 1);
    if (!navSeq.isCurrent(token)) return;
    const id = ids[0];
    if (!id) return;
    if (alnOpen) {
      alnOg = id;
      currentOg = id;
      alnTitle.textContent = id;
      updateDownloads();
      updateNwkLink();
      history.replaceState(null, '', alnHash(id));
      void loadAlignment(id);
    } else {
      openDrawer(id, selectedGenome);
    }
  } catch {
    /* ignore */
  }
}

function badgeHtml(m: MsaResult): string {
  return `<span class="badge src">${esc(m.source)}</span>` + (m.trimmed === true ? ' <span class="badge">trimmed</span>' : '');
}

function rememberMsa(og: string, m: MsaResult): void {
  msaCache.delete(og);
  msaCache.set(og, m);
  while (msaCache.size > ALN_MSA_CACHE_MAX) {
    const oldest = msaCache.keys().next().value as string | undefined;
    if (oldest === undefined) break;
    msaCache.delete(oldest);
  }
}

async function loadAlignment(og: string): Promise<void> {
  const my = beginAlignment();
  // a pending retry notice belongs to the OG it was raised for
  if (alnRetryToast && alnRetryToast.og !== og) {
    toasts().dismiss(alnRetryToast.id);
    alnRetryToast = null;
  }
  alnStatusEl.textContent = 'loading…';
  alnBadges.innerHTML = '';
  alnRetry.hidden = true;
  alnTitle.textContent = og;
  // drop the previous OG's rows/hover immediately; the header shows "loading…"
  if (alignmentView.contentOg !== og) alignmentView.clear();
  const spec = specField();
  alignmentView.setSpec(spec);
  alignmentView.setIdMode(alnIdMode.value === 'raw' ? 'raw' : 'display');
  alignmentView.setMode(alnMode.value as 'phylogram' | 'cladogram');

  let tree: GeneTreeArrays | null = null;
  let treeMessage: string | null = null;
  if (treeState.get(og) === 'missing') {
    treeMessage = treeError.get(og) ?? 'no gene tree (fewer than 4 genes)';
  } else {
    try {
      tree = treeData.get(og) ?? (await loadTree(og, alnSignal));
    } catch {
      tree = null;
      treeMessage = treeError.get(og) ?? null;
    }
  }
  if (!alnGen.isCurrent(my) || alnOg !== og) return;
  updateNwkLink();

  const build = (msa: MsaResult | null, msg: string | null): AlnContent => {
    const {rows, hasTree} = orderAlignmentRows(msa, tree);
    const cons = msa ? computeConservation(msa.rows, msa.length) : {conservation: null, gap: null};
    return {
      og,
      rows,
      hasTree,
      treeMessage,
      msaLength: msa?.length ?? 0,
      conservation: cons.conservation,
      gap: cons.gap,
      msaMessage: msg,
      genomes,
      spec,
      idMode: alnIdMode.value === 'raw' ? 'raw' : 'display',
      selectedGenomes: [...selectionNames(selectedGenome)],
    };
  };
  const apply = (msa: MsaResult | null, msg: string | null): void => {
    alignmentView.setTree(tree);
    alignmentView.setData(build(msa, msg));
  };

  const cached = msaCache.get(og) ?? null;
  apply(cached, cached ? null : 'loading alignment…');
  if (cached) {
    alnStatusEl.textContent = `${cached.rows.length} rows × ${cached.length} cols`;
    alnBadges.innerHTML = badgeHtml(cached);
    updateNwkLink();
    return;
  }

  try {
    const msa = await getMsa(og, alnSignal);
    if (!alnGen.isCurrent(my) || alnOg !== og) return;
    rememberMsa(og, msa);
    apply(msa, null);
    alnStatusEl.textContent = `${msa.rows.length} rows × ${msa.length} cols`;
    alnBadges.innerHTML = badgeHtml(msa);
    updateNwkLink();
  } catch (err) {
    if (!alnGen.isCurrent(my) || alnOg !== og) return;
    if ((err as Error).name === 'AbortError') return;
    apply(null, (err as Error).message);
    alnStatusEl.textContent = (err as Error).message;
    if (err instanceof HttpError && err.status === 503) {
      alnRetry.hidden = false;
      const id = toasts().notice('alignment limit busy, try again', {
        action: {
          label: 'Retry',
          onClick: () => {
            if (alnOpen && alnOg === og) void loadAlignment(og);
          },
        },
      });
      alnRetryToast = {id, og};
    } else {
      toastError(err, 'alignment');
    }
  }
}

function pct(x: number): string {
  return `${Math.round(x * 100)}%`;
}

function onAlnHover(h: AlnHover | null, clientX: number, clientY: number): void {
  if (!h) {
    hideTooltip();
    return;
  }
  const loc = h.col >= 0
    ? `<br>col ${h.col + 1} · ${esc(h.residue || '-')}${h.ungapped > 0 ? ` · ungapped ${h.ungapped}` : ''}<br>conservation ${pct(h.conservation)} · gap ${pct(h.gap)}`
    : '';
  tooltip.innerHTML = `<b>${esc(h.gene || '—')}</b>${h.genome ? ` · <span class="tip-dim">${esc(h.genome)}</span>` : ''}${loc}`;
  tooltip.style.display = 'block';
  placeTooltip(clientX, clientY);
}

function onAlnBarHover(h: AlnBarHover | null, clientX: number, clientY: number): void {
  if (!h) {
    hideTooltip();
    return;
  }
  tooltip.innerHTML = `col ${h.col + 1} · conservation ${pct(h.conservation)} · gap ${pct(h.gap)}`;
  tooltip.style.display = 'block';
  placeTooltip(clientX, clientY);
}

function closeDrawer(): void {
  beginDrawer();
  navSeq.invalidate(); // discard any in-flight prev/next
  panel.classList.remove('visible');
  currentOg = null;
  selectedGenome = null;
  clearDrawerPanes();
  updateHash();
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
  } catch (err) {
    host.innerHTML = '';
    toastError(err, 'pattern');
  }
}

async function orderByPattern(og: string): Promise<void> {
  try {
    const {order: id} = await postOrder({by: 'pattern', og});
    setOrder(id);
  } catch (err) {
    toastError(err, 'pattern order');
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
    toastError(err, 'contrast sort');
  } finally {
    sortBtn.disabled = false;
  }
}

// --- URL hash (order, viewport, mode, search, drawer, alignment view) --------

function scheduleHash(): void {
  window.clearTimeout(hashTimer);
  hashTimer = window.setTimeout(updateHash, 200);
}

/** Shared matrix part of the hash (order, viewport, mode, search, selected genome). */
function hashParams(): string {
  const q = searchBox?.query ?? '';
  let h =
    `order=${encodeURIComponent(order)}&cx=${view.cx.toFixed(2)}&cy=${view.cy.toFixed(2)}` +
    `&zx=${view.zx.toFixed(3)}&zy=${view.zy.toFixed(3)}&mode=${mode}`;
  if (q) h += `&q=${encodeURIComponent(q)}`;
  h += encodeSelection(selectedGenome);
  return h;
}

function updateHash(): void {
  let h = hashParams();
  if (alnOpen && alnOg) h += `&og=${encodeURIComponent(alnOg)}&view=aln`;
  else if (currentOg && panel.classList.contains('visible')) h += `&og=${encodeURIComponent(currentOg)}`;
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
  const view_ = p.get('view');
  const sel = decodeSelection(p.get('genome'), p.get('rows'), genomes.map((g) => g.name));
  if (og && view_ === 'aln') showAlignment(og, panel.classList.contains('visible'), sel);
  else if (og) openDrawer(og, sel);
}

boot().catch((err) => {
  statusEl.textContent = `fatal: ${(err as Error).message}`;
  toastError(err, 'startup');
});
