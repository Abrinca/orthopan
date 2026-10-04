// OrthoPan MVP viewer: wires the tree, heatmap, tracks, labels and interaction.
import {
  HttpError,
  getGenomes,
  getMatrix,
  getMeta,
  getOrthogroup,
  getOrthogroups,
  getTracks,
  postOrder,
  type Genome,
  type MatrixResult,
  type Meta,
  type TracksResult,
} from './api';
import {Heatmap, colorForValue, toImageData} from './matrix';
import {
  bufferAt,
  cellRect,
  clampView,
  colRange,
  computeWindow,
  fitView,
  makeTransform,
  rowRange,
  type View,
  type Window,
} from './transform';
import {Tracks, fmt3} from './tracks';
import {Tree} from './tree';

const $ = <T extends HTMLElement>(id: string): T => document.getElementById(id) as T;

const matrixWrap = $('matrix-wrap');
const treeCanvas = $<HTMLCanvasElement>('tree');
const tracksCanvas = $<HTMLCanvasElement>('tracks-canvas');
const labelsCol = $('labels-col');
const labelsCanvas = $<HTMLCanvasElement>('labels');
const overlay = $<HTMLCanvasElement>('overlay');
const tooltip = $('tooltip');
const panel = $('panel');
const orderSelect = $<HTMLSelectElement>('order');
const sortBtn = $<HTMLButtonElement>('sortbtn');
const statusEl = $('status');
const noticeEl = $('notice');
const legendEl = $('legend');

let meta: Meta;
let genomes: Genome[] = [];
let tree: Tree;
let heatmap: Heatmap;
let tracks: Tracks;
let order = 'tree';
let view: View = {cx: 0, cy: 0, zx: 0, zy: 0};
let fitted = true;
let activeTracks: string[] = [];
let lastWin: Window | null = null;
let lastMat: MatrixResult | null = null;
// the order the currently displayed matrix/tracks belong to
let dataOrder: string | null = null;
let lastTracks: TracksResult | null = null;
let lastTracksWorld: {og0: number; og1: number} | null = null;
let hover: {row: number; col: number} | null = null;

// fetch / request lifecycle
let gen = 0;
let inflight: AbortController | null = null;
let fetchTimer = 0;
let hashTimer = 0;
let detailSeq = 0;

// OG-id cache: the active window plus a small bounded LRU.
let activeOgs: {order: string; og0: number; ids: string[]} | null = null;
const ogCache = new Map<string, string[]>();

function dims(): {pw: number; ph: number} {
  return {pw: matrixWrap.clientWidth, ph: matrixWrap.clientHeight};
}

const esc = (s: string) =>
  s.replace(/[&<>"]/g, (c) => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c] as string);

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

  activeTracks = pickTracks(meta);
  setupLegend();

  const {pw, ph} = dims();
  view = fitView(meta.n_orthogroups, meta.n_genomes, pw, ph);
  populateOrders(meta);
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

  window.addEventListener('resize', onResize);
  new ResizeObserver(onResize).observe(matrixWrap);

  bindInteractions();
  updateHash();
  scheduleFetch();
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

function setOrder(o: string): void {
  order = o;
  if (!hasOption(o)) addOption(o);
  orderSelect.value = o;
  activeOgs = null;
  ogCache.clear();
  // The displayed data now belongs to the old order: make it non-interactive
  // immediately and cancel any pending OG lookup / panel request.
  dataOrder = null;
  lastTracks = null;
  lastTracksWorld = null;
  hover = null;
  detailSeq++;
  hideTooltip();
  scheduleRender();
  scheduleHash();
  scheduleFetch();
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
  // Invalidate immediately so any in-flight response is discarded.
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
  const {pw, ph} = dims();
  if (pw < 2 || ph < 2) return;
  const win = computeWindow(view, meta.n_orthogroups, meta.n_genomes, pw, ph);
  const ac = new AbortController();
  inflight = ac;
  statusEl.textContent = 'loading…';
  try {
    const m = await getMatrix({order: snapOrder, ...win}, ac.signal);
    if (myGen !== gen || ac.signal.aborted) return;
    lastMat = m;
    lastWin = win;
    dataOrder = snapOrder;
    heatmap.setData(toImageData(m.data, m.width, m.height), [win.og0, win.g1, win.og1, win.g0]);
    statusEl.textContent = `${m.width}×${m.height} cells · ${snapOrder}`;
    render();

    if (activeTracks.length) {
      const tr = await getTracks(
        {order: snapOrder, og0: win.og0, og1: win.og1, w: win.w, names: activeTracks},
        ac.signal,
      );
      if (myGen !== gen || ac.signal.aborted) return;
      lastTracks = tr;
      lastTracksWorld = {og0: win.og0, og1: win.og1};
      render();
    }
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

function ogAt(pos: number): string | null {
  if (activeOgs && activeOgs.order === order && pos >= activeOgs.og0 && pos < activeOgs.og0 + activeOgs.ids.length) {
    return activeOgs.ids[pos - activeOgs.og0] ?? null;
  }
  return null;
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
  if (lastTracks && lastTracksWorld) {
    tracks.draw(tr, lastTracks, lastTracksWorld, meta.tracks, activeTracks, tracksCanvas.clientWidth, tracksCanvas.clientHeight);
  } else {
    const ctx = tracksCanvas.getContext('2d');
    ctx?.clearRect(0, 0, tracksCanvas.width, tracksCanvas.height);
  }
  drawLabels(tr, labelsCanvas.clientWidth, labelsCanvas.clientHeight);
  drawOverlay(tr, pw, ph);
  updateSortButton();
}

function syncLabelColumn(): void {
  // Collapse the genome-label column while rows are too dense for text; the
  // matrix width then changes and the shared transform/deck resize with it.
  const collapsed = 2 ** view.zy < 11;
  labelsCol.classList.toggle('collapsed', collapsed);
  panel.style.right = `${labelsCol.offsetWidth}px`;
}

function fitCanvas(c: HTMLCanvasElement, pw: number, ph: number): CanvasRenderingContext2D | null {
  const dpr = window.devicePixelRatio || 1;
  if (c.width !== Math.round(pw * dpr) || c.height !== Math.round(ph * dpr)) {
    c.width = Math.round(pw * dpr);
    c.height = Math.round(ph * dpr);
  }
  const ctx = c.getContext('2d');
  if (!ctx) return null;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, pw, ph);
  return ctx;
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
  const ctx = fitCanvas(labelsCanvas, pw, ph);
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
  const ctx = fitCanvas(overlay, pw, ph);
  if (!ctx) return;

  // clade tint across the matrix rows (requirement: extend across the matrix)
  drawCladeTint(ctx, tr, pw);

  // orthogroup labels at high zoom, rotated downward so nothing is clipped
  if (activeOgs && activeOgs.order === order && tr.sx >= 24) {
    ctx.font = '10px system-ui, sans-serif';
    ctx.fillStyle = '#475467';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    for (let i = 0; i < activeOgs.ids.length; i++) {
      const x = tr.xOf(activeOgs.og0 + i + 0.5);
      if (x < 0 || x > pw) continue;
      ctx.save();
      ctx.translate(x + 2, 2);
      ctx.rotate(Math.PI / 2);
      ctx.fillText(activeOgs.ids[i], 0, 0);
      ctx.restore();
    }
    ctx.textAlign = 'left';
  }

  // hovered buffer cell outline
  if (hover && lastWin) {
    const [x, y, w, h] = cellRect(lastWin, tr, hover.col, hover.row);
    ctx.strokeStyle = 'rgba(17,24,39,0.7)';
    ctx.lineWidth = 1;
    ctx.strokeRect(x + 0.5, y + 0.5, w, h);
  }
}

function setupLegend(): void {
  const stops: string[] = [];
  const n = 16;
  for (let i = 0; i <= n; i++) {
    const [r, g, b] = colorForValue(Math.round((i / n) * 255));
    stops.push(`rgb(${r},${g},${b}) ${((i / n) * 100).toFixed(0)}%`);
  }
  legendEl.style.background = `linear-gradient(90deg, ${stops.join(',')})`;
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
    if (Math.hypot(e.clientX - down.x, e.clientY - down.y) > 4) return; // drag, not click
    void onClick(e);
  });

  treeCanvas.addEventListener('click', (e) => {
    const rect = treeCanvas.getBoundingClientRect();
    const idx = tree.hitTest(e.clientX - rect.left, e.clientY - rect.top);
    if (idx !== null) {
      tree.selected = tree.selected === idx ? null : idx;
      render();
    }
  });

  sortBtn.onclick = () => void sortByContrast();
}

function ogRangeText(a: string, b: string): string {
  return a === b ? a : `${a} – ${b}`;
}

function onHover(e: MouseEvent): void {
  if (!lastMat || !lastWin || !meta || dataOrder !== order) {
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
    // pointer is outside the fetched window (stale data): not interactive
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
    detail = `${ogRangeText(ogA, ogB)} · present ${Math.round((value / 255) * 100)}%`;
  } else {
    nameText = genomes[hit.genome]?.name ?? `row ${hit.genome}`;
    const og = ogAt(hit.og) ?? `#${hit.og}`;
    detail = `${og} · ${value > 0 ? 'present' : 'absent'}`;
  }
  tooltip.innerHTML = `<b>${esc(nameText)}</b><br>${esc(detail)}`;
  tooltip.style.display = 'block';
  tooltip.style.left = `${px + 14}px`;
  tooltip.style.top = `${py + 14}px`;
  drawOverlay(tr, pw, ph);
}

async function onClick(e: MouseEvent): Promise<void> {
  if (!lastMat || !lastWin || !meta || dataOrder !== order) return;
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
    const members = (d.members ?? [])
      .map((m) => `<li><code>${esc(m.genome)}</code> — ${m.genes.map(esc).join(', ')}</li>`)
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
        ${trackRows}
      </table>
      <h2>members (${d.members?.length ?? 0})</h2>
      <ul>${members}</ul>`;
    panel.querySelector('.close')!.addEventListener('click', closePanel);
  } catch (err) {
    if (my !== detailSeq) return;
    panel.innerHTML = `<button class="close">×</button><h2>error</h2><p>${esc((err as Error).message)}</p>`;
    panel.querySelector('.close')!.addEventListener('click', closePanel);
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

// --- URL hash (order + viewport) -------------------------------------------

function scheduleHash(): void {
  window.clearTimeout(hashTimer);
  hashTimer = window.setTimeout(updateHash, 200);
}

function updateHash(): void {
  const h = `order=${encodeURIComponent(order)}&cx=${view.cx.toFixed(2)}&cy=${view.cy.toFixed(2)}&zx=${view.zx.toFixed(3)}&zy=${view.zy.toFixed(3)}`;
  history.replaceState(null, '', `#${h}`);
}

function applyHash(): void {
  const p = new URLSearchParams(location.hash.slice(1));
  const o = p.get('order');
  if (o) {
    order = o;
    if (!hasOption(o)) addOption(o);
  }
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
}

boot().catch((err) => {
  statusEl.textContent = `fatal: ${(err as Error).message}`;
});
