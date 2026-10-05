# OrthoPan web viewer (MVP)

Vite + TypeScript + deck.gl frontend. It talks only to the API in `docs/api.md`.

## Layout

`src/` modules (plain TypeScript, no framework; deck.gl draws the heatmap, the rest is Canvas2D/DOM):

- `main.ts` — wiring and UI state: boot, fetch lifecycle, render loop, interaction, URL hash, keyboard.
- `theme.ts` — theme resolution (auto/light/dark), the light/dark palette, and the live theme store.
- `transform.ts` — the shared world↔pixel transform and viewport math (zoom/clamp, window and cell
  maths, DPR canvas setup). Every matrix-aligned panel uses it.
- `api.ts` — typed client for every endpoint, plus `HttpError`/`NetworkError`.
- `matrix.ts` — deck.gl heatmap: colour maths, the `TiledBitmapLayer` rounded-cell shader and the
  `Heatmap` with its flat fallback.
- `tree.ts` — genome dendrogram on Canvas2D, clade selection and hit testing.
- `tracks.ts` — track strips above the matrix (numeric bars, categorical colours).
- `annotations.ts` — the categorical annotation palette.
- `metadata.ts` — genome metadata columns between tree and heatmap.
- `search.ts` — the search box and hit viewport.
- `selection.ts` — the selected-genome model and hash encoding.
- `members.ts` — member-list pagination with pinned genomes.
- `pending.ts`, `generation.ts` — the shared-promise cache and the stale-response guards.
- `requestkey.ts` — request dedupe keys for the matrix window and its track layers.
- `dropdown.ts` — position-fixed dropdowns anchored to a top-bar control.
- `toast.ts` — the toast queue and container (hint/notice/error).
- `alignment.ts` — the alignment view (row order, virtual band, renderer).
- `msa.ts` — MSA column maths, residue colours, conservation.
- `genetree.ts` — gene-tree layout, `treeNodeX`, leaf colours.
- `decode.ts` — little-endian binary decoders and copy-number byte maths.
- `mock/server.ts` — synthetic 50 genomes × 2,000 OGs API (Bun) serving this contract.

## Run

```sh
cd web
bun install

# development: mock API + vite dev server (proxy /api -> 127.0.0.1:8000)
bun run mock        # terminal 1
bun run dev         # terminal 2, http://127.0.0.1:5173

# or against the real backend (built with scripts/build_datasets.sh)
uv run orthopan serve ../scratch/datasets/demo --port 8000

# production build (type check + bundle) into ../src/orthopan/static
bun run build
bun run mock        # serves the built UI at http://127.0.0.1:8000/
```

## Interaction

- wheel zooms the OG axis, shift+wheel zooms the genome axis, drag pans.
- order dropdown: `tree`, `frequency`, `track:<name>`, `track:<name>:desc`.
- hover a cell for a tooltip; click for orthogroup details in the side panel.
- click a tree node to select its clade, then "sort by contrast".
- order and viewport are stored in the URL hash.
