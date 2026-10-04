# OrthoPan web viewer (MVP)

Vite + TypeScript + deck.gl frontend. It talks only to the API in `docs/api.md`.

## Layout

- `src/api.ts` — typed client for every endpoint.
- `src/matrix.ts` — deck.gl `BitmapLayer` heatmap, viewport math and colormap.
- `src/tree.ts` — rectangular dendrogram on Canvas2D, aligned with matrix rows.
- `src/tracks.ts` — top track strips aligned with matrix columns.
- `src/main.ts` — wiring: layout, fetches, interaction, URL hash.
- `mock/server.ts` — synthetic 50 genomes x 2000 OGs API (Bun) serving this contract.

## Run

```sh
cd web
bun install

# development: mock API + vite dev server (proxy /api -> 127.0.0.1:8000)
bun run mock        # terminal 1
bun run dev         # terminal 2, http://127.0.0.1:5173

# or against the real backend
uv run orthopan serve ../scratch/import-smoke/demo-import --port 8000

# production build (type check + bundle) into web/dist
bun run build
bun run mock        # serves web/dist at http://127.0.0.1:8000/
```

## Interaction

- wheel zooms the OG axis, shift+wheel zooms the genome axis, drag pans.
- order dropdown: `tree`, `frequency`, `track:<name>`, `track:<name>:desc`.
- hover a cell for a tooltip; click for orthogroup details in the side panel.
- click a tree node to select its clade, then "sort by contrast".
- order and viewport are stored in the URL hash.
