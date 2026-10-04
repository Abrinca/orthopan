# OrthoPan viewer API (draft v0)

The FastAPI backend serves one imported pangenome directory (see `orthopan import`) and the built frontend. Everything lives under `/api`. JSON is used for metadata. Matrices and tracks are raw little-endian binary (`application/octet-stream`).

## Coordinates

- **Genome axis**: rows in tree order (`genomes.row_order`), `0 .. n_genomes`.
- **Orthogroup axis**: positions `0 .. n_orthogroups` in the selected **order**.
- **Windows** are half-open ranges, `[og0, og1)` × `[g0, g1)`.
- **Resolution**: the client asks for `w` × `h` output cells, with `w <= og1 - og0` and `h <= g1 - g0`.
  - Output column `i` covers OG positions `[og0 + floor(i*(og1-og0)/w), og0 + floor((i+1)*(og1-og0)/w))`.
  - Rows are binned the same way.
  - When `w == og1 - og0` and `h == g1 - g0`, the response is at full detail.

## Orders

- `order` names an orthogroup permutation:
  - Built in: `tree` (default) and `frequency`.
  - `track:<name>` and `track:<name>:desc` sort by a track column. Missing values go last.
  - Ids returned by `POST /api/orders`.
- The server caches permutations, using an LRU for the dynamic ones.

## Endpoints

### `GET /api/meta`
```json
{"n_genomes": 24, "n_orthogroups": 5657, "tree_source": "input",
 "orders": ["tree", "frequency"],
 "tracks": [{"name": "copy_number_mean", "min": 1.0, "max": 12.3, "missing": 0}]}
```

### `GET /api/genomes`
A JSON list in row order. Each entry has `{"row", "genome_idx", "name", "n_orthogroups", "n_genes"}`.

### `GET /api/tree`
The tree as pre-order node arrays, so the client needs no Newick parser:
```json
{"parent": [-1, 0, 1], "length": [0, 0.2, 0.05], "label": ["", "N1", "FAM1"], "row": [-1, -1, 0]}
```
`row` holds the genome row for leaves and -1 for internal nodes.

### `GET /api/matrix?order=&og0=&og1=&g0=&g1=&w=&h=`
- Returns a `Uint8Array` of `h*w` values, row-major (genome rows, then OG columns).
- Each value is the fraction of present cells in the bin, scaled to 0–255. At full detail it is 0 or 255.
- Response headers: `X-Width`, `X-Height`.
- The total requested output is capped at 4,000,000 cells (`w * h`); larger requests are rejected with 400.
- Copy-number mode (`value=copies`) is planned for later.

### `GET /api/tracks?order=&og0=&og1=&w=&names=a,b`
- Returns a `Float32Array` of `len(names)*w` values, track-major.
- Each bin holds the mean of its non-missing values. NaN means no values.
- Response header: `X-Tracks` (the comma-separated names, in order).

### `GET /api/orthogroups?order=&og0=&og1=`
The orthogroup IDs at those positions, used for labels at high zoom. At most 5,000 per request.

### `GET /api/orthogroup/{og}`
Details for one orthogroup, addressed by ID:
```json
{"orthogroup": "OG0000123", "og_idx": 123, "n_genomes": 20, "n_genes": 21, "mrca_label": "N5", "purity": 0.83,
 "tracks": {"copy_number_mean": 1.05},
 "members": [{"genome": "FAM1079-i1-1.1", "genes": ["gnl|extdb|FAM1079-i1-1.1_000123"]}]}
```

### `POST /api/orders`
```json
{"by": "contrast", "genomes": ["FAM1079-i1-1.1", "FAM3257-i1-1.1"]}
```
- Returns `{"order": "<id>"}`, which can be used as the `order` parameter.
- A contrast score is the fraction present inside the selection minus the fraction present outside it. Orthogroups are sorted by this score, descending.

## Performance targets

Reference dataset: 2,193 × 50,515. Long-term target: 10k genomes.

- `/api/matrix` for a full-screen window (about 2000 × 1000 output cells): under 100 ms.
- `/api/tracks`: under 50 ms.

## Notes

- The server is intended to run on localhost; it has no authentication or global request
  limiter.
- The default `tree` order is packed at startup. Dynamic orders (`POST /api/orders`) are
  packed lazily on first use, which costs a few hundred milliseconds on the reference
  dataset; later requests with the same order id are fast.
