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
- Caching: the `tree` and `frequency` **permutations** are precomputed at import; `frequency`'s
  packed matrix is built in the background at startup. `track:<name>` permutations are **recomputed
  on every request** (no cache). The packed matrices for `tree`, `frequency` and every `track:*`
  share one 8-entry LRU, so heavy track sorting can evict `tree` (it is rebuilt transparently on the
  next request). Dynamic orders (`POST /api/orders`) have their own 8-entry LRU; an evicted id
  returns 404.

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
- Copy-number mode is `value=copies` (see the v1 additions below).

### `GET /api/tracks?order=&og0=&og1=&w=&names=a,b`
- Returns a `Float32Array` of `len(names)*w` values, track-major.
- Each bin holds the mean of its non-missing values. NaN means no values.
- Response header: `X-Tracks` (the comma-separated names, in order).
- Without `names`, every track is returned. An unknown `names` entry returns 404.

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
- A request body larger than 64 KiB returns 413 (checked against `Content-Length` and while streaming).

## Performance targets

Reference dataset: 2,193 × 50,515. Long-term target: 10k genomes.

- `/api/matrix` for a full-screen window (about 2000 × 1000 output cells): under 100 ms.
- `/api/tracks`: under 50 ms.

## Notes

- The server is intended to run on localhost; it has no authentication or global request
  limiter.
- The default `tree` order is packed at startup, and `frequency` is packed in the background.
  Dynamic orders (`POST /api/orders`) are packed lazily on first use, which costs a few hundred
  milliseconds on the reference dataset; later requests with the same order id are fast.

## v1 additions: annotations, genome metadata, search, copy number

These additions are backwards compatible. Existing endpoints keep their behaviour unless a change is noted here.

### `GET /api/meta` (extended)
New fields:
```json
{"annotations": [{"name": "description", "kind": "text"},
                 {"name": "cog_category", "kind": "categorical", "categories": ["J", "K", "L"]}],
 "genome_metadata": [{"name": "clade", "kind": "categorical", "categories": ["A", "B"]},
                     {"name": "gc", "kind": "numeric", "min": 0.35, "max": 0.41}],
 "detail": {"msa": true, "gene_trees": true, "sequences": true}}
```
- A text column is free text: tooltip, panel and search only.
- A categorical column is a string column with at most 64 distinct values. `categories` lists them by frequency, most common first.

### `GET /api/genomes` (extended)
Each entry gains `"metadata": {"clade": "A", "gc": 0.38}`. A missing value is `null`.

### `GET /api/matrix` (extended): `value=presence|copies`
- `presence` (default) is unchanged.
- `copies` returns a `Uint8` per bin: 0 if the bin has no present cell, otherwise `min(255, round(10 * mean copy number over present cells))`. So 10 means 1.0 copies and 25 means 2.5.

### `GET /api/annotation-track?order=&og0=&og1=&w=&name=`
- For a categorical annotation, returns a `Uint16Array` of `w` values.
- Each value is `1 + index` into `categories` for the most common category in the bin. 0 means missing.
- Ties go to the lower index.
- Response headers: `X-Width` and `X-Categories` (the annotation name).
- A text annotation, or a name that is not an annotation, returns 400.

### `GET /api/search?q=&order=&limit=50`
- Case-insensitive substring search over orthogroup IDs and all text and categorical annotations.
- Returns `[{"orthogroup", "og_idx", "position", "field", "text"}]`, where `position` is the OG's position in `order`.
- Results are sorted by position. `limit` is capped at 500. The response also gives the total hit count as `{"hits": [...], "total": N}`; `total` counts **(orthogroup, field) pairs**, so one OG matched in two fields counts twice.

### `POST /api/orders` (extended)
`{"by": "search", "q": "transporter"}` returns an order with the hits first, keeping their relative order from `tree`, followed by the rest in `tree` order. The same id rules and LRU apply as for contrast orders.

### `GET /api/orthogroup/{og}` (extended)
- Adds `"annotations": {"description": "..."}`, `"pattern_id": 123` and `"n_same_pattern": 4`.
- Each member gene is `{"id": "FAM1_000123", "raw": "gnl|extdb|FAM1_000123"}`. `id` is the display ID with the prefix stripped. `raw` is the ID as in the input, used to match OrthoFinder files.
- Note: `members[].genes` changes from a list of strings to a list of these objects.

### `GET /api/orthogroup/{og}/pattern?limit=200`
- Returns the OGs whose presence/absence pattern is identical to this one: `{"pattern_id", "total", "orthogroups": [{"orthogroup", "og_idx", "annotations"}]}`.
- The list **includes the OG itself**, in `og_idx` order. `limit` is capped at 5,000. The frontend filters the OG out when it displays the list.
- `POST /api/orders {"by": "pattern", "og": "OG0000123"}` returns an order with the pattern group first.

## v2 additions: single-orthogroup detail view

Data sources:
- OrthoFinder results recorded at import (`--orthofinder`): `MultipleSequenceAlignments/`, `Orthogroup_Sequences/` and `Resolved_Gene_Trees/Resolved_Gene_Trees.txt`.
- Protein FASTA given at import with `--fasta` (files or directories). These are indexed by gene ID with byte offsets, not copied.
- If an orthogroup has no MSA but its sequences are available, the server aligns it on demand. It runs the same command OrthoFinder uses, `famsa INPUT OUTPUT -t 1`, with the binary set by `orthopan serve --famsa PATH` (default: `famsa` on `PATH`). Results are cached in `<data_dir>/cache/msa/<OG>.fa`. At most 2 alignments run at once, and only for orthogroups with at most 2,000 sequences.

### `GET /api/meta` (extended)
`"detail": {"msa": bool, "gene_trees": bool, "sequences": bool, "famsa": bool}`. `famsa` is true when the binary is usable.

### `GET /api/orthogroup/{og}/msa`
```json
{"orthogroup": "OG0000123", "source": "orthofinder", "length": 412, "trimmed": false,
 "rows": [{"genome": "FAM1079-i1-1.1", "row": 5, "gene": "FAM1079-i1-1.1_000123",
           "raw": "gnl|extdb|FAM1079-i1-1.1_000123", "seq": "MKT-LV..."}]}
```
- `source` is `"orthofinder"` or `"famsa"`.
- Rows are sorted by genome row (tree order), then by gene ID.
- `row` is the genome row; it is -1 if the gene's genome can't be resolved.
- `trimmed` uses the same rule as the `trimmed` track. It is null if unknown.
- Errors: 404 with `{"detail": "<reason>"}` when the orthogroup has no MSA and none can be made (no sequences, no famsa, or too many sequences). 503 when the alignment concurrency limit is busy. A FAMSA failure or timeout (120 s) returns 500.
- `aligned=1` downloads use the same MSA path, so they can also return 404, 500 or 503.

### `GET /api/orthogroup/{og}/gene-tree`
Pre-order node arrays, like `/api/tree`:
```json
{"orthogroup": "OG0000123", "parent": [-1, 0], "length": [0, 0.1], "label": ["n0", "FAM1_x"],
 "gene": [null, "FAM1079-i1-1.1_000123"], "genome_row": [-1, 5],
 "duplication": [true, false], "support": [0.8, null]}
```
- Leaves are matched to genes and genomes using the OG's membership.
- `duplication` and `support` come from `Gene_Duplication_Events/Duplications.tsv`, matched by the gene tree's node label.
- 404 with a reason when there is no tree. OrthoFinder only builds trees for orthogroups with at least 4 genes.

### Downloads
- `GET /api/orthogroup/{og}/fasta?aligned=0|1&ids=display|raw`: `text/plain` FASTA with `Content-Disposition: attachment; filename="<OG>[.aligned].fa"`. `aligned=1` uses the same source as `/msa`.
- `GET /api/orthogroup/{og}/gene-tree.nwk?ids=display|raw`: Newick, as an attachment.
