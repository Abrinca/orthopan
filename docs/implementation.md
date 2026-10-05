# OrthoPan implementation guide

This document explains how OrthoPan works internally. It is written for developers who maintain or extend it. Other documents cover related ground:

- [`README.md`](../README.md): installation and usage.
- [`docs/development.md`](development.md): the environment, test datasets and day-to-day workflow.
- [`docs/design-decisions.md`](design-decisions.md): why things are built this way.
- [`docs/api.md`](api.md): the HTTP contract (parameters, response shapes, status codes).
- [`docs/anvio-internals.md`](anvio-internals.md): the homogeneity maths and the anvi'o analysis behind them.
- [`docs/orthofinder-outputs.md`](orthofinder-outputs.md) and [`docs/ortholog_file_formats_documentation.md`](ortholog_file_formats_documentation.md): the input formats.

Code is referenced by file path and function or class name. Line numbers are deliberately not used because they drift. Statements about behaviour were checked against the code at commit `61d3b2b`. Unless a measurement says otherwise, numbers come from the code's constants.

## Contents

1. [Overview](#1-overview)
2. [Input parsing (`io.py`)](#2-input-parsing-iopy)
3. [The import pipeline (`importer.py`)](#3-the-import-pipeline-importerpy)
4. [The dataset directory format](#4-the-dataset-directory-format)
5. [Server internals (`server.py`)](#5-server-internals-serverpy)
6. [Track preprocessors (`tracks/`)](#6-track-preprocessors-tracks)
7. [Frontend architecture (`web/src`)](#7-frontend-architecture-websrc)
8. [Testing](#8-testing)
9. [Performance and scale](#9-performance-and-scale)
10. [Extending OrthoPan](#10-extending-orthopan)
11. [Known limitations and gotchas](#11-known-limitations-and-gotchas)

---

## 1. Overview

### 1.1 Components and data flow

OrthoPan has two Python commands and a browser UI. The entry point is `orthopan.main` in `src/orthopan/__init__.py`, which registers the `import` and `serve` subcommands from `importer.add_parser` and `server.add_parser`.

```
 inputs                         orthopan import                 dataset directory
 ─────────────────────────      ──────────────────────────      ─────────────────────────
 orthologs TSV (A or B) ──┐     io.read_matrix (streaming)      manifest.json
 --tree NEWICK ───────────┤     tree: parse / UPGMA             genomes.parquet
 --orthofinder RESULTS ───┤──▶  presence bit-packing      ──▶   orthogroups.parquet
 --track TSV ...  ────────┤     membership batches              membership.parquet
 --annotation FILE ... ───┤     orderings, pattern_id           presence.npy
 --genome-metadata TSV ───┤     track preprocessors             multicopy.npz
 --fasta PATH ... ────────┘     FASTA byte-offset index         annotations.parquet  (opt.)
                                                                genome_metadata.parquet (opt.)
                                                                sequences.parquet    (opt.)
                                                                tree.nwk, cache/msa/ (lazy)
                                                                         │
                                                                         ▼
 browser (Vite + TS + deck.gl)       HTTP /api/*           orthopan serve (FastAPI/uvicorn)
 ─────────────────────────────  ◀──────────────────────▶  ──────────────────────────────
 main.ts wiring, transform.ts       JSON + raw LE binary    Pangenome: data in memory,
 heatmap (BitmapLayer)                                      packed per-order matrices,
 tree / tracks / metadata canvases                          popcount aggregation,
 drawer, alignment view                                     on-demand MSA / gene tree,
                                                            static files at /
```

OrthoFinder result files (MSAs, `Orthogroup_Sequences/`, gene trees, `Duplications.tsv`) and `--fasta` files are **not** copied. The manifest stores their absolute paths, and the server reads them on demand.

### 1.2 Repository map

| Path | What it does |
| --- | --- |
| `src/orthopan/__init__.py` | `main()`: argparse entry point for `orthopan import` and `orthopan serve`. |
| `src/orthopan/io.py` | Reads orthologs TSVs: Format A/B detection, padding, prefix stripping, gzip detection and the closeable row iterators. |
| `src/orthopan/importer.py` | `orthopan import`: the Newick tokenizer and parser, `_Tree` (MRCA queries), track, annotation and metadata readers, the streaming import in `run()`, and the dataset writer. |
| `src/orthopan/tree.py` | `orthopan-tree` and `upgma_newick()`: UPGMA on Jaccard distances computed from packed presence, plus an iterative Newick writer and `_quote()`. |
| `src/orthopan/sequences.py` | The protein FASTA byte-offset index (`index_files`, `resolve`, `read_record`), used by the importer and the server. |
| `src/orthopan/server.py` | `orthopan serve`: the `Pangenome` class (all data and algorithms) and `create_app()` (routes, validation, static mount). |
| `src/orthopan/tracks/__init__.py` | `write_tracks()`: the shared track TSV writer. |
| `src/orthopan/tracks/basic.py` | Paralog-fraction and length-CV tracks for any tool (copy number is computed at import). Also hosts `fasta_files()` and `_stem()`, which `sequences.py` reuses. |
| `src/orthopan/tracks/orthofinder.py` | Duplication counts and gene-tree length tracks from an OrthoFinder results directory. |
| `src/orthopan/tracks/msa.py` | Geometric and functional homogeneity, gap fraction and the trimmed flag from MSAs. Also provides the MSA header matching (`_expected_headers`, `_match_headers`) and `read_fasta`, which the server reuses. |
| `src/orthopan/static/` | The built UI, written by `bun run build`. It is gitignored but included in wheels. |
| `web/index.html` | Page layout, CSS and every DOM element that `main.ts` looks up by id. |
| `web/src/*.ts` | The frontend modules (see [§7.1](#71-module-map)). |
| `web/mock/server.ts` | A Bun mock of the whole API: synthetic data, 50 genomes × 2,000 OGs. |
| `web/test/*.test.ts` | Bun unit tests for the pure frontend logic. |
| `scripts/simulate_pangenome.py` | A seeded simulator for large, tree-structured pangenomes. It writes a gzipped Format B TSV, a tree, annotations and metadata. |
| `scripts/orthofinder/*.sh` | Run OrthoFinder 3.1.5 on `resources/orthofinder-demo` (podman image, or a core + `--assign` run from the `scratch/.venv` environment). |
| `tests/` | pytest suite (see [§8](#8-testing)). `tests/data/` is downloaded on first use by `conftest.py`. |

### 1.3 Identifiers used throughout

Four index spaces appear everywhere. Mixing them up is the most common source of bugs.

| Name | Meaning | Defined by |
| --- | --- | --- |
| `genome_idx` | Genome column order in the input matrix (Format A: order of first appearance). | `genomes.parquet.genome_idx` |
| row, `row_order` | Genome position in the tree: the pre-order rank of its leaf (0 = top row). | `genomes.parquet.row_order`, `_Tree.leaf_pos` |
| `og_idx` | Orthogroup row order in the input matrix. | `orthogroups.parquet.og_idx` |
| position | Orthogroup position (column) under the current *order*. | `perm[position] = og_idx` on the server |

The API speaks in rows and positions ([`docs/api.md`](api.md) §Coordinates). Storage uses `genome_idx` and `og_idx`. The only exception is `multicopy.npz`, which stores tree **rows** (see [§4.6](#46-multicopynpz)).

---

## 2. Input parsing (`io.py`)

`read_matrix(path, strip_prefix=DEFAULT_GENE_PREFIX)` returns `(labels, rows)`:

- `labels` is the list of genome names.
- `rows` is an iterator of `(orthogroup_id, cells)`. It always yields exactly `len(labels)` cells, and each cell is a comma-separated gene list exactly as written (the empty string means absent).

The default prefix is `DEFAULT_GENE_PREFIX = r"^gnl\|[^|]+\|"`.

### 2.1 Gzip detection

`open_text()` reads the first two bytes and opens the file with `gzip.open` if they are the gzip magic number `1f 8b`. The suffix does not matter. `tree._matrix` and the importer's track, annotation and metadata readers use the same function, so every TSV input may be gzipped. FASTA input is the exception (see [§3.5](#35-the---fasta-byte-offset-index)).

### 2.2 Format detection

Blank lines are skipped. Detection looks at the first non-blank line, split on tabs:

| First line | Result |
| --- | --- |
| First cell is `Orthogroup` or `HOG` | **Format B with header.** The labels are the remaining header cells. |
| Exactly two cells, no header | **Format A** (two columns, `OG<TAB>g1, g2, ...`). |
| More than two cells, no header | **Format B without header.** The labels are `genome_1 .. genome_{n}`, and the first line is kept as data (`chain([first], numbered)`). |
| Empty file, or a header with no genome columns | `ValueError("input has no genome columns")` |

A two-column file *with* an `Orthogroup`/`HOG` header is a one-genome Format B matrix, not Format A (`test_io.py::test_format_a_two_column_with_header_is_matrix`).

### 2.3 Format B: padding and validation

`_MatrixRows._generate` splits each line on tabs. Rows shorter than the header are padded with `""`. Writers often drop trailing empty cells, so this is common. A row with **more** cells than genomes raises `ValueError("line N: ...")` with the 1-based line number in the file.

### 2.4 Format A: two passes and genome derivation

Format A has no genome columns, so the genome is derived from each gene ID by `genome_from_gene()`: strip the prefix, then take everything before the last `_` (`rsplit("_", 1)[0]`). A gene without `_` is its own genome label (`test_format_a_no_underscore_gene`).

1. `_format_a_labels()` makes a first full pass over the file and collects labels in order of first appearance.
2. The file is reopened. `_FormatARows` buckets each row's genes into per-genome lists and yields them as `", "`-joined cells, so downstream code sees Format B shaped rows. A gene whose label was not seen in pass 1 raises an error, which cannot happen for a file that is unchanged between passes.

The cells keep the **raw** gene IDs with the prefix still attached. Prefix stripping for display happens later, in the importer.

### 2.5 Prefix stripping and raw IDs

`compile_prefix(pattern)` returns `None` for `None` or `""`, which disables stripping. `strip_gene_prefix()` applies `pattern.sub("", gene)`. `split_genes()` splits on commas and trims whitespace.

The importer stores both forms (`membership.genes` and `membership.genes_raw`). The display ID is for people. The raw ID is needed to match OrthoFinder MSA headers, gene-tree leaves and FASTA headers, which all contain the prefix. See [design-decisions](design-decisions.md) §Annotations.

### 2.6 The closeable iterator

`_OwnedRows` is an iterator that **owns the open file**:

- Iterating to exhaustion closes it, because the generator's `finally` calls `fh.close()`.
- `close()`, `with rows:` and `__del__` close it deterministically when the caller stops early, for example after a validation error.

The importer wraps its whole streaming phase in `with rows:` for this reason. Tests: `test_io.py::test_close_before_iterating_closes_file`, `test_with_closes_file_when_body_raises`, `test_format_a_close_semantics`.

---

## 3. The import pipeline (`importer.py`)

### 3.1 Step order in `run(args)`

```
resolve inputs (--orthofinder?)                      ─┐
tree text: --tree | OrthoFinder species tree | UPGMA  │ before the matrix is opened
_Tree(_parse_newick(text))                           ─┘
read_matrix() ── with rows: ───────────────────────────────────────────────┐
   validate labels vs tree leaves                                          │
   --fasta: index files (in memory)                                        │ streaming
   read track sources (+ OrthoFinder preprocessors, in-process)            │ phase
   for each OG row: presence bits, stats, MRCA, membership, multicopy,     │
                    FASTA resolution, track lookup                         │
   flush membership batches                                                │
   write sequences.parquet; parse annotations and genome metadata ─────────┘
warnings for unknown track / annotation / metadata IDs
orderings (_ranks), tree.nwk, presence.npy, pattern_id, multicopy.npz,
orthogroups.parquet, genomes.parquet, annotations.parquet, genome_metadata.parquet,
membership.parquet (merge parts), manifest.json
```

Outside the streaming loop the per-OG and per-genome lists are O(n_og)/O(n_genomes), but three things are not and stay in memory: the `--fasta` byte-offset index and `seq_rows` (O(FASTA records) / O(resolved genes)), and the `_pattern_ids` dense unpack (O(n_genomes × n_og) bytes). The per-cell data is otherwise touched only inside the loop and in 10k-row membership batches.

### 3.2 Input resolution and `--orthofinder`

- With `--orthofinder RESULTS_DIR`, the orthologs path defaults to `RESULTS_DIR/Orthogroups/Orthogroups.tsv`. An explicit positional `orthologs` overrides it, but the results directory is still used for the tree, the tracks and the detail paths.
- Without either input, the command exits with `"import: provide ORTHOLOGS or --orthofinder RESULTS_DIR"`.
- The output directory is created if it is missing. Files from a previous import in the same directory are overwritten, and a stale optional output (`sequences.parquet`, `annotations.parquet`, `genome_metadata.parquet`) is removed when the new import has none. `cache/` is **not** cleared (see [§11](#11-known-limitations-and-gotchas)).

What `--orthofinder` turns on:

| Helper | Effect |
| --- | --- |
| `_orthofinder_tree_text()` | Uses `Species_Tree/SpeciesTree_rooted_node_labels.txt` if present, otherwise `SpeciesTree_rooted.txt`, unless `--tree` is given. `tree_source` becomes `"orthofinder"`. |
| `_orthofinder_track_sources()` | Runs `tracks.orthofinder.build_tracks(results_dir)` and, if `MultipleSequenceAlignments/` exists, `tracks.msa.build_tracks(msa_dir, seqs=Orthogroup_Sequences or None, orthogroups=Orthogroups.tsv)`. Both run **in-process** (no subprocess or temporary TSV). Their `{og: {col: value}}` dicts become track sources through `_dict_source()`. The MSA step uses a process pool (see [§6.4](#64-msapy-parallelism)). `tracks.basic` is not run. |
| `_orthofinder_detail()` | Records `msa`, `gene_trees` and `sequences` as **resolved absolute paths**, or `null` if the file or directory does not exist. These go into `manifest.detail`. |

### 3.3 Tree handling

#### Tokenizer: `_tokens(text)`

`_tokens` yields `("punct", c)` for `( ) , ; :` and `("text", value)` otherwise:

- Whitespace is skipped.
- `[...]` comments are skipped, including nested ones. Unterminated comments and a stray `]` raise errors.
- `'quoted labels'` become one text token, even when they contain punctuation. `''` inside quotes is an escaped `'`.
- Unquoted text runs until one of `(),:; \t\r\n[]'`.

#### Parser state machine: `_parse_newick(text)`

The parser is non-recursive, so deep caterpillar trees work (`test_deep_caterpillar_tree_is_iterative`). It keeps a `stack` of open internal nodes, the last completed node `last`, and a `state`:

```
             label                      ':'                 number
 subtree ─────────────▶ after_label ──────────▶ length ───────────▶ after_length
   │  ▲                     ▲                                         │
 '(' │  │ ','                 │ label (only after ')')                 │
   ▼  │                     │                                         │
 (push node, stay)      after_close ◀────────── ')' (pop) ◀───────────┘
                                         (')' and ',' are also legal
                                          from after_label and after_close)
 ';' in any state other than subtree/length ends parsing.
```

- In `subtree`, a text token creates a leaf and `(` pushes a new internal node. Any other token is an error.
- In `length`, the token must be text that parses as a float.
- In the "node complete" states (`after_label`, `after_length`, `after_close`):
  - `;` stops parsing.
  - `:` starts a branch length; a second length raises "repeated branch length".
  - `,` returns to `subtree`; outside a group it raises "comma outside a group".
  - `)` pops the stack.
  - Text is accepted only in `after_close`, where it becomes an internal node label.

After the loop the parser rejects an open stack, a dangling `:`, an empty input and any token after `;`. A missing `;` is accepted (`test_valid_newick_without_trailing_semicolon`). `add()` raises "more than one root" when a second top-level node appears.

The server reuses `_parse_newick` for the genome tree and for gene trees.

#### `_Tree`: pre-order arrays, leaf intervals and MRCA

`_Tree(root)` flattens the node graph iteratively in **pre-order**. Children keep their file order. It then builds:

- `label`, `children` and `parent` arrays indexed by pre-order index.
- `leaf_order`: leaf pre-order indices in left-to-right order. `leaf_pos[i]` is the leaf rank, which **is the genome row** (`genomes.row_order`).
- `lo[i]` and `hi[i]`: the leaf-rank interval covered by node `i`. Leaves get their own rank, and internal nodes take `lo` of the first child and `hi` of the last, filled in reverse pre-order. Every subtree covers a contiguous range of rows because rows are leaf order.
- `up`: binary-lifting ancestor tables. `up[0]` is the parent (the root points to itself) and `up[k][i] = up[k-1][up[k-1][i]]`, for `max(1, n.bit_length())` levels.

`mrca(min_leaf, max_leaf)` finds the deepest node whose interval contains `[min_leaf, max_leaf]`:

1. Start at the leaf `leaf_order[min_leaf]`. If its `hi` already reaches `max_leaf`, which only happens for a single-genome OG, return it.
2. Otherwise climb from the highest power of two down, jumping to an ancestor whenever `hi[ancestor] < max_leaf`. This lands on the highest ancestor that still does **not** cover `max_leaf`.
3. Return its parent.

This is correct because the ancestors of a leaf have `lo <= min_leaf`, and their `hi` is non-decreasing going up. The MRCA of all present genomes therefore depends only on the minimum and maximum present leaf rank, which the streaming loop tracks per OG. Each query costs O(log n). `test_mrca_matches_bruteforce_on_random_trees` and `test_polytomy_mrca` cover it.

#### Label validation

Inside `with rows:` the importer checks three things:

1. The matrix has no duplicate genome labels.
2. The tree has no duplicate leaf labels.
3. The set of leaf labels equals the set of matrix labels. On failure the error lists `tree only` and `matrix only` names.

It then maps each `genome_idx` to its leaf rank (`genome_leafpos`).

#### UPGMA fallback: `tree.upgma_newick(path, strip_prefix)`

When neither `--tree` nor an OrthoFinder tree is available, the importer calls `upgma_newick` before it opens the matrix. `tree_source` becomes `"upgma"`. This reads the input separately from the main import.

1. `_matrix()` builds a genomes × OGs boolean matrix.
   - Format A goes through `read_matrix`.
   - Format B is parsed with a byte-level fast path: tab positions come from `np.frombuffer(line) == 9`, and a cell counts as present when it is non-empty. Short rows are zero-padded, and over-wide rows raise an error.
2. The matrix is bit-packed (`np.packbits(axis=1, bitorder="big")`) and the boolean copy is dropped.
3. `_jaccard_packed()` accumulates `inter = X @ X.T` in float32. It unpacks one chunk of 8,192 OG bits at a time, and multiplies in row blocks of 2,048 into a reusable `(block, n)` buffer. No float copy of the full matrix is ever made.
4. The distance is `1 − inter/union` with `union = |a| + |b| − inter`. Two empty rows get distance 0, matching scipy. The result is condensed with `squareform`.
5. `scipy.cluster.hierarchy.linkage(method="average")` and `to_tree()` build the tree, and `_newick()` writes it iteratively. Branch lengths are `(parent height − child height) / 2`, formatted with `.6g`. `_quote()` quotes labels that contain whitespace or `(),:;'[]`.

A single genome produces `"<label>;"`.

The `orthopan-tree` command (`tree.main`) runs the same function. It has no `--strip-prefix` option, so it always uses the default prefix.

### 3.4 Streaming the matrix

For each `(og, cells)` row in `run()`, `og_idx` is the running row counter. The loop performs the following steps.

**Presence bit-packing.** `presence` is a list of one `bytearray` per genome, grown in 4,096-byte steps for every genome at once (`new_len = ((byte >> 12) + 1) << 12`). OG `og_idx` lives in byte `og_idx >> 3` under the mask `1 << (7 - (og_idx & 7))`. That is MSB-first, the same layout as `np.packbits(..., bitorder="big")`. At the end, each buffer is truncated to `nbytes = ceil(n_og / 8)`, the buffers are concatenated, and the result is saved as a `(n_genomes, nbytes)` uint8 array. Rows are in `genome_idx` order.

**Per-cell work.** A cell is skipped if it is empty or contains no gene IDs. Otherwise:

- `raw_genes` are the split IDs as written, and `genes` are the prefix-stripped display IDs.
- `copies = len(raw_genes)`.
- The presence bit is set, and the genome's `n_orthogroups` and `n_genes` counters are updated.
- The minimum and maximum leaf rank of the OG are tracked for the MRCA.
- A membership row `(og_idx, genome_idx, copies, genes, genes_raw)` is appended.
- If `copies > 1`, a multicopy entry `(tree row, og_idx, copies − 1)` is appended.
- If `--fasta` is set, each gene is resolved against the FASTA index (see [§3.5](#35-the---fasta-byte-offset-index)).

**Per-OG statistics.** For each OG the loop records `n_genomes` (present genomes) and `n_genes` (total copies). If any genome is present it also records:

- `copy_number_mean = n_genes / n_genomes` and `copy_number_max`;
- `mrca_node = tree.mrca(min_leaf, max_leaf)` (a pre-order index), and `mrca_label` from it;
- `purity = n_genomes / (hi − lo + 1)`, the fraction of the MRCA clade's genomes that have the OG.

An OG with no genes gets null for all of these.

**Track values.** Track values are looked up per OG from the parsed track sources (see [§3.6](#36-tracks-annotations-and-genome-metadata)).

**Membership batches.** When the buffered membership reaches `_MEMBERSHIP_BATCH = 10_000` rows, `flush_membership()` writes a part file named `_membership.<k>.parquet` in the output directory. The gene lists are built through a flat `(row, gene)` frame and a `group_by("row", maintain_order=True).agg(...)` join, so the `List(Utf8)` columns are built by polars and not as Python lists of lists. At the end:

- One part is renamed to `membership.parquet`.
- Several parts are merged with `pl.scan_parquet(parts).sink_parquet(..., row_group_size=_MEMBERSHIP_ROW_GROUP)`, which is 10,000 rows, and then deleted.
- If there were no members, an empty frame with the full schema is written.

The rows come out sorted by `og_idx` and then `genome_idx`. The server's row-group pruning depends on this (see [§5.9](#59-orthogroup-detail-and-row-group-pruning)).

### 3.5 The `--fasta` byte-offset index

`--fasta` takes one or more paths (`nargs="+"`). The steps are:

1. `seqmod.fasta_files()` expands directories into their sorted `*.fa`, `*.faa` and `*.fasta` files (each optionally `.gz`).
2. `warn_gzip()` drops `.gz` files with a warning, because gzip cannot be read at random offsets.
3. The remaining paths are resolved to absolute paths and stored in `manifest.fasta`, where their list index is the file index.

`sequences.index_files(files, prefix)` reads each file in binary, line by line. For every record it stores `(file_idx, offset, byte_length, seq_len)`:

- `offset` is the byte position of the `>` line.
- `byte_length` covers the header and sequence lines up to the next record.
- `seq_len` counts the non-whitespace bytes of the sequence.

Records are keyed by the first word of the header and, when the prefix pattern changes that word, also by the stripped word. The whole index is a Python dict held in memory during the import.

During streaming, each `(raw, genome)` pair is resolved once with `sequences.resolve()`:

1. Look up the raw ID, then the display ID. The first key that has any records decides the result.
2. A single record wins.
3. With several records, prefer the one whose file stem (`_stem`, without the FASTA extension) maps to this genome through `stem_to_genome`.
4. Otherwise, if all candidates have the same length **and** identical sequences (read back from disk), take the first.
5. Otherwise the gene is ambiguous and unresolved.

`sequences.parquet` gets one row per resolved `(gene=raw, genome)`. The import warns with the count of unresolved membership genes. The server re-reads a record with `read_record()`: it seeks, reads `length` bytes, and joins the non-header lines.

### 3.6 Tracks, annotations and genome metadata

**Track files.** `_read_track_sources()` always treats the first line as the header and reads columns `header[1:]`. Cells pass through `_float_or_none`: an empty cell, NaN or ±inf becomes `None`. Short rows are padded with `None`. Non-numeric text raises a `ValueError` from `float()`.

`_merge_sources(sources, reserved=_ORTHOGROUP_COLUMNS)` concatenates the column lists. A name that is already used, either by another track or by a built-in orthogroup column, raises `duplicate track column name`. The reserved set is every column the importer writes to `orthogroups.parquet`, including `copy_number_mean` and `copy_number_max`.

During streaming, each OG gets `values.get(og)[j]` or `None`. After streaming, the importer warns about track IDs that are not in the matrix. They are ignored.

**Annotations.** `_read_annotation()` uses the first non-blank line to decide the layout:

- If its first cell, lowercased, is `orthogroup`, `hog` or `og`, the file has a header and named columns.
- Otherwise there is no header, and the single column is named by `_annotation_name()`: the file stem without `.gz`, or `description` when the stem contains `descriptions`.

Typing is done in `_annotation_specs()` and uses **only OG IDs that are in the matrix**, so stray IDs cannot change the typing or the category order:

| Layout | Kind |
| --- | --- |
| Headerless | always `text` |
| Header, ≤ 64 distinct non-empty values | `categorical`, with `categories` sorted by descending count, then by value |
| Header, > 64 distinct values | `text` |

Names must be non-empty, unique across all annotation files, and must not collide with `_ORTHOGROUP_COLUMNS`. An annotation may share a name with a track, because tracks and annotations are stored in separate files.

**Genome metadata.** `_read_genome_metadata()` always expects a header with the genome name in the first column. Column names must be non-empty, unique, and not in `_GENOME_COLUMNS` (`genome_idx`, `name`, `row_order`, `n_orthogroups`, `n_genes`). Cells that parse as a non-finite float become missing.

A column is `numeric`, with `min` and `max`, if every non-empty value for a known genome parses as a float. Otherwise it is `categorical`, with categories by frequency. Unlike annotations, metadata has **no 64-value cap**: every non-numeric column is categorical. Unknown genome names produce a warning.

### 3.7 Orderings and rank columns

Both built-in orders are stored as **rank** columns (`order_tree[og_idx] = position`) computed by `_ranks(n, key)`:

| Column | Sort key per OG (ascending) |
| --- | --- |
| `order_tree` | `(0, mrca_node, −purity, −n_genomes, og_idx)`; OGs with no genes get `(1, …, og_idx)` and go last. |
| `order_frequency` | `(−n_genomes, og_idx)` |

`mrca_node` is a pre-order index, so OGs are grouped by the DFS position of their MRCA. Core OGs (MRCA = root, index 0) come first, then each clade's OGs follow the tree from top to bottom. Within one MRCA, purer OGs come first, then more frequent ones, and `og_idx` breaks the remaining ties (`test_same_clade_rank_ties_break_by_index`). The server turns each rank column into a permutation with `np.argsort(kind="stable")`. For the reasoning see [design-decisions](design-decisions.md) §Orthogroup ordering.

### 3.8 `pattern_id`

`_pattern_ids(packed, n_og)` groups OGs whose presence columns are identical:

1. Unpack `presence` to a dense genomes × OGs uint8 matrix.
2. Repack it along the **genome** axis (`packbits(axis=0)`), so each OG column becomes a `ceil(n_genomes/8)`-byte key.
3. Walk the OGs in `og_idx` order and give each new key the next id.

Pattern ids are therefore numbered in order of first appearance. `n_same_pattern` is the size of the OG's group, including the OG itself. Grouping is exact: there is no hashing, so collisions are impossible. The dense unpack is the import's largest temporary allocation (see [§9](#9-performance-and-scale)).

### 3.9 Writing outputs and the manifest

`tree.nwk` gets the tree text exactly as read, with a trailing newline added if missing. For an OrthoFinder or `--tree` input this is the original file, labels and all. The Parquet schemas are spelled out explicitly; [§4](#4-the-dataset-directory-format) lists them. Annotation columns are forced to `Utf8`, because an all-null column would otherwise be inferred as `Null` and break the server's string search (`test_search_nullable_annotations`).

`FORMAT_VERSION = 2`. Version 2 introduced `genes_raw`, `multicopy.npz` and `pattern_id`/`n_same_pattern` together with the v1 viewer features. The server refuses any directory with `format_version < 2`, and a missing field counts as version 1 (`test_outdated_format_version_fails`). The v2 detail-view additions (`sequences.parquet`, `manifest.detail`, `manifest.fasta`) are optional and did not bump the version. For the rule on when to bump it, see [§10.6](#106-when-to-bump-format_version).

---

## 4. The dataset directory format

These schemas were checked against `scratch/datasets/demo` (an `--orthofinder` import of the 24-genome demo with `--genome-metadata`) and `scratch/datasets/hog` (a Format A import with `--annotation` and `--fasta`). Both are rebuilt by `scripts/build_datasets.sh`.

| File | Always present | Written by |
| --- | --- | --- |
| `manifest.json` | yes | importer |
| `genomes.parquet` | yes | importer |
| `orthogroups.parquet` | yes | importer |
| `membership.parquet` | yes | importer |
| `presence.npy` | yes | importer |
| `multicopy.npz` | yes (may be empty) | importer |
| `tree.nwk` | yes | importer |
| `annotations.parquet` | only with annotations | importer |
| `genome_metadata.parquet` | only with `--genome-metadata` | importer |
| `sequences.parquet` | only with `--fasta` and ≥ 1 resolved gene | importer |
| `cache/msa/<OG>.fa` | after the first on-demand FAMSA alignment | server |
| `_membership.<k>.parquet` | only while an import is running | importer |

### 4.1 `manifest.json`

| Field | Type | Meaning |
| --- | --- | --- |
| `format_version` | int | `2`. |
| `orthologs` | str | Orthologs path **as given** (may be relative). Informational only. |
| `orthofinder` | str \| null | Results directory as given. |
| `tree` | str \| null | `--tree` path as given. |
| `tree_source` | `"input"` \| `"orthofinder"` \| `"upgma"` | Where the tree came from. Shown in `/api/meta`. |
| `tracks` | [str] | `--track` paths as given. |
| `track_columns` | [str] | Extra track columns in `orthogroups.parquet`, in order: user tracks first, then the OrthoFinder and MSA preprocessor tracks. The server adds `copy_number_mean` and `copy_number_max` in front. |
| `fasta` | [str] | Absolute FASTA paths. The list index is `sequences.parquet.file`. |
| `annotations` | [spec] | `{"name", "kind": "text"\|"categorical", "categories"?}`. |
| `genome_metadata` | [spec] | `{"name", "kind": "categorical", "categories"}` or `{"name", "kind": "numeric", "min", "max"}`. |
| `detail` | object | `{}`, or `{"msa", "gene_trees", "sequences"}`, each an absolute path or `null` (OrthoFinder imports only). |
| `n_genomes`, `n_orthogroups` | int | Matrix size. |
| `n_membership` | int | Number of present cells (membership rows). |
| `n_genes` | int | Total genes (sum of copies). |

The server reads `format_version`, `track_columns`, `annotations`, `genome_metadata`, `detail`, `fasta` and `tree_source`. The other fields are provenance.

### 4.2 `genomes.parquet`

There is one row per genome in `genome_idx` order.

| Column | dtype | Meaning |
| --- | --- | --- |
| `genome_idx` | UInt32 | 0..n−1, matrix column order. |
| `name` | Utf8 | Genome label, equal to the tree leaf label. |
| `row_order` | UInt32 | Tree row (leaf rank in pre-order). |
| `n_orthogroups` | UInt32 | OGs present in this genome. |
| `n_genes` | UInt32 | Genes of this genome across all OGs. |

### 4.3 `orthogroups.parquet`

There is one row per OG in `og_idx` order.

| Column | dtype | Meaning |
| --- | --- | --- |
| `og_idx` | UInt32 | 0..n−1, input row order. |
| `orthogroup` | Utf8 | OG ID as in the input. |
| `n_genomes` | UInt32 | Present genomes. |
| `n_genes` | UInt32 | Total genes. |
| `copy_number_mean` | Float64 | `n_genes / n_genomes`; null when absent everywhere. |
| `copy_number_max` | Float64 | Largest copy number; null when absent everywhere. |
| `mrca_node` | UInt32 | Pre-order index (in `tree.nwk`) of the MRCA of present genomes; null when absent everywhere. |
| `mrca_label` | Utf8 | That node's label, which may be null or empty for unlabelled nodes. |
| `purity` | Float64 | `n_genomes / clade size`. |
| `order_tree` | UInt32 | Position in the `tree` order (rank). |
| `order_frequency` | UInt32 | Position in the `frequency` order. |
| `pattern_id` | UInt32 | Identical-presence group id. |
| `n_same_pattern` | UInt32 | Size of that group. |
| *track columns* | Float64 | One per `manifest.track_columns`; null means missing. |

The MSA preprocessor's `trimmed` column is stored as Float64 1.0/0.0/null like any other track. The server's `/msa` endpoint converts it to a boolean.

### 4.4 `membership.parquet`

This is the sparse matrix: one row per present `(og_idx, genome_idx)` cell, sorted by `og_idx`, then `genome_idx`. When the import produced more than one batch, the file is written in row groups of 10,000 rows.

| Column | dtype | Meaning |
| --- | --- | --- |
| `og_idx` | UInt32 | |
| `genome_idx` | UInt32 | |
| `copies` | UInt32 | Number of genes in the cell (≥ 1). |
| `genes` | List(Utf8) | Display IDs (prefix stripped), in input order. |
| `genes_raw` | List(Utf8) | IDs as written in the input, parallel to `genes`. |

### 4.5 `presence.npy`

A NumPy `.npy` array, `uint8`, with shape `(n_genomes, ceil(n_og / 8))`:

- **Rows** are in `genome_idx` order (not tree order).
- **Bits** are MSB-first within each byte (`np.packbits(..., bitorder="big")`). OG `j` is bit `7 − (j % 8)` of byte `j // 8`.
- **Padding:** the unused low bits of the last byte are 0. No extra bytes are stored. The demo has 5,657 OGs and 708 bytes per row.

```
byte j//8:   bit7  bit6  bit5  bit4  bit3  bit2  bit1  bit0
OG:          8k    8k+1  8k+2  8k+3  8k+4  8k+5  8k+6  8k+7
```

Unpacking with `np.unpackbits(p, axis=1, bitorder="big")[:, :n_og]` gives the dense `genome_idx × og_idx` matrix (`test_presence_bits_roundtrip`).

### 4.6 `multicopy.npz`

This file is written with `np.savez_compressed`. It holds three parallel `uint32` arrays with one entry per present cell that has `copies > 1`:

| Key | Meaning |
| --- | --- |
| `row` | **Tree row** of the genome (`row_order`), not `genome_idx`. |
| `og` | `og_idx`. |
| `extra` | `copies − 1`. |

Single-copy cells are implied by `presence.npy`. Together the two files give the full copy-number matrix: `copies = presence + extra`. The 10k simulation has 731,498 entries.

### 4.7 `annotations.parquet`

The first column is `og_idx` (UInt32, 0..n−1, every OG present). After it comes one `Utf8` column per annotation, in manifest order. Empty strings and IDs missing from the file are stored as null. The server relies on row index == `og_idx`.

### 4.8 `genome_metadata.parquet`

The first column is `genome_idx` (UInt32, every genome). After it comes one column per field: `Float64` for numeric fields and `Utf8` for categorical ones. Missing values are null.

### 4.9 `sequences.parquet`

| Column | dtype | Meaning |
| --- | --- | --- |
| `gene` | Utf8 | **Raw** gene ID from the membership. |
| `genome` | Utf8 | Genome name. Together with `gene` this is the lookup key. |
| `file` | UInt32 | Index into `manifest.fasta`. |
| `offset` | Int64 | Byte offset of the record's `>` line. |
| `length` | UInt32 | Record length in bytes (header and sequence lines). |

### 4.10 `tree.nwk` and `cache/`

`tree.nwk` is the genome tree that was actually used, as text. Both the importer's `_Tree` and the server's `Pangenome._build_tree` walk it in pre-order with children in file order. That shared walk is why the server's leaf counter equals `row_order`.

`cache/msa/<OG>.fa` holds FAMSA alignments made on demand (see [§5.10](#510-msa-and-famsa-on-demand)), with headers rewritten into OrthoFinder style. The server creates the directory when first needed, so the dataset directory must be writable.

---

## 5. Server internals (`server.py`)

### 5.1 Process model

`server.run()` calls `create_app(data_dir, static, famsa)` and then `uvicorn.run(app, host, port)`. This is a single process. The defaults are `--host 127.0.0.1` and `--port 8000`.

All endpoints except `POST /api/orders` are synchronous `def` functions, so FastAPI runs them in its worker thread pool, and several requests can run at once. NumPy releases the GIL in the heavy operations. A module-level `ThreadPoolExecutor` (`_matrix_pool()`, `min(cpu_count, 16)` workers, created lazily under a lock) parallelises inside a single request.

### 5.2 The `Pangenome` class and its startup work

`create_app` builds one `Pangenome` and stores it in `app.state.pg`. `Pangenome.__init__` loads everything the endpoints need:

1. It reads and checks `manifest.json` (`FileNotFoundError` when missing; `ValueError` when `format_version < 2`).
2. It reads `genomes.parquet` and `orthogroups.parquet`, sorted by their index, and builds dicts for name → `genome_idx` and OG ID → `og_idx`. It keeps `row_order` (genome_idx → row) and its inverse.
3. It loads `presence.npy` and reorders it to **tree row order**: `presence_tree = presence[inv]`. `padded_nbytes = nbytes + 1` (see [§5.3](#53-matrix-aggregation)).
4. It builds `perm_tree` and `perm_frequency` (position → `og_idx`) by stable argsort of the rank columns.
5. It sets the track list to `DEFAULT_TRACKS` followed by `manifest.track_columns`.
6. Annotations: it reads `annotations.parquet` and precomputes `_ann_codes[name]`, a uint16 array per categorical annotation where 0 means missing and `1 + i` is category `i`. It also builds `_search_fields`, a list of `(column, frame)` pairs: the `orthogroup` column of the orthogroups frame plus every annotation column.
7. Genome metadata: it builds `_genome_meta`, a list of dicts indexed by `genome_idx`.
8. It loads `multicopy.npz` into int64 arrays `mc_row`, `mc_og` and `mc_extra`.
9. Detail sources:
   - If `manifest.fasta` is set and `sequences.parquet` exists, it loads the whole file into `_seq_lookup[(raw, genome)] = (file, offset, length)`.
   - It sets `msa_dir`, `seq_dir` and `gene_tree_file` from `manifest.detail`.
   - `results_dir` is the parent of `msa_dir` (or of `seq_dir`). It is needed to find `Duplications.tsv`.
10. FAMSA: `_famsa_ok(famsa)` resolves the binary with `shutil.which` and **runs it** on two 4-residue sequences. The result counts only if it comes back as two equal-length, non-empty records. A broken binary therefore shows up as `detail.famsa = false` and is never used.
11. `_build_tree()` converts `tree.nwk` into the pre-order arrays of `/api/tree` (`parent`, `length`, `label`, `row`).
12. The genomes JSON cache: `_build_genomes_json()` serialises the `/api/genomes` list once, sorted by row, using compact separators and `ensure_ascii=False`. The endpoint returns these bytes directly. At 10k genomes the body is about 1.26 MB, and serving it takes about 1.5 ms instead of re-encoding it each time.
13. Prewarm:
    - `_packed("tree")` builds the default order's packed matrix **synchronously**, so the first screen is fast.
    - A daemon thread builds `frequency` in the background (`_prewarm`), so startup does not wait for it.

### 5.3 Matrix aggregation

The request is `/api/matrix?order&og0&og1&g0&g1&w&h&value`. The window and binning semantics are defined in [api.md](api.md) §Coordinates. The bin edges are computed by `_edges(n, k)`, which returns `floor(i·n/k)` for `i = 0..k`.

#### Per-order packed copies

`_build_packed(perm)` builds a `(n_genomes, padded_nbytes)` uint8 matrix in which column **position** `p` holds the bit of OG `perm[p]`, with rows in tree order. That turns any order into contiguous bit ranges.

The rows are split into chunks of `_ROW_CHUNK = 48`. Each chunk is unpacked, its columns are gathered with `dense[:, perm]`, and it is repacked. The chunks run on `_matrix_pool()`. Small chunks keep the peak memory per task at about 48 × n_og bytes and balance the load across threads.

The extra trailing **zero pad byte** per row lets the last bin edge index byte `og1 >> 3` even when `og1 = n_og` is a multiple of 8.

#### Popcount prefix sums: `matrix_counts()`

For output column edges `e_i = og0 + floor(i·(og1−og0)/w)` (absolute bit positions), let `b_i = e_i >> 3` and `r_i = e_i & 7`. Only the byte span `[bmin, bmax] = [b_0, b_w]` of the window is touched:

```
 packed row (bytes bmin..bmax)       cs = exclusive prefix sum of popcount(byte)
 ┌────────┬────────┬────────┬───     cs[k] = set bits in bytes bmin .. bmin+k−1
 │10110010│01110000│11111111│ ...
 └────────┴────────┴────────┴───     bits before edge e (from bit 0 of byte bmin):
      ▲ e_i = 8·b + r                  prefix(e) = cs[b−bmin] + popcount(byte[b] & MASKS[r])
                                       MASKS[r] = the top r bits (MASKS[0] = 0)
 count in column bin i = prefix(e_{i+1}) − prefix(e_i)
```

For a row slice this is `np.bitwise_count`, `np.cumsum` (int32), one gather at `w+1` byte positions, and a difference. The per-row column sums are then reduced into row bins with `np.add.reduceat(colsum, row_edges[:-1], axis=0)`. No Python loop touches individual cells.

#### Threading across row bins

When `h > 1` and the window has at least 512 genome rows, the `h` output rows are split into `min(16, h)` contiguous groups. Each group's `work(ha, hb)` processes only its genome rows `[row_edges[ha], row_edges[hb])` and does its own cumsum. Threads therefore write to separate output slices, and the result is bit-identical to the single-threaded version (`test_matrix_aggregation_matches_bruteforce`, `..._wide_...`).

#### Encoding

- **Presence** (`matrix()`): `frac = count / (rows_in_bin × cols_in_bin)`, then `rint(frac × 255)` as `<u1`. At full detail the value is 0 or 255.
- **Copies** (`matrix_copies()`): see [§5.4](#54-copy-number-mode).

The response is `out.tobytes()` with the headers `X-Width` and `X-Height`.

### 5.4 Copy-number mode

`matrix_copies()` reuses `matrix_counts()` for the number of present cells per bin. If the window has no present cells it returns zeros. Otherwise it adds the sparse extras:

1. Map each multicopy entry's `og_idx` to its position under the order through the inverse permutation (`_inverse_perm`, O(n_og) per request).
2. Keep the entries inside the window. `mc_row` is already a tree row, so it can be compared with `g0` and `g1` directly.
3. Find each entry's bin with `searchsorted(edges, side="right") − 1` on both axes.
4. Sum the `extra` values per bin with `np.bincount(weights=...)`.

`_encode_copies(counts, extra)` returns `rint(10 × (counts + extra) / counts)`, clipped to 255. Bins with no present cell are 0. The arithmetic is float64 on purpose: float32 moves values across `.5` rounding boundaries (`test_encode_copies_rounding_boundary`).

### 5.5 Orders and caches

`order_perm(order)` resolves an order name to a permutation:

| Order | Source | Cached? |
| --- | --- | --- |
| `tree`, `frequency` | `perm_tree`, `perm_frequency` from the rank columns | Permanently in memory. |
| `track:<name>` and `track:<name>:desc` | `_track_perm()`: stable argsort of the non-NaN values, reversed for `:desc`, then the NaN OGs in `og_idx` order | **Not** cached; recomputed on every call (an O(n log n) argsort). |
| `contrast:<id>`, `search:<id>`, `pattern:<id>` | The dynamic LRU `_dynamic` | Evicted entries return 404. |

With `:desc`, ties come out in reverse `og_idx` order because the stable ascending result is simply reversed.

**Packed matrices** live in two coordinated LRUs, both of size `_LRU_SIZE = 8`:

- `_packed_cache` holds the non-dynamic orders: `tree`, `frequency` and every `track:*` order. These are ordinary LRU entries, so `tree` itself can be evicted after eight other orders have been used. It is then rebuilt on the next request.
- `_dynamic` holds one entry per dynamic order: `{"perm", "packed", "lock"}`. The permutation and its packed matrix are evicted **together**, so a perm can never be alive while its packed matrix is gone, or the reverse ("no 200/404 skew").

Locking:

- One global `_lock` guards all cache bookkeeping (`OrderedDict` lookups, `move_to_end`, inserts and evictions). It is never held during a build.
- **Non-dynamic orders:** `_packed()` first checks the cache under `_lock`. On a miss it computes the permutation, which raises 404 for unknown names before any lock is created, and takes the per-order build lock from `_order_lock()`. It checks the cache again, builds, and inserts. Concurrent requests for the same cold order therefore build it once. Because unknown names are rejected first, `_order_locks` stays bounded by the built-ins plus the real track orders (`test_order_locks_bounded_under_churn`).
- **Dynamic orders:** `_packed_dynamic()` holds the **entry's own** lock while building. If the entry is evicted and re-created during the build, the waiter keeps using its own reference, and the finished matrix is also copied into the live entry when one exists (`test_dynamic_order_eviction_recreation`).

**Order IDs** are minted by `POST /api/orders`. Each body is canonicalised before hashing, and the ID is the prefix plus the first 12 hex digits of a SHA-1:

| `by` | Hashed value | Permutation |
| --- | --- | --- |
| `contrast` | `json.dumps(sorted(set(genomes)))`, so order and duplicates don't matter | `contrast_perm()`: score = fraction present inside the selection − fraction present outside it, then `argsort(−score, kind="stable")`. Ties come out in **`og_idx`** order, not tree order. The selection must contain 1..n−1 distinct genomes. |
| `search` | `json.dumps(q)` | `search_perm()`: hits first, then the rest, each part in tree order. |
| `pattern` | the OG ID | `pattern_perm()`: the OG's pattern group first, then the rest, each part in tree order. |

If the ID already exists, the server returns it and touches the LRU without recomputing (`test_repost_reuses_warm_cache`). Otherwise it computes the permutation outside the lock, then inserts the entry under the lock with a re-check, and evicts the oldest entry if needed. The packed matrix is built lazily on the first `/api/matrix` call. That first call costs about 230 ms at 10k genomes (see [§9](#9-performance-and-scale)).

`contrast_perm` unpacks the **selected rows only**, `nsel × n_og` bytes, and takes the outside count as `n_genomes − inside`, using the `n_genomes` column. It never unpacks the full matrix.

### 5.6 Tracks and the annotation track

`track_means(order, og0, og1, w, names)` handles the track-major response. For each track it:

1. orders the column, `values[perm][og0:og1]`;
2. masks the NaNs;
3. runs `np.add.reduceat` over the bin starts for both sums and counts (float64);
4. divides where `count > 0` and leaves NaN elsewhere.

The output is float32 and is sent as `<f4` with an `X-Tracks` header. Without `names`, every track is returned. An unknown name returns 404.

`annotation_track()` takes the precomputed codes in order, then loops over the `w` bins in Python. Each bin runs `np.bincount(seg, minlength=n_cat+1)[1:]` and `argmax`, so ties go to the lower category index, which is the more frequent category overall. A bin with no codes is 0. The response is `<u2` with the headers `X-Width` and `X-Categories: <name>`. A name that is not a categorical annotation returns 400, whether it is a text annotation or does not exist at all.

### 5.7 Search

`_field_mask(frame, name, q)` runs a polars **expression**, `str.to_lowercase().str.contains(q, literal=True).fill_null(False)`, over one column. The query is lowercased once. Using expressions on a DataFrame rather than Series methods is a leftover of the Python 3.15 issue (see [§11.3](#113-python-314-and-the-polars-on-315-history)).

`search(q, order, limit)` goes through every `(field, frame)` in `_search_fields` and adds one hit dict **per matching (OG, field) pair**. Each hit contains `position = inv[og_idx]` and the matched text. The hits are sorted by position (Python's stable sort keeps the field order for equal positions) and sliced to `limit`. The endpoint clamps `limit` to 1..500. `total` counts (OG, field) pairs, so an OG that matches in two fields counts twice.

`_hit_idxs(q)` ORs the masks into one boolean array of OGs. `search_perm` uses it.

### 5.8 Pattern groups

`/api/orthogroup/{og}/pattern` returns the first `limit` OGs (clamped to 1..5,000) with the same `pattern_id`, in `og_idx` order. The list includes the OG itself, and the frontend filters it out. Each entry carries the OG's annotations. `total` is `n_same_pattern`.

### 5.9 Orthogroup detail and row-group pruning

`GET /api/orthogroup/{og}` works as follows:

1. It finds the row with a polars `filter(orthogroup == og)` on the orthogroups frame. This is a linear scan, not the `_og_to_idx` dict.
2. It reads tracks (NaN becomes null), annotations (`_annotation_values`, using row index == `og_idx`), the pattern fields and the members.
3. `og_members(og_idx)` runs `pl.scan_parquet(membership).select(...).filter(og_idx == x).collect()`. Polars pushes the predicate into the Parquet reader and skips row groups whose `og_idx` statistics exclude `x`. Because membership is sorted by `og_idx` and written in 10k-row groups, a lookup reads only the one or two row groups that contain the OG and not the whole file (290 MB at 10k).
4. Members are returned sorted by `genome_idx`, which is the input column order and **not** tree order.
5. The body is serialised with `json.dumps` and returned as a `Response`. For a 10k-member OG, FastAPI's `jsonable_encoder` was far slower.

`_member_context(og_idx)` is shared by the MSA, gene-tree and download paths. It returns the members plus:

- `member_map`: genome → raw genes;
- `raw_to_display`;
- `expected`: the header map from `tracks.msa._expected_headers` (`{"<genome>_<raw>" and "<genome with . → _>_<raw>": {(genome, raw)}}`).

This is the same exact, never-by-prefix matching that the MSA tracks use (see [§6.3](#63-msapy)).

### 5.10 MSA and FAMSA on demand

`og_msa(og)` calls `_msa_text()`, which tries the sources in this order:

1. `msa_dir/<OG>.fa`, the OrthoFinder alignment. Source `"orthofinder"`.
2. `cache/msa/<OG>.fa`, an earlier FAMSA result. Source `"famsa"`.
3. A new alignment, which runs these checks and steps in order:
   - No sequences → 404.
   - More than `_FAMSA_MAX_SEQS = 2000` sequences → 404.
   - FAMSA unusable → 404.
   - `_famsa_sem.acquire(blocking=False)` fails → **503** (`_FAMSA_CONCURRENCY = 2`). Requests never wait for a semaphore slot; a same-OG request can still wait on the per-OG lock after acquiring one.
   - Take the **per-OG** lock (`_famsa_lock`) and check the cache again, so two requests for the same OG run FAMSA once (`test_famsa_same_og_single_run`). A same-OG request that waits here already holds a semaphore slot, so it can occupy one of the two slots while blocked rather than freeing it.
   - `_run_famsa()` writes a temporary input in `cache/msa/` with **synthetic unique headers** `m0..m{n−1}`. Real gene IDs can contain characters FAMSA might mangle, and duplicates must stay distinct. It runs `famsa IN OUT -t 1` (OrthoFinder's command) with `timeout=_FAMSA_TIMEOUT = 120` s. A timeout or non-zero exit becomes a 500 with a short stderr excerpt. The temporary files are always deleted.
   - `_canonicalize_famsa()` validates the output. The headers must be an exact permutation of `m0..m{n−1}`, no record may be empty, and all records must have the same length. If any check fails it returns 500 and **nothing is cached**. Otherwise it rewrites the headers to OrthoFinder's `<genome with . → _>_<raw>` form.
   - `_write_cache()` writes to a unique `mkstemp` file and calls `os.replace` onto `<OG>.fa`, so readers never see a partial file.

`og_sequences(og_idx)` gives the unaligned input. It reads `Orthogroup_Sequences/<OG>.fa` if that exists, whose headers are raw IDs. Otherwise it reads each member's `_seq_lookup` record from the `--fasta` files. Rows are sorted by genome row, then by raw ID.

**Header mapping.** The MSA text is parsed, and `_match_headers(headers, expected)` assigns each header to `(genome, raw)` if the match is unique. Each row gets the genome, `row = _row_of(genome)`, the display `gene`, `raw` and `seq`. Unmatched or ambiguous headers keep `genome = null`, `row = −1`, and the header as gene and raw. Rows are sorted by `(row, or 2³⁰ when unknown; display gene)`. `trimmed` comes from the `trimmed` track column (`_trimmed_value`): null when the column is missing or NaN, otherwise a boolean.

### 5.11 Gene trees and duplication flags

- **Index:** `_gene_tree_index()` scans `Resolved_Gene_Trees.txt` once, lazily, under `_gt_lock`. It records `{og: (byte offset, line length)}` for every line containing `:`. `_gene_tree_raw(og)` seeks to the line and returns the Newick after the first `:`.
- **Duplications:** `_duplications(og)` reads `results_dir/Gene_Duplication_Events/Duplications.tsv` once, lazily, under `_dup_lock`, using `read_matrix` (its header starts with `Orthogroup`). It builds `{og: {gene-tree node label: support}}` from the `Gene Tree Node` and `Support` columns. The first event per node wins, and a non-numeric support becomes null.
- **Arrays:** `_gene_tree_arrays()` walks the parsed tree in pre-order, like `/api/tree`.
  - `duplication[i]` is true **iff the node's label appears in the duplications table**. The support value may still be null.
  - A leaf whose label maps uniquely through `expected` gets the display `gene` and its genome row. Other leaves get `gene = null` and `genome_row = −1`.

### 5.12 Downloads and `Content-Disposition`

- `download_fasta(og, aligned, ids)`:
  - `aligned=1` uses `og_msa`, so it may trigger FAMSA and may return 503.
  - `aligned=0` uses `og_sequences`; with no sequences it returns 404.
  - `ids=raw|display` chooses the header.
- `download_newick(og, ids)` re-serialises the parsed gene tree with `_to_newick(root, name_of)`. Leaves are renamed to raw or display IDs when they map uniquely. Internal labels are kept (`test_newick_round_trip_preserves_internal_labels`). Branch lengths are written with `:g` formatting, which keeps six significant digits.
- `_content_disposition(filename)` rejects control characters with 400. It emits `attachment; filename="<ASCII fallback>"; filename*=UTF-8''<percent-encoded>`; in the fallback, non-ASCII characters become `?`, and `\` and `"` become `_` (`test_download_filename_encoding`).

### 5.13 Request limits

| Limit | Value | Where |
| --- | --- | --- |
| `w`, `h` per axis | 1..`MAX_CELLS = 4096`, and `≤` the window size | `_window()`, tracks, annotation-track |
| `w × h` | ≤ `MAX_CELLS_TOTAL = 4_000_000` | `_window()` |
| OG ID list | ≤ `MAX_LIST = 5000` per request | `/api/orthogroups` |
| Pattern list | `limit` clamped to 1..5000 | `/api/orthogroup/{og}/pattern` |
| Search | `limit` clamped to 1..500 | `/api/search` |
| `POST /api/orders` body | ≤ `MAX_BODY = 64 KiB`, checked against `Content-Length` **and** while streaming chunks; otherwise 413 | `create_order` |
| Contrast selection | 1..n−1 distinct genomes, no more names than genomes | `create_order`, `contrast_perm` |
| FAMSA | ≤ 2,000 sequences, 2 concurrent alignments, 120 s timeout | `_msa_text`, `_run_famsa` |

Out-of-range windows and unknown `value` modes return 400 (`test_bad_requests_400`). There is no authentication and no global rate limit (see [§11.6](#116-localhost-only-server)).

### 5.14 The static mount

At the end of `create_app` the server picks a static directory:

1. `--static PATH` if given;
2. otherwise `src/orthopan/static/` inside the installed package, if it exists;
3. otherwise `<repo>/web/dist`, a legacy fallback;
4. otherwise none, and the server is API only.

The chosen directory is mounted at `/` with `StaticFiles(html=True)`. The mount is registered **after** all `/api` routes, so API paths always win. Vite builds with `base: './'` (relative asset URLs) into `src/orthopan/static`, so the UI works from any mount path.

---

## 6. Track preprocessors (`tracks/`)

Each preprocessor turns its inputs into `{og: {column: float | None}}` and writes a track TSV with `tracks.write_tracks()`. Values are formatted with `.6g`, and None or non-finite values become empty cells. Each module also exposes a `build_tracks` or `basic_tracks` function that can be called in-process. The file format is described in [design-decisions](design-decisions.md) §Per-orthogroup tracks.

### 6.1 `basic.py`

`basic_tracks(path, fasta=())` reads the matrix with `read_matrix` and outputs `paralog_fraction` and `length_cv`:

- **Copies per present cell** are `cell.count(",") + 1`.
- **Paralog fraction** is the share of present genomes with more than one copy.
- **Length CV** is computed by `_length_cv`. It is the population CV of member protein lengths, and is None with fewer than two known lengths or a zero mean.
- **Lengths** come from `load_lengths()`. It counts characters in each FASTA record, ignoring whitespace and `*`, and keys them by the first word of the header.
  - An ID that appears with **different** lengths in different files is ambiguous. It is re-keyed as `(genome, gene)` using the file stem → genome label map.
  - Ambiguous IDs that cannot be mapped are dropped, with a warning.

Copy number is deliberately absent: `orthopan import` always writes `copy_number_mean`/`copy_number_max`, so a `basic.tsv` carrying those names would collide (`duplicate track column name`). With them gone, `basic.tsv` can be passed straight to `--track`.

### 6.2 `orthofinder.py`

`build_tracks(results_dir)` takes the OG universe from `Orthogroups/Orthogroups.tsv`.

- **`duplications` and `duplications_supported`** come from `Gene_Duplication_Events/Duplications.tsv`. They count events per OG, and the supported count uses `Support >= 0.5`. Rows that are too short or have an empty `Support` are skipped, with a warning.
- An OG that is missing from the duplications table is **0** if it has a gene tree, because OrthoFinder only writes events for OGs with trees. It is **empty** otherwise.
- **`tree_length` and `tree_length_per_leaf`** come from `Resolved_Gene_Trees.txt`. The branch lengths are the regex matches of `:number` (signed, scientific notation allowed), and the leaf count is `commas + 1`.
- Missing files leave their columns empty.

### 6.3 `msa.py`

`msa_metrics(seqs, weights=None)` computes all four MSA tracks for one alignment.

**Encoding.** `_NONSTANDARD` maps every byte to a code:

- `A`–`Z` and lowercase letters become 0..25, except `O` and `U`;
- `-` becomes 26 (the gap);
- everything else becomes `X`.

All-gap columns are dropped first, because HOG pruning creates them. `gap_fraction` is the share of gap cells in the cleaned alignment. With fewer than two sequences, only `gap_fraction` is set.

**Geometric homogeneity: the closed form with weights.** Let `w_g` be the weight of sequence g, `G = Σ w_g`, `S₂ = Σ w_g²`, `k_g` the gap count of sequence g, `W_c = Σ_g w_g·gap_{g,c}` the weighted gap count of column c, and `L` the number of columns (needs L ≥ 2):

```
by_residue = Σ_g w_g [k_g(k_g−1) + (L−k_g)(L−k_g−1)] / (G · L · (L−1))
by_gene    = Σ_c [W_c² + (G−W_c)² − S₂] / (L · (G² − S₂))
geometric  = (by_residue + by_gene) / 2
```

With `w = 1` these reduce exactly to anvi'o's pairwise-XOR result: `W_c² + (n−W_c)² − n = m_c(m_c−1) + (n−m_c)(n−m_c−1)`. The derivation and the anvi'o reference are in [anvio-internals.md §1](anvio-internals.md#1-geometric-homogeneity). The by-gene term weights each ordered pair (a, b) with a ≠ b by `w_a·w_b`.

**Functional homogeneity.** `_SCORE` is a 26 × 26 pair-score table built from anvi'o's conserved groups (`_CONSERVED`):

- 3 for identical standard residues;
- 2 for a conserved pair, and for `B/B`, `Z/Z` and `J/J`;
- 0 for `X/X` and for everything else.

Per column, `T[c, x] = Σ w` and `U[c, x] = Σ w²` over the sequences with residue x. Then:

```
score_c = ½ (Tᵀ·S·T − Σ_x U_x·S_xx)        # weighted unordered non-gap pairs
max_c   = 1.5 ((Σ_x T_x)² − Σ_x U_x)        # = 3 × weighted pair count
functional = Σ score_c / Σ max_c   over columns with ≥ 2 non-gap sequences
```

The non-gap totals come from the residue counts directly, rather than "total minus gaps", to avoid catastrophic cancellation with many small paralog weights. The result is None when no column has a pair, and it is clamped to [0, 1] only to absorb float noise. The anvi'o algorithm is in [anvio-internals.md §2](anvio-internals.md#2-functional-homogeneity). With all weights equal to 1, both scores equal anvi'o's (`test_pooled_*_matches_anvio`), using anvi'o's own `homogeneityindex.py` loaded from `resources/anvio`.

**Genome weighting.** With `weighting="genome"`, the default, each MSA header must match exactly one `(genome, gene)` from that OG's `Orthogroups.tsv` row, through `_expected_headers` and `_match_headers`. A genome with k copies then gives each copy weight `1/k`. If any header is unmatched or ambiguous, the OG falls back to pooled weights (all 1), and one summary warning reports the number of fallbacks. Without `--orthogroups` the whole run falls back to pooled. The reason for per-genome weighting is in [design-decisions](design-decisions.md) §Homogeneity.

**The trimmed tri-state.** This needs `--seqs` (`Orthogroup_Sequences/`). `_trimmed()` compares each aligned sequence's non-gap count with its unaligned length. The unaligned length comes from `_unaligned_lengths`, with trailing `*` stripped; conflicting duplicate IDs become None.

- The gene comes from the exact header match. Otherwise `_unique_suffix_gene` tries a **unique** suffix match after any `_`.
- **1.0** if any matched sequence is shorter in the MSA than unaligned.
- **0.0** only if every sequence matched and none is shorter.
- **None** otherwise, and also when there is no `--seqs` file.

### 6.4 `msa.py` parallelism

`build_tracks(msa_dir, seqs, orthogroups, threads=cpu_count, weighting)` streams the `Orthogroups.tsv` rows. It builds members only when they are needed: genome weighting, or `--seqs` for exact matching. Then it appends the MSA files that are not in `Orthogroups.tsv`, sorted, with empty members.

With more than one thread it uses a `ProcessPoolExecutor`:

- The platform default start method is used, because forking after native threads can deadlock. `_init` passes the config explicitly to each worker.
- Work is submitted in batches of 32 OGs (`_process_batch`).
- At most `2 × threads` futures are in flight. `flush()` waits for `FIRST_COMPLETED` whenever the limit is reached.

At no point is the per-OG input held beyond a batch, but `build_tracks` keeps the result dict for **every** orthogroup in memory (that is its return value), so its footprint is O(n_og) plus the bounded worker batches. The sequential path calls `_process` directly. `test_parallel_matches_sequential` checks that both give the same result. The module prints a timing summary to stderr.

---

## 7. Frontend architecture (`web/src`)

This is plain TypeScript with no framework. deck.gl draws the heatmap. Everything else is Canvas2D or DOM.

### 7.1 Module map

| Module | Responsibility |
| --- | --- |
| `main.ts` | Wiring and all mutable UI state: boot, fetch lifecycle, render loop, interaction, tooltips, drawer, alignment-view orchestration, URL hash, keyboard. |
| `theme.ts` | Theme resolution (`auto`/`light`/`dark` → `light`/`dark`), the light/dark `Palette`, `palette()`, and the live store (`initTheme`, `setThemeSetting`, `onThemeChange`). The single source of colour for the Canvas2D panels and deck.gl. |
| `transform.ts` | The shared world↔pixel transform, view clamping and zoom, window and cell math, `prepareCanvas` (DPR). Pure functions apart from `prepareCanvas`. |
| `api.ts` | Typed client for every endpoint, plus binary payload size checks, `HttpError` and `NetworkError`. |
| `decode.ts` | Little-endian `Uint16`/`Float32` decoders and the copy-number byte maths. |
| `matrix.ts` | The themed presence colormap (gamma LUT), the copies palette, `toImageDataFor`, the tile maths, the `TiledBitmapLayer` rounded-cell shader, and the `Heatmap` class (deck.gl `Deck` + layer, with a flat fallback). |
| `tree.ts` | `Tree`: the genome dendrogram on Canvas2D (phylogram or cladogram), clade selection and hit testing. |
| `tracks.ts` | `Tracks`: numeric bar strips and categorical colour strips above the matrix, `robustDomain`, strip hit testing, `fmt3`. |
| `annotations.ts` | The categorical palette (12 colours, then grey) and value → colour and label helpers. |
| `metadata.ts` | `Metadata`: the genome metadata columns between tree and heatmap, numeric gradient. |
| `search.ts` | `SearchBox` (debounced search with stale-result guards) and `hitViewport`. |
| `selection.ts` | The selected-genome model (single genome or binned row range), hash encoding, member summary. |
| `members.ts` | `memberSlice`: pagination of the member list, with selected genomes pinned. |
| `pending.ts` | `PendingCache`: a bounded memo cache that shares in-flight promises. |
| `generation.ts` | `Generation` (counter + `AbortController`) and `Sequence` (bare token) for stale-response guards. |
| `requestkey.ts` | `requestKeys`/`needsFetch`: dedupe keys for the matrix window and its track/annotation layers. |
| `dropdown.ts` | `placeDropdown`/`trackDropdown`: position-fixed dropdowns anchored to a top-bar control. |
| `toast.ts` | The toast queue (dedupe, max-3 eviction, pause/resume) and the bottom-centre container: `hint`/`notice`/`error`. |
| `alignment.ts` | Alignment-view maths (row order, virtual band, neighbour stepping) and the `AlignmentView` renderer. |
| `msa.ts` | MSA column maths (zoom and scroll clamping, visible columns), residue colours, conservation. |
| `genetree.ts` | Gene-tree layout (`layoutGeneTree`), `treeNodeX`, leaf colours. |

### 7.2 The shared world↔pixel transform

Everything aligned with the matrix uses one coordinate system:

- **world X** is an OG position, and OG `p` covers `[p, p+1)`;
- **world Y** is a genome row, and row `r` covers `[r, r+1)` with its centre at `r + 0.5`.

A `View {cx, cy, zx, zy}` holds the world point at the centre of the matrix viewport and the **log2** pixels per unit for each axis, so the axes zoom independently. `makeTransform(view, pw, ph)` returns `xOf`, `yOf`, `ogAt` and `rowAt` with `sx = 2^zx` and `sy = 2^zy`, where `(pw, ph)` is the CSS size of `#matrix-wrap`.

```
 tree canvas      meta canvas   tracks canvas (x only) / #matrix-wrap        labels canvas
 yOf(r) rows      yOf(r) rows   xOf(p), yOf(r)                               yOf(r+0.5)
      └──────────────┴─────────────── makeTransform(view, pw, ph) ───────────────┘
                                          ▲
                deck.gl OrthographicView {target:[cx,cy], zoomX:zx, zoomY:zy}
```

- **The heatmap:** `Heatmap` drives deck.gl's `OrthographicView` with the same numbers (`viewStateOf`: `target = [cx, cy]`, `zoomX = zx`, `zoomY = zy`). Each fetched bitmap is placed at its own world rectangle, `bounds = [og0, g1, og1, g0]`. Old data therefore stays in its correct place while the view moves and a new window loads. Textures use `nearest` filtering.
- **The Canvas2D panels** (`Tree.draw`, `Tracks.draw`, `Metadata.draw`, `drawLabels`, `drawOverlay`) all take the same `Transform`.
- **Pointer picking:** `bufferAt(win, tr, px, py, …)` maps a pointer to a buffer cell **within the fetched window** and returns null outside it, so stale data is not interactive. `cellRect` uses the bitmap's continuous texel edges for the hover outline. `colRange`, `rowRange` and `binRange` use the server's floor-rounded bin edges for tooltip text.

Why everything must use it: the heatmap and each Canvas2D layer are separate surfaces. If one panel computed its own mapping, for example with rounding to integer pixels, its rows and columns would drift against the bitmap, especially at fractional zoom or fractional DPR. `transform.test.ts` (`bitmap <-> label alignment`, `transform round-trip`) pins this down. New panels must take a `Transform` and must not keep their own scale.

`clampView` and `zoomedView` enforce zoom between "whole axis fits" and 40 px per unit (min 0.01 px per unit). `zoomedView` clamps the scale **before** it derives the anchor offset, so zooming at a limit does not move the view. `computeWindow` returns the visible half-open window and the `w × h` request: about one output cell per CSS pixel, capped at `MAX_CELLS = 4096` per axis and scaled down by `capCellCounts` to stay within `MAX_CELLS_TOTAL`. These constants mirror the server's limits.

### 7.3 DPR handling

`prepareCanvas(canvas, cssW, cssH)` is the **only** way Canvas2D panels set up drawing:

1. It sets the backing store to `round(css × devicePixelRatio)` (`backingSize`), and only touches `width`/`height` when the size changed.
2. It calls `ctx.setTransform(dpr, 0, 0, dpr, 0, 0)`.
3. It clears the canvas.

Callers then draw in CSS pixels. Each canvas's CSS box is set to 100% by `index.html`.

A DPR change at the same CSS size, for example moving the window to another monitor or zooming the browser, does not fire `resize` or the `ResizeObserver`. `main.ts` handles it in two ways:

- `watchDevicePixelRatio()` registers a `matchMedia("(resolution: <dpr>dppx)")` change listener and re-arms it with the new value on every change.
- `checkDevicePixelRatio()` polls every 400 ms as a safety net, because some environments, notably CDP device-metric emulation, update the query without firing `change`.

Both call `onDprChange()`, which runs `heatmap.resyncDpr()` (`useDevicePixels` plus a redraw), then `render()`, then `alignmentView.draw()`. It does not refit the view or fetch again.

### 7.4 Data-fetch lifecycle

```
view/order/mode/layers change
   └─ scheduleFetch(): gen++, abort `inflight`, debounce 140 ms
        └─ fetchData(): snapshot (gen, order, mode, window), new AbortController
             0. requestKeys(order, value, window, tracks, annos) ─▶ needsFetch(lastMatrixKey, lastTracksKey, hasMatrix)
                (unchanged key with a displayed matrix ⇒ nothing to fetch)
             1. in parallel, each with its own stale check (needed ones only):
                ├─ getMatrix ──▶ (stale? return) ──▶ lastMat/lastWin/dataOrder, setData, render
                │       404 on a non-tree order ──▶ notice toast + setOrder('tree')
                ├─ getTracks ∥ getAnnotationTrack×N ──▶ (stale? return) ──▶ lastTracks/lastAnnos, render
                └─ ensureOgNames (only if og1−og0 ≤ 5000) ──▶ activeOgs (labels, ogAt)
```

**Request dedupe** (`requestkey.ts`). `requestKeys(order, value, win, tracks, annos)` builds `matrix = "<order>|<value>|og0|og1|g0|g1|w|h"` and `tracks = matrix + "|<tracks>|<annos>"`. `fetchData` keeps `lastMatrixKey`/`lastTracksKey`; `needsFetch` skips the matrix fetch when the key is unchanged and a matrix is already displayed (`hasMatrix = lastMat !== null && dataOrder === order`), and skips the track fetch only when the composite track key is unchanged. A layers-only change fetches tracks and annotations without refetching the matrix. `setOrder` clears `dataOrder` and both keys, so the next call always refetches. This stops a tracks-only change from re-requesting the matrix; a theme change never calls `fetchData` at all.

**Transport errors** (`api.ts`). `request()` wraps `fetch` and throws `NetworkError(url, cause)` when there was no HTTP response; aborts are re-thrown unchanged. `toastError` turns an `HttpError` into `"<context> failed: <detail>"` keyed by `errorIdentity(context, status, message)`, a non-HTTP client-side failure into a context-keyed error, and **collapses every `NetworkError` to one `"server unreachable"` toast** keyed `network`, with the URL and message in the details.

**Stale-response rules.**

- Every async result is dropped unless its captured generation is still current and its signal was not aborted (`myGen !== gen || ac.signal.aborted`).
- The OG names in `ensureOgNames` also require `snapOrder === order`.
- Hover and click do nothing while `dataOrder !== order` or `lastMat.value !== mode`, so a cell is never resolved against the previous order's bitmap.
- `setOrder()` clears every order-dependent value (`activeOgs`, `ogCache`, `dataOrder`, the strips and hover), invalidates the search box and `navSeq`, and schedules a render, a hash update and a fetch.

**Generations and sequences** (`generation.ts`):

| Guard | Type | Protects |
| --- | --- | --- |
| `gen` + `inflight` | counter + `AbortController` in `main.ts` | Matrix, track and label fetches. |
| `drawerGen` | `Generation` | The Members pane and the pattern list: `beginDrawer()` aborts the previous OG's requests. |
| `alnGen` | `Generation` | The alignment view's MSA and gene-tree loads (`beginAlignment()`). |
| `navSeq` | `Sequence` | Queued ←/→ steps. It is invalidated by close, route change and order change, and the steps are chained through `navChain` so they never overlap. |
| `hoverSeq`, `clickSeq` | counters | The async hover annotation and click OG lookup. |
| `SearchBox.seq` | counter + controller | Search results. Both query and order must still match when the response arrives, and a click on a stale result list is ignored. |

**Caches and their bounds.**

| Cache | Bound | Notes |
| --- | --- | --- |
| `ogCache` | 8 windows | OG ID lists per `order:og0:og1`. Cleared on order change. |
| `ogAnnoCache` (`PendingCache`) | 64 OGs | Hover annotations, which come from `/api/orthogroup/{og}`. Only one hover request is in flight at a time (`annoHover`): hovering a different OG aborts it and calls `forget()`, so the aborted promise is never cached. Errors other than an abort are cached as `{}`. |
| `msaCache` | 8 OGs (`ALN_MSA_CACHE_MAX`) | Successful MSA responses, in LRU order. |
| `treeState` / `treeData` / `treeError` | 16 OGs (`TREE_CACHE_MAX`), evicted together | Gene trees. Only a **404** is cached as `missing`. Aborts and transient errors stay retryable. |
| `treeInflight` | cleared on every new drawer or alignment generation | Shares an in-flight gene-tree request between the drawer's Newick-link probe (`probeTree`) and the alignment view. It is cleared because those requests were aborted. |

`PendingCache.load(key, loader)` returns the cached value, or the pending promise, or starts the loader. A loader result is cached only if its promise is still the registered pending one, so `forget()` wins. The oldest entries are evicted past `limit`.

### 7.5 Heatmap rendering

`toImageDataFor(mode, values, w, h, theme)` expands the server's one-byte bins into an RGBA `ImageData`, using the theme's palette.

- **Presence:** a 256-entry LUT built from `displayT(v)`. Absent (`v = 0`) is the first stop — near-white in light, near-background dark navy in dark. Non-zero values map to `0.12 + 0.88·(v/255)^0.45` on the theme's five-stop ramp. The gamma and the floor keep sparse bins visible at low zoom. `colorForValue` drives the legend gradient.
- **Copies:** `copiesBucket(v) = clamp(round(v/10), 1, 5)` for `v > 0` selects one of six discrete colours (absent, 1, 2, 3, 4 and 5+) from `palette().copiesColors`. The tooltip shows `formatCopiesMean` for binned cells and `formatCopiesExact` for single cells. A saturated byte (255) is shown as `≥25.5`.
- The image's **alpha channel is only a presence mask** for the tile shader; the flat path ignores it.

**Rounded tiles.** At high zoom the heatmap draws rounded cells via `TiledBitmapLayer`, a `BitmapLayer` subclass whose fragment shader (`TILE_FS`) applies a per-cell rounded box (`sd_round_box`) from a `tile` uniform block (`tileModule`). `heatmapTileParams()` is passed to `Heatmap.setData` and read on every deck.gl frame:

- Tiles are used only when the window is at full detail on **both** axes (`!binnedX && !binnedY`); any binned window stays flat.
- `tileFade(minCellPx)` is the activation: 0 at ≤ `TILE_MIN_PX = 6` CSS px, 1 at ≥ `TILE_FULL_PX = 10`, smoothstepped between. `mix = 0` reproduces the flat rendering exactly (the shader early-returns).
- Gap ≈ `TILE_GAP_FRAC = 0.12` of the cell, clamped to `[TILE_GAP_MIN, TILE_GAP_MAX] = [1, 3]` CSS px; corner radius ≈ `TILE_RADIUS_FRAC = 0.25` of the smaller side, clamped so it never collapses into the gap. Sizes are multiplied by `devicePixelRatio` for framebuffer pixels.
- Absent tiles are multiplied by `palette().absentDarken` (0.955 in light, 1.0 in dark).
- The hover outline in `drawOverlay` follows the same rounded shape (`roundRectPath`) once `mix > 0.5`.

**Flat fallback.** luma.gl reports shader compile/link failures only through `console.error` (deck's `onError` never fires). The inline `<head>` script in `index.html` hooks `console.error` **before** the bundle captures it and forwards shader-looking messages to `window.__onShaderError`, which `Heatmap` installs. `Heatmap.handleError` sets `flat = true` permanently (exposed as `usingFlat`), warns once, and rebuilds the layer as a plain `BitmapLayer` over an opaque copy of the image. The layer's `onError` calls the same handler.

`Heatmap.setData` replaces a single layer (`id: 'matrix'`). The deck controller only enables `dragPan`. Wheel zoom is handled in `bindInteractions` (plain wheel zooms the OG axis, Shift+wheel the genome axis) and pushed back with `heatmap.setView`. `onViewStateChange` clamps the controller's state and writes it back (controlled mode).

### 7.6 Theming

`theme.ts` is the single source of colour. The user setting is `auto` | `light` | `dark`; `resolveTheme` maps it to `light`/`dark`, with `auto` following `prefers-color-scheme`.

- `document.documentElement.dataset.theme` drives the CSS custom properties in `index.html`, so all DOM chrome (top bar, panels, tooltips, toasts) follows the theme.
- `LIGHT_PALETTE` / `DARK_PALETTE` hold the same values in TypeScript for the Canvas2D panels and deck.gl; `palette()` returns the current one. Each theme has its own presence/numeric ramp and categorical colours (dark uses brighter categories and a navy → cyan → pale-yellow presence ramp).
- **Monotonic ramps.** Both presence ramps are monotonic in luminance, in opposite directions: light runs white → dark blue as presence rises, dark runs near-background → pale yellow. The copies ramps are monotonic too. `theme.test.ts` pins this, so a bin value maps to one colour within a theme.
- **Live switching, no refetch.** `setThemeSetting` persists to `localStorage` (`orthopan.theme`), updates the `data-theme` attribute and notifies listeners. `main.ts` registers `onThemeChange(() => applyTheme())`; `applyTheme` re-rasterises the current matrix image with the new theme (`refreshHeatmapRaster`), re-renders the panels and redraws the alignment view. It never calls `fetchData`.
- **First paint.** An inline `<head>` script reads the same `localStorage` key and `prefers-color-scheme` and sets `data-theme` (plus a `color-scheme` meta) before the bundle loads, so there is no light flash. `initTheme()` re-reads the setting at boot and starts a `matchMedia` listener that re-resolves `auto` on OS changes.

### 7.7 Notifications (toasts)

`toast.ts` is a small bottom-centre message queue. The pure `ToastQueue` (clock-injected) does dedupe, cap and timing; `ToastContainer` does the DOM. `toasts()` lazily creates one container.

- Three kinds with different lifetimes (`TOAST_DURATIONS`): `hint` 3 s, `notice` 6 s, `error` **sticky** until closed (`×`).
- **Hints** with the same text restart the existing toast and its timer instead of stacking.
- **Errors** dedupe by key (`errorIdentity`): `<context>|<status>|<message>` for HTTP errors, `network` for any transport failure, `client|<context>|<message>` for client-side ones. A repeat bumps a visible `×N` count and updates the details in place. Errors carry collapsible `details` (status/url/message) and use `role="alert"`.
- **Notices** may carry an action, for example the 503 alignment toast with a **Retry** button.
- At most `MAX_TOASTS = 3` are shown; the oldest auto-dismissing toast is evicted first, otherwise the oldest. Sticky errors survive while auto toasts remain.
- Hovering a toast **pauses** its countdown; leaving resumes it with at least `RESUME_MIN_MS = 1500` left.
- No keyboard shortcuts: toasts are mouse-only (hover to pause, `×` to close), and the container is `pointer-events: none` except on the toasts, so it never blocks the matrix.

Usage in `main.ts`: `toasts().hint('zoom in until the cell covers a single orthogroup to open it')` on a binned click; `toasts().hint('copied <id>')` / `toasts().error('copy failed')`; `toasts().notice("unknown order '<x>' — using tree")` on a 404 order; `toasts().notice('alignment limit busy, try again', {action: {label: 'Retry', …}})` on a 503 MSA; and `toastError(err, context)` for every failed request (the NetworkError collapse is in §7.4).

### 7.8 Tree, tracks, annotation strips and metadata columns

- **Tree** (`Tree.load` / `draw`): rebuilds children, cumulative branch-length depth, edge depth and leaf-row intervals from `/api/tree` arrays. Only strictly earlier parents are accepted, which guards against cycles. Node y is the centre of its row interval. Node x uses depth (phylogram), or edge count with tips aligned (cladogram, also used automatically when branch lengths are degenerate). `hitTest` picks the nearest horizontal branch within 7 px. A selected node tints its clade on the tree, labels and overlay canvases, and shows "sort by contrast" when the clade has more than one genome.
- **Tracks** (`Tracks.draw`): the strip area (`--track-h`, 68 px) is divided evenly between the active strips.
  - Numeric strips draw bars, scaled by `robustDomain`: the 1st–99th percentile of the **current window's** bin values, falling back to the global min and max from meta.
  - Category strips draw `annotationStripColor(code)`.
  - Each strip remembers the world span it was fetched for (`world`), so stale strips stay in place and `hit()` returns null outside it.
  - By default, `pickTracks` turns on `copy_number_mean`, the first `*homogeneity*` track and every categorical annotation.
- **Metadata** (`Metadata.draw`): one 16 px column (`META_COL_W`) per field, drawn only for visible rows, with rotated headers. Categorical values use the palette by category index. Numeric values use a five-stop gradient between the meta `min` and `max`.
- **Labels and OG names:** genome names are drawn when rows are at least 11 px tall. Below that the label column collapses to 14 px. OG IDs are drawn rotated on the overlay when columns are at least 24 px wide and the names for the window are loaded.

### 7.9 Search

`SearchBox`:

- runs `/api/search` 220 ms after typing stops, or at once on Enter;
- shows at most 50 hits plus the total;
- clears on Esc.

Selecting a hit calls `selectHit`, which centres the hit at 40 px per OG with `hitViewport` and `clampView` and highlights the column for 2 s. "Order by hits" posts `{by: "search"}` and switches to the returned order.

`setQuery()` restores the query from the hash without causing a hash update.

### 7.10 The drawer and Members

A click on a cell (a pointer move of at most 4 px counts as a click) works like this:

1. If the cell aggregates several OGs, a hint toast asks the user to zoom in.
2. Otherwise the OG ID comes from `activeOgs` or from a one-item `/api/orthogroups` request.
3. The genome selection is `singleSelection(name)`, or `rangeSelection(names, [r0, r1])` when the row bin covers several genomes.
4. `openDrawer(og, sel)` runs. For a new OG, it starts a drawer generation, sets the download links, probes for a gene tree (which disables the Newick link if there is none) and renders the members. For the same OG with a new selection, it only re-renders the selection block and the list.

`renderMembers` fetches `/api/orthogroup/{og}` and renders the summary table, then:

- **The selection block** (`summarizeSelection`):
  - a single genome shows "In genome X: N genes" or "absent in X";
  - a range shows how many of the selected genomes have the OG, and lists genes only up to `SELECTION_GENE_CAP = 50`.
- **The member list:** `memberSlice(groups, membersLimit, pinned)` renders the first `MEMBERS_PAGE = 400` genes across genome groups, in server order (`genome_idx`). Selected genomes are always shown in full. "show more" adds another 400.
- **The pattern section:** `loadPattern` lists the other OGs with the same pattern and has an "Order by pattern" button.

Copy buttons copy the **raw** ID through one delegated handler (`toasts().hint` on success, `toasts().error` on failure).

### 7.11 The alignment view

`openAlignment(og)` pushes `#…&og=<OG>&view=aln` and calls `showAlignment`. `loadAlignment(og)` then:

1. starts `alnGen` and clears the previous OG's content;
2. loads the gene tree: from the cache, or as a cached "missing", or with `loadTree`;
3. renders immediately with the cached MSA or a "loading alignment…" placeholder;
4. fetches the MSA. A 503 shows an in-view **retry** button **and** a notice toast with a **Retry** action.

Every await re-checks `alnGen.isCurrent(my) && alnOg === og`.

**Layout and virtualisation** (`AlignmentView`):

```
 #aln-view (fixed, full window)
 ├─ header controls, legend
 ├─ #aln-bar  canvas: labels | conservation bar | gap bar | overview strip (MSA part only)
 └─ #aln-scroll (overflow-y: auto)                        ← the only vertical scroller
     └─ #aln-rows  height = rows × 16 px                  ← provides scroll extent
         └─ #aln-canvas  position: sticky; top: 0; height = viewport height
              x: 0 ─ tree (300) ─ colour (22) ─ name (210) ─ MSA (rest, own horizontal scroll)
```

- The canvas is only as tall as the viewport. `visibleBand(scrollTop, viewH, n)` selects the rows to draw, with 2 rows of overscan. A 2,000-row alignment therefore never needs a canvas larger than the screen.
- Horizontal MSA scrolling is virtual: `scrollCol` and `pxPerCol` are clamped by `clampScrollCol` and `clampPxPerCol` (0.25–40 px per column).
  - Shift+wheel or a horizontal wheel pans.
  - Ctrl/⌘+wheel zooms around the pointer column.
  - The overview strip can be dragged.
  - The default fit is `min(12, width/length)` px per column until the user zooms.
  - Letters are drawn at 8 px per column and above.
- Drawing happens in two passes (`drawBody`):
  1. Row chrome: selection tint for the drawer's selected genomes, hover tint, leaf colour, gene name and separator.
  2. MSA residue blocks inside one clip region, with Clustal-like `residueColor`.
- `computeConservation` (client side) gives each column's most common residue share among non-gap residues, plus its gap fraction. Both feed the bar and the tooltip.

**Row order** (`orderAlignmentRows`):

- With a gene tree, rows follow its leaf order, so tree tips, colour cells, names and MSA rows line up. Each leaf is matched to an MSA row by `gene`, `raw` or label.
- Without a tree, rows follow the server's MSA order (genome row, then gene).

**Gene-tree drawing:** `layoutGeneTree` computes children, depth, edge depth, leaf intervals and leaf indices from the pre-order arrays. `treeNodeX` places nodes as a phylogram or as a cladogram with tips aligned. Duplication nodes are drawn as white diamonds, with dotted leaders from each tip to the colour column. `leafColor` colours by the selected genome metadata field (missing = grey). With no field selected, it uses a stable identity colour: an FNV-1a hash of the genome name, modulo 12 palette entries.

Stepping (←/→ or the arrow buttons) calls `stepOg`. It finds the current OG's position in the current order with `positionOfOg`, which looks in `activeOgs` first and otherwise does a search for the ID. It then fetches the neighbouring position's ID and `replaceState`s the hash. On close, `restoreMatrix` restores the drawer if the view was opened from it. If the user stepped to another OG, it recentres the matrix on that OG with `jumpToOg`.

### 7.12 URL hash state and routing

`hashParams()` writes `order`, `cx` and `cy` (2 decimals), `zx` and `zy` (3 decimals), `mode`, `q` (if any) and the selection, which is `genome=<name>` or `rows=r0-r1` for a binned range. `updateHash()` adds `og=<OG>` when the drawer is open, or `og=<OG>&view=aln` in the alignment view.

| Parameter | Meaning |
| --- | --- |
| `order` | Order name or dynamic ID. An unknown ID is added to the dropdown; if the server returns 404 it shows a notice toast and falls back to `tree`. |
| `cx`, `cy`, `zx`, `zy` | The view. If present, the view is not refitted on resize (`fitted = false`). |
| `mode` | `presence` or `copies`. |
| `q` | Search query, restored and executed. |
| `genome` / `rows` | Selected genome or row range. `rows` wins (`decodeSelection`). |
| `og` | Open drawer for this OG. |
| `view=aln` | Show the full-screen alignment view for `og`. |

Most updates use `history.replaceState`, debounced by 200 ms through `scheduleHash`, so they do not create history entries. Only opening the alignment view uses `pushState`, which makes Back return to the matrix.

`applyHash()` runs once at boot and applies everything. Afterwards, `onLocationChange()` (bound to both `popstate` and `hashchange`) handles **only** the routing parameters `og`, `view` and `genome`/`rows`:

- open, step or close the alignment view;
- open the drawer, or re-select within it, for selection-only changes;
- close the drawer when `og` is gone.

It always invalidates `navSeq`. Hand-edited `order`, view or `mode` values are not re-applied until reload.

### 7.13 Keyboard handling

`onKeyDown` (on the document):

- **Esc** closes the alignment view if it is open, otherwise the drawer. It is ignored while focus is in an `INPUT`, `SELECT` or `TEXTAREA`.
- **←/→**, and also **[** and **]**, step to the previous or next OG in the current order while the drawer or the alignment view is open, again not while typing.

The search input handles its own Esc (clear) and Enter (run).

### 7.14 The mock server

`web/mock/server.ts` (`bun run mock`, port 8000, `--port N`) implements the API with a deterministic PRNG: 50 genomes in 5 clades and 2,000 OGs, plus annotations, metadata, MSAs, gene trees with duplications, and downloads. It also simulates edge cases:

- some MSAs come from `famsa` with a 250 ms delay;
- OG index 1500 returns **503** at most once every 3 s, to test the retry path;
- OGs with no MSA or tree return 404;
- a synthetic 2,000-row OG, `OG0002000`, stresses the virtualisation.

Non-API paths are served from `src/orthopan/static`, so a production build can be checked without Python. The mock follows the same limits as the server. When an endpoint changes, update the mock too.

---

## 8. Testing

The day-to-day commands and the prebuilt test datasets are described in [`docs/development.md`](development.md). This section covers what each test file checks.

### 8.1 Python (`uv run pytest`)

At commit `61d3b2b` the suite had 229 tests and took about 4 s on this machine; it is about that size now. `tests/conftest.py` defines the session fixture `test_data`, which downloads `orthologs.tsv` and `orthologs-descriptions.tsv` from demo.abrinca.com into `tests/data/` (gitignored) on first use, and skips if the download fails. No current test requests it; `scripts/build_datasets.sh hog` downloads the same files.

| File | Covers |
| --- | --- |
| `test_io.py` | Header and headerless Format B, gzip detection by magic bytes, padding, the over-wide row error with its line number, empty input, the close semantics of the iterator, and Format A genome derivation (prefix, gz, no underscore, header). |
| `test_scripts.py` | `make_test_metadata.read_genomes_and_ogs` genome/label detection for Format A, headed Format B and headerless Format B (the first line is data). |
| `test_tree.py` | UPGMA topology, gz input, branch lengths and tip-to-tip distance equal to Jaccard, iterative deep trees, quoting, a single genome, short and long rows, empty inputs. |
| `test_importer.py` | Row order, MRCA, purity and `order_tree`, frequency order, the presence bit round trip, membership sorting, sparsity and multi-batch merge, track joins, collisions and non-finite values, label validation, UPGMA fallback, the Newick parser (quotes, comments, malformed input, missing `;`), polytomy MRCA, tie breaking, brute-force MRCA on random trees, CLI. |
| `test_server.py` | Meta, genome and tree arrays, matrix full detail and aggregation against brute force (including windows wider than a byte), track means, OG lists and detail, contrast orders, order-ID collisions, LRU eviction consistency, lock boundedness, eviction during a build, re-post reuse, body limits (including chunked bodies), validation errors, the cell cap, 404s. |
| `test_tracks_basic.py` | Paralog fraction, the length CV with and without FASTA, length loading (`*` stripping, directories, ambiguous IDs), padding, headerless and gz input, `write_tracks` formatting. |
| `test_tracks_orthofinder.py` | Duplication counts and the support threshold, the zero-versus-empty rule, tree metrics with unusual labels and signed branch lengths, missing files, short rows. |
| `test_tracks_msa.py` | Pooled results equal to anvi'o (handmade and random cases, against `resources/anvio/anvio/homogeneityindex.py`), edge cases, genome weighting, weighted fragments, the trimmed tri-state, header matching with tricky genome names and fallbacks, parallel equals sequential. |
| `test_v1.py` | Prefix stripping and raw IDs, Format A import, annotation and metadata typing, name clashes and unknown IDs, pattern grouping, the pattern endpoint and order, copies against brute force, search (case, IDs, limit, nulls, order), annotation-track ties and rejection, the format-version check, `test_orthofinder_end_to_end`. |
| `test_v2.py` | MSA header mapping and row order, gene-tree leaf mapping and duplications (fake and real OrthoFinder data), downloads, FAMSA (aligns and caches, 503, unavailable, too many sequences, probe, invalid output not cached, single run per OG, timeout cleanup), the FASTA index (duplicates, gz policy, prefix, equal-length duplicates, stale `sequences.parquet` cleanup on reimport), Newick round trip, filename encoding. |

`resources/` is gitignored, and several tests depend on it. Every such test **skips with a clear reason** when the dependency is absent, so a clean clone passes:

- `test_orthofinder_end_to_end` (test_v1) and the `demo`/`demo_client` tests (test_v2) are skipped when the demo results are missing.
- `test_tracks_msa.py`'s reference-comparison tests skip (through an `anvio` fixture) when `resources/anvio/anvio/homogeneityindex.py` is missing; its synthetic tests still run.
- The `requires_famsa` tests are skipped when `scratch/.venv/lib/python3.12/site-packages/orthofinder/bin/famsa` is missing.
- `ORTHOPAN_TEST_NO_RESOURCES=1` forces all of the above to skip, which is how the skip paths are tested.
- `tests/conftest.py`'s `test_data` fixture downloads from demo.abrinca.com on first use and skips if the download fails.

Server tests use `fastapi.testclient.TestClient(create_app(dir))` on small datasets imported in `tmp_path`. FAMSA failure paths monkeypatch `Pangenome._run_famsa`.

### 8.2 Frontend (`cd web && bun test`)

There are 135 tests in 15 files. They exercise only pure modules; nothing needs a DOM except the small stubs in the tests.

| File | Covers |
| --- | --- |
| `transform.test.ts` | DPR backing size, transform round trip, bitmap ↔ label alignment, clamping and zoom limits, window and cap maths, `bufferAt`. |
| `matrix.test.ts` | Tile activation fade (thresholds, smoothstep, monotonicity), gap/radius clamps, `tileParams` (flat when binned or small, scaled by DPR), `isShaderError`. |
| `theme.test.ts` | Theme resolution (`auto` vs explicit, OS changes), `paletteFor`, the per-theme presence/copies endpoints and legends, categorical/residue palettes, and per-theme ramp luminance monotonicity. |
| `toast.test.ts` | The toast queue (hint restart, max-3 eviction, sticky errors), error keying and the network collapse, pause/resume timing with the 1.5 s floor, helpers, and `NetworkError` classification. |
| `requestkey.test.ts` | `requestKeys`/`needsFetch`: unchanged state fetches nothing, order/track keys change independently, and a cleared displayed state forces a refetch. |
| `tracks.test.ts` | Strip hit testing (including stale spans), `binRange`. |
| `v1.test.ts` | LE decoders, copies maths and saturation, the categorical palette, `hitViewport`. |
| `pending.test.ts` | `PendingCache` sharing, `forget`, bounds. |
| `search.test.ts` | `SearchBox` staleness on query and order changes. |
| `generation.test.ts` | `Generation` abort and current checks, `Sequence`. |
| `selection.test.ts`, `members.test.ts` | Selection summaries and hash encoding; member slicing with pinning. |
| `alignment.test.ts`, `msa.test.ts`, `genetree.test.ts` | Row order, visible band, stepping; column maths, conservation, colours; tree layout, `treeNodeX`, leaf colours. |

`bun run build` runs `tsc --noEmit` (strict, with no unused locals or parameters) before `vite build`, so type errors fail the build.

### 8.3 Browser QA

Rendering, alignment and interaction are checked in a real browser with **playwright-core** driving the system **`/usr/bin/chromium-browser`** from bun scripts. The Playwright MCP cannot launch Chrome in this environment. Most scripts call `chromium.launch({executablePath: '/usr/bin/chromium-browser', headless: true, args: [... '--use-angle=swiftshader', '--force-device-scale-factor=<dpr>']})` for software WebGL. `lib.ts` instead connects over CDP to an already running Chromium on port 9444.

The scripts live in `scratch/vision-qa/` (`v1-check.ts`, `v2-check.ts`, `dpr-probe.ts`, `stale-check.ts`, `hash-check.ts`, `large-check.ts`, `scale10k.ts`, …), have their own `node_modules` (`playwright-core`, `pngjs`), and are **not committed** (`scratch/` is gitignored). They:

- check pixels in screenshots (for example, that heatmap cells and labels line up at DPR 1, 1.25 and 2);
- check the tooltip content against the expected genome and OG;
- check 404 and 503 states and the absence of console errors.

A typical run:

```sh
(cd web && bun run build)                       # build the UI
uv run orthopan serve scratch/datasets/demo     # or: (cd web && bun run mock)
(cd scratch/vision-qa && bun run v2-check.ts http://localhost:8000/)
```

### 8.4 Running everything

```sh
uv run pytest                      # Python suite
(cd web && bun test)               # frontend unit tests
(cd web && bun run build)          # type check + bundle into src/orthopan/static
uv build                           # wheel including the built UI
```

---

## 9. Performance and scale

### 9.1 Measured numbers

The reference measurements are in [design-decisions](design-decisions.md) §Scale. They use a 10,000 × 50,000 dataset simulated with `scripts/simulate_pangenome.py --genomes 10000 --orthogroups 50000 --seed 42` (25.1M genes, 24.4M present cells). The dataset is `scratch/datasets/sim10k` (input in `scratch/datasets/sim10k-input/`), rebuilt by `scripts/build_datasets.sh sim10k`.

The import with `--tree`, two annotations and metadata took 52.2 s wall time, with a peak RSS of 963 MB. These were re-checked on 2026-10-05 against `scratch/datasets/sim10k`: the server measurements below, plus `/api/orthogroup/{og}` at about 23 ms for a core OG (10k members) and about 5 ms for a rare one, and `/api/search` at 10–36 ms.

| Measurement | Result |
| --- | --- |
| Server ready (process start → `/api/meta` answers) | 0.75 s |
| Steady RSS after prewarm (tree + frequency packed) | ≈ 0.47 GB |
| `/api/matrix`, full view, 2000 × 1000, presence | 73–88 ms |
| same, `value=copies` | 131 ms |
| `frequency` order (prewarmed) | 104 ms |
| `track:copy_number_mean`, cold / warm | 239 / 105 ms |
| `/api/tracks`, 2000 bins | 1.7 ms |
| `/api/genomes` (cached bytes, 1.26 MB) | 1.5 ms |
| `/api/tree` (0.54 MB JSON) | 35 ms |
| Peak RSS (VmHWM) across these requests | ≈ 1.03 GB |

### 9.2 Memory model at 10k genomes × 50k OGs

Let `B = n_genomes × ceil(n_og/8) ≈ 62.5 MB`, the size of one packed matrix.

**Importer**

| Item | Size | Lifetime |
| --- | --- | --- |
| Presence bytearrays | ≈ B (grown in 4 KiB steps per genome) | whole import |
| Membership buffer | ≤ 10k rows of gene tuples | per batch |
| Per-OG lists (stats, tracks, MRCA) | O(n_og) Python objects | whole import |
| `--fasta` index dict + `seq_rows` | O(total FASTA records) Python tuples | whole import (large at scale) |
| Concatenated `raw` + `packed` | ≈ B | end |
| `_pattern_ids` dense unpack | **n_genomes × n_og bytes ≈ 500 MB** | end; the dominant peak |
| UPGMA (no `--tree`): bool matrix, packed copy, `inter` and `union` (n × n float32, 400 MB each), condensed distances, scipy linkage | ≈ 1 GB, per [design-decisions](design-decisions.md) | before the import loop |

**Server**

| Item | Size |
| --- | --- |
| `presence_tree` | B |
| Each cached packed order | B (+ one byte per row) |
| `_packed_cache` (≤ 8) + `_dynamic` (≤ 8) | worst case 16 × B ≈ 1 GB |
| Multicopy arrays (731k entries × 3 int64) | ≈ 17.5 MB |
| Orthogroups/annotations frames, genome JSON, tree dict | tens of MB |
| `_seq_lookup` (with `--fasta`) | one Python dict entry per indexed gene, which is large at 25M genes |
| Per full-window matrix request: popcount (uint8), cumsum (int32 rows × span) and gathers, across ≤ 16 threads | several hundred MB, transient (seen as the 0.47 → ≈ 1 GB HWM) |
| `contrast_perm` | `nsel × n_og` bytes, transient (up to ≈ 500 MB for a near-total selection) |

**Browser:** about 1.3 MB of genome JSON, the tree arrays, one `ImageData` of at most 4M cells × 4 bytes = 16 MB plus its texture, member JSON of about 1 MB for a 10k-member OG (rendered 400 genes at a time), and up to 8 cached MSAs.

### 9.3 Known hotspots

- **Import:** `_pattern_ids` (dense unpack), the pure-Python per-cell loop, and the UPGMA `n²` matrices with scipy linkage when there is no tree.
- **First use of an order:** `_build_packed` (about 230 ms at 10k). Track orders also re-sort on every request, because `_track_perm` has no cache.
- **Full-window matrix requests:** the cumsum dominates (75–130 ms). Copies mode adds an O(n_og) inverse permutation and a filter over all multicopy entries.
- **`annotation_track`:** a Python loop over `w` bins. Fine at w ≤ 4096.
- **`search`:** builds a dict for every (OG, field) hit before it slices to `limit`, so very broad queries cost O(hits).
- **`/api/orthogroup/{og}`:** a linear polars filter by ID, plus the Parquet scan.
- **`/api/tree`:** it is not pre-serialised like `/api/genomes` and takes 35 ms at 10k.
- **Startup with `--fasta` at scale:** loading all of `sequences.parquet` into a dict.

---

## 10. Extending OrthoPan

### 10.1 A new track preprocessor

1. Add `src/orthopan/tracks/<name>.py` with a `COLUMNS` list, a `build_tracks(...)` function that returns `{og: {column: float | None}}`, and a `main()` that calls `write_tracks(rows, COLUMNS, out)`. Read TSVs with `orthopan.io.read_matrix` and `open_text` so gzip and both formats work.
2. Choose column names that are not in `importer._ORTHOGROUP_COLUMNS`, and not used by other preprocessors that may run in the same import. Otherwise the import fails with `duplicate track column name`.
3. Users can run it themselves and pass `--track out.tsv`. To run it automatically for OrthoFinder input, append `_dict_source(COLUMNS, build_tracks(results_dir))` in `importer._orthofinder_track_sources()`.
4. Nothing else is needed. The server lists every `manifest.track_columns` entry in `/api/meta`, serves it from `/api/tracks`, offers `track:<name>[:desc]` orders and shows it in the OG detail. The frontend lists it in the Layers menu and the order dropdown. To show it by default, change `pickTracks()` in `main.ts`.
5. Add `tests/test_tracks_<name>.py`. No `format_version` bump is needed.

### 10.2 A new annotation type

Annotations have two kinds today, `text` and `categorical`, and are stored as `Utf8`. A new kind, for example a numeric annotation, touches these places:

1. **Importer:**
   - choose the kind in `_annotation_specs()` and add any spec fields;
   - set the dtype for it in the `annotations.parquet` writer in `run()`.
2. **Server:**
   - `Pangenome.__init__`: precompute codes like `_ann_codes`, and decide whether the column joins `_search_fields`. `_field_mask` uses string operations, so exclude non-string columns.
   - Add a binning method next to `annotation_track()` and its route.
   - `_annotation_values()` converts every value with `str()`; adjust it if needed.
3. **API docs:** [docs/api.md](api.md).
4. **Frontend:**
   - widen `AnnotationInfo.kind` in `api.ts`;
   - add a fetch function;
   - extend `populateLayers` (it currently lists only categorical annotations), `buildStrips` and the `Strip` union in `tracks.ts`, with a new `draw*` method and hit text in `onStripHover`.
5. **Mock server.**
6. **Tests:** `test_v1.py`-style importer and server tests, plus bun tests for the drawing maths.

Bump `format_version` only if the server can no longer read older directories correctly (see [§10.6](#106-when-to-bump-format_version)).

### 10.3 A new API endpoint plus its frontend client

1. Put the logic in a `Pangenome` method. Keep `create_app` routes thin: parse, validate against the module constants or new ones, call the method, build the response.
   - Binary arrays: `Response(content=arr.astype("<f4"/"<u2"/"<u1").tobytes(), media_type="application/octet-stream", headers={...})`.
   - Large JSON: `json.dumps(..., separators=(",", ":"))` in a `Response`, as `/api/orthogroup/{og}` does.
   - Raise `HTTPException` with a `detail` string.
   - If the method touches shared caches, take `self._lock` only for bookkeeping and use per-key locks for builds, as `_packed` and `_famsa_lock` do.
2. Document the endpoint in [docs/api.md](api.md).
3. Add a typed function to `web/src/api.ts`:
   - JSON: `getJSON` (or `getJSONDetail` to surface the `detail` text).
   - Binary: `fetch` + `res.ok` check + **byte-length check** + `decode.ts` LE decoder.
   - Always accept an `AbortSignal`.
4. Call it under the right guard: `gen`/`inflight` for view-dependent data, `drawerGen` or `alnGen` for OG detail, and a cache with a bound (§7.4). Re-check staleness after every `await`.
5. Add the route to `web/mock/server.ts`.
6. Add tests: `TestClient` in `tests/test_server.py` (or a v1/v2 file), and bun tests for any pure logic you add.

### 10.4 A new order type

**Dynamic (computed on request):**

1. Add the ID prefix to `_DYNAMIC_PREFIX`.
2. Add a branch in `create_order`. Validate the body, **canonicalise** it, and hash it into a deterministic ID such as `"<type>:" + sha1(...)[:12]`.
3. Add a `<type>_perm()` method that returns an int64 permutation of length `n_og` (position → `og_idx`), and call it in the compute branch.
4. The LRU, locking, lazy packing and 404-on-eviction then work as for the existing dynamic types.
5. Frontend: `postOrder({by: '<type>', ...})`, then `setOrder(id)`. Add it to the mock.

**Built-in (precomputed at import):**

1. Add an `order_<x>` rank column with `_ranks()` in `run()`, and add it to `_ORTHOGROUP_COLUMNS` and the schema.
2. In `Pangenome`, add `perm_<x>`, a branch in `order_perm` and the name in `meta()["orders"]`. Optionally prewarm it.
3. Bump `format_version`, because the server now requires the column.

### 10.5 A new detail-view column (alignment view)

1. Add the width constant in `alignment.ts` and include it in `ALN_MSA_X`, which is the left edge of the MSA. Hit testing (`colAt`, `mousemove`, `barMousemove`, `barPointerdown`) and the bar header use `ALN_MSA_X`, so they follow automatically.
2. Carry the data in `AlnRow` and/or `AlnContent` and fill it in `orderAlignmentRows` or in `loadAlignment`'s `build()` in `main.ts`.
3. Draw it in `drawBody` pass 1, using the row loop's `y` and `ALN_ROW_H`. Add a divider in `drawGuides` and a header label in `drawBar`.
4. If the data must come from the server, extend `/msa` or `/gene-tree` (§10.3) and the mock.
5. Add a case to `web/test/alignment.test.ts`.

### 10.6 When to bump `format_version`

`FORMAT_VERSION` in `importer.py` is written to the manifest. `Pangenome.__init__` rejects anything older than the minimum it supports, currently `< 2`.

- **Bump it, and raise the server's minimum,** when a server running the new code would misread or crash on a directory written by older code. Examples: a new **required** file or column, a changed bit layout or row orientation, or a changed meaning of an existing field (`multicopy.row` switching from tree row to `genome_idx` would be one).
- **Don't bump it** for additions that the server treats as optional when absent: new track columns, new optional files such as `sequences.parquet`, new manifest keys read with `.get(...)`. Version 2 already covers the v2 detail-view additions this way.

When you bump, update the check message ("re-run orthopan import") and `test_outdated_format_version_fails`.

---

## 11. Known limitations and gotchas

### 11.1 Input and import quirks

- **`tracks.basic` emits only `paralog_fraction` and `length_cv`.** Copy number is a built-in import column (`copy_number_mean`/`copy_number_max`), so a `basic.tsv` can be passed straight to `--track`; if you add a copy-number column back it will collide and raise `duplicate track column name`.
- `--fasta` takes `nargs="+"`. In `orthopan import --fasta a.faa orthologs.tsv -o out` the orthologs path is swallowed by `--fasta`. Put `--fasta` after the positional argument or after other options.
- Track files must contain numbers. A non-numeric cell aborts the import with a bare `ValueError`.
- `orthopan-tree` has no `--strip-prefix` option, so it always uses the default prefix. It accepts both Format A and Format B inputs.
- Without `--orthogroups`, `python -m orthopan.tracks.msa DIR` silently falls back to pooled (anvi'o) weighting, with only a warning. Without `--seqs`, `trimmed` is empty.
- Manifest paths for `detail` and `fasta` are absolute. Moving or deleting the OrthoFinder results or the FASTA files breaks the detail views, but moving the dataset directory itself is fine. Re-importing into the same directory does not clear `cache/msa/`, so stale FAMSA alignments from a different run can be served for reused OG IDs. Delete `cache/` when you re-import.
- Without `resources/` (gitignored), the tests that depend on it **skip with a reason** instead of erroring. `ORTHOPAN_TEST_NO_RESOURCES=1` forces those skips (see §8.1).

### 11.2 OrthoFinder quirks

- **Missing MSAs:** in OrthoFinder 3.1.5, about 2% of OGs with 2–3 genes have no MSA. The cause is unknown ([design-decisions](design-decisions.md) §Per-orthogroup tracks). Their MSA tracks are empty. The viewer aligns them with FAMSA on demand if it can, and otherwise returns 404 with a reason.
- **Trimming:** OrthoFinder trims MSAs by default (unless run with `-z`). Trimming inflates geometric homogeneity. The `trimmed` track flags it, the `/msa` response carries `trimmed`, and the alignment view shows a "trimmed" badge.
- **Gene trees** exist only for OGs with at least 4 genes. Otherwise the endpoint returns 404 and the frontend caches it as `missing`.
- **Unstable OG IDs:** IDs are renumbered between OrthoFinder runs, including `--assign` runs ([orthofinder-outputs §2.2](orthofinder-outputs.md)). Treat each results folder and dataset directory as a sealed snapshot. Bookmarked URLs (`og=…`, `pattern:` orders) are only valid for the same dataset.
- **The deleted `N0.tsv`:** OrthoFinder 3.1.5 renumbers the N0 HOGs into `Orthogroups.tsv` and deletes `Phylogenetic_Hierarchical_Orthogroups/N0.tsv`, so `Orthogroups.tsv` is the matrix unit. The demo results contain `N1.tsv`–`N22.tsv` but no `N0.tsv`. HOG files can still be imported as Format A or B, but there are no hierarchy features.
- **MSA headers** are `<species with . → _>_<gene ID>`. Mapping them back to genomes must be exact (`_expected_headers`), because prefix matching is ambiguous for names like `FAM1` and `FAM10`, or `A` with gene `B_g` versus `A_B` with gene `g`. Ambiguous headers lead to pooled weighting in the tracks, and to `row = −1` in the viewer.

### 11.3 Python 3.14 and the polars-on-3.15 history

The project began on Python 3.15 (`requires-python = ">=3.15"`, `.python-version` set to 3.15). On that interpreter the installed polars build silently returned `None` from Series methods such as `Series.str.*`, `fill_null`, `is_null` and `list.len`, because of a broken bytecode check in polars' expression dispatch. Only DataFrame **expressions** worked. A guard test (`tests/test_polars_guard.py`) scanned `src/` for such Series calls.

Commit `84d6d75` moved the project to Python 3.14 (`requires-python = "==3.14.*"`, `.python-version` = 3.14) and removed the guard test. The code still uses the expression style, for example `Pangenome._field_mask`, which works on both versions. Keep using expressions and don't widen `requires-python` until polars works on 3.15 again.

Separately, OrthoFinder and FAMSA live in their own Python 3.12 environment, `scratch/.venv`. Leave it alone (see `AGENTS.md`).

### 11.4 FAMSA is an external binary

- It is found with `shutil.which(--famsa)` at startup and probed by running a tiny alignment. If it is missing or broken, `meta.detail.famsa` is false and MSA requests for OGs without an alignment return 404. The probe runs on every server start.
- OGs with more than 2,000 sequences are never aligned. Alignments run single-threaded (`-t 1`, as in OrthoFinder) and two at a time; a third concurrent request gets 503 and the UI shows **retry**. An alignment that takes more than 120 s returns 500.
- Alignments are cached in `<dataset>/cache/msa/`, which must be writable. The cache has no invalidation (see §11.1).
- `_famsa_locks` gains one lock per OG that is ever aligned. It is not pruned.

### 11.5 Canvas size limits

Browsers cap canvas and texture sizes. These limits shaped several choices:

- The matrix request is capped at 4,096 cells per axis and 4M cells in total, matching the server and typical WebGL texture limits. The bitmap is always one viewport-sized window and never the whole matrix.
- Canvas2D panels are sized to their CSS boxes × DPR, so the backing stores scale with the window, not with the data.
- The alignment view uses one sticky, viewport-tall canvas and draws only the visible band and columns. A tall alignment (16 px × 2,000 rows = 32,000 px) would otherwise exceed common canvas height limits, which are around 32k px and lower in some browsers.

### 11.6 Localhost-only server

`orthopan serve` binds to `127.0.0.1` by default and has no authentication, CORS configuration, TLS or global request limiter ([api.md](api.md) §Notes). The per-request caps in §5.13 bound the cost of any single request, but not the number of concurrent requests. Requests can trigger heavy work: a 64 KiB contrast order costs up to n_genomes × n_og bytes, and an MSA request can start FAMSA. `--host 0.0.0.0` exposes all of this, including FASTA and alignment downloads, to the network. To share a viewer, put it behind an authenticating reverse proxy.

### 11.7 Other behaviours worth knowing

- `tree`, `frequency` and `track:*` packed matrices share one 8-entry LRU, so heavy track sorting can evict `tree`. It is rebuilt transparently, which costs a cold-order delay.
- Contrast sorts break ties by `og_idx`, while search and pattern orders keep tree order within groups.
- `positionOfOg` and `jumpToOg` locate an OG that is outside the loaded label window by **substring search** with `limit=5`, taking the exact match if it is among those hits and otherwise the first hit. IDs that are substrings of many other IDs or annotations can therefore make stepping or "jump" land on the wrong OG.
- Members are listed in `genome_idx` (input) order, not tree order. The selected genome is shown separately at the top.
- Downloaded Newick branch lengths are rewritten with six significant digits (`:g`).
