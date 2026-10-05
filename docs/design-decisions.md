# OrthoPan design decisions

Decisions taken so far (2026-10-04). See `VISION.md` for context and `docs/anvio-internals.md` for the anvi'o analysis behind them.

## Scope

- **MVP is the viewer first**: presence/absence heatmap plus genome tree, built from an orthologs TSV alone. Homogeneity tracks follow later.
- **Delivery**: local server (`orthopan serve`). A Python backend serves slices of the matrix to a browser UI.
- **Frontend**: Deck.gl.

## Inputs

- The orthologs TSV is required. Both layouts in `docs/ortholog_file_formats_documentation.md` count as flat "group -> genes" mappings, whether they come from `Orthogroups.tsv` or HOG `N*.tsv`. There are no HOG-hierarchy features.
- The MSA directory is optional. Homogeneity tracks appear only when MSAs are present. Warn if the alignments look gap-trimmed, because OrthoFinder trims by default unless run with `-z`, and trimming inflates geometric homogeneity.
- The genome tree comes from OrthoFinder's `SpeciesTree_rooted.txt` if available. Otherwise `orthopan-tree` builds it with UPGMA on Jaccard distances between presence/absence profiles.

## Scale

- Reference dataset: 2,193 genomes × 50,515 orthogroups. The long-term target is 10k+ genomes.
- Measured on a simulated 10,000 × 50,000 dataset (`scripts/simulate_pangenome.py`, 25M genes):
  - import with `--tree`: about 52 s at about 1 GB peak memory;
  - server: 0.8 s startup, about 0.6–1 GB memory;
  - full-screen matrix: about 80 ms for presence and about 100 ms for copies;
  - first use of an order: ≤230 ms.
  Without `--tree`, the UPGMA fallback needs about 1 GB.
- The **data layer is designed for 10k+ from the start**: sparse Parquet, bit-packed presence/absence, and a backend that serves windowed slices and zoomed-out aggregates.
- The **renderer starts simple**: it draws only the visible window and requests aggregates from the backend when zoomed out. Server-side tiling can be added later behind the same API.

## Storage and rendering

- The processed pangenome is a directory of Parquet files queried with DuckDB or polars. There is no server database. Planned files:
  - `orthogroups.parquet`: per-orthogroup keys and scores.
  - `genomes.parquet`: one row per genome.
  - `membership.parquet`: sparse membership (orthogroup, genome, copies, gene IDs).
  - A bit-packed presence matrix.
- When many cells share one pixel, show the fraction present.
- Cells encode presence by default. Copy number is an optional color mode.
- Backend: FastAPI. It also serves the built frontend.
- Frontend: Vite + TypeScript + deck.gl, without React.
- Matrix windows travel as raw binary arrays. The client applies colors. Revisit this if something better turns up.
- Zoomed-out aggregates are computed on the fly with a numpy block mean over `presence.npy`. Add caching or an import-time pyramid only if 10k genomes turn out too slow.
- The API contract is in `docs/api.md`.
- The first build step is `orthopan import`: TSV plus optional tree in, Parquet directory with orderings out.

## Annotations, metadata and detail data (v1)

- `--annotation FILE` can be repeated.
  - A 2-column file without a header is a free-text column named after the file stem (`description` if the stem is `*descriptions*`).
  - A file with a header (first column = orthogroup ID) gives named columns. A string column with at most 64 distinct values is categorical; otherwise it is text.
  - Text and categorical annotations appear in the tooltip and panel, and are searchable. Categorical annotations are also drawn as a colour strip.
- `--genome-metadata TSV` has a header, and its first column is the genome name. Each field is drawn as a narrow coloured column between the tree and the heatmap, using a categorical palette or a numeric gradient.
- Gene ID prefixes (default regex `^gnl\|[^|]+\|`, configurable) are stripped at import for display. The raw ID is kept, because it's needed to match OrthoFinder MSA and FASTA headers.
- `--orthofinder RESULTS_DIR` auto-detects `Orthogroups.tsv` and the species tree, runs the track preprocessors, and records the MSA, gene-tree and sequence paths in the manifest. The server reads those files on demand for the single-orthogroup view.
- Format A (2-column orthogroup → genes, with the genome taken from the gene-ID prefix) is supported, so `tests/data` can be used.
- Each orthogroup gets a `pattern_id` (identical presence/absence rows share one) for "same pattern" queries.
- Copy-number colour mode: a zoomed-out bin shows the mean copy number over its present cells. Absent cells stay white.
- Order of work: annotations, metadata, search and copy number first; then the single-orthogroup view (MSA viewer, gene tree, members table).

## Single-orthogroup view (v2)

- The detail drawer (Members tab) shows the orthogroup's members grouped by genome. The clicked genome, or the binned row range, appears at the top as "In genome X: N genes" or "absent in X". Lists with more than 400 genes are paginated.
- "Open alignment view" opens a full-screen view (`#og=…&view=aln`). Columns are gene tree | metadata colours | gene name | MSA, one 16 px row per gene, scrolled with the page. Only the visible band is drawn, so the canvases stay under browser size limits. Rows follow the gene tree's leaf order, or genome tree order when there is no gene tree.
- MSA: our own canvas viewer with no dependencies. It has residue colours, a conservation/gap bar and an overview strip, and only the MSA scrolls horizontally.
- Gene tree: leaves are coloured by a selectable genome metadata field, with a legend. Duplication nodes come from `Duplications.tsv`. The view switches between phylogram and cladogram.
- ←/→ step to the previous or next orthogroup in the current matrix order, keeping the selected genome. Esc and Back return to the matrix.
- Missing MSAs are aligned on demand with the same FAMSA command OrthoFinder uses, and the results are cached. Non-OrthoFinder inputs supply sequences via `--fasta`, indexed by byte offset.
- Downloads: unaligned FASTA, aligned FASTA and the Newick gene tree.

## Per-orthogroup tracks

- Tracks are the extension point. A track file is a TSV with the header `orthogroup<TAB><track>[<TAB><track>...]`, float values, and an empty cell for missing values. Users can supply their own track files as custom tracks.
- Copy number is a default track, computed at import from the orthologs TSV.
- Preprocessors live in `src/orthopan/tracks/` and write track TSVs (`python -m orthopan.tracks.<name>`):
  - `basic`, for any tool: copy number mean/max, paralog fraction, and protein length CV (needs FASTA).
  - `orthofinder`: duplication counts and gene-tree length per leaf.
  - `msa`: anvi'o-style geometric and functional homogeneity, gap fraction, and a trimmed flag.
- Rationale (`docs/orthofinder-outputs.md`): cluster quality has three nearly independent axes, namely structure, divergence and paralogy. Length CV and copy number cover structure and paralogy without alignments.
- Orthogroup IDs are not stable across OrthoFinder runs, so each results folder is a sealed snapshot. `Orthogroups.tsv` is the matrix unit, because OrthoFinder 3.1.5 already renumbers N0 HOGs into it and deletes `N0.tsv`.
- In OrthoFinder 3.1.5, about 2% of orthogroups with 2–3 genes have no MSA. The cause is unresolved. MSAs are optional per orthogroup.

## Homogeneity (MSA track)

- Compute geometric homogeneity in closed form from per-gene and per-column gap counts. It is numerically identical to anvi'o's pairwise XOR formulation.
- **Paralogs: weight per genome.** Each genome contributes equal total weight and all copies are kept. This differs from anvi'o, which pools all genes. An anvi'o-exact mode is useful for validation against fixtures.

## Orthogroup ordering

- **Default: tree-guided.** Each orthogroup is placed at the DFS position of the MRCA of the genomes containing it. Ties break by clade purity, then frequency. This produces clade-specific blocks above their clades.
- **Alternate: frequency gradient** (core to accessory).
- Both are precomputed at import as permutation columns. Per-orthogroup sort keys (frequency, MRCA node, purity, copy number, homogeneity, annotation) are stored as well.
- On-demand sorts by any track and clade-contrast sorts ("present in selection, absent elsewhere") run on the server.
- There is no full hierarchical clustering of orthogroups. 91% of orthogroups are rare accessory, and there are 26k unique patterns, so it costs too much memory for too little benefit.
