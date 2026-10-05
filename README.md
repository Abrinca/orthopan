# OrthoPan

OrthoPan is a fast, interactive viewer for pangenomes. It imports orthogroups (for example an
OrthoFinder results directory), builds a processed dataset (genome tree, presence/absence and
copy-number matrix, per-orthogroup tracks) and serves a browser UI. The UI is a heatmap with
genomes as rows and orthogroups as columns. It has search, annotations, genome metadata, and
per-orthogroup alignments and gene trees.

## Install

Requires Python 3.14 and [uv](https://docs.astral.sh/uv/).

```sh
git clone <repo> && cd orthopan
uv sync
cd web && bun install && bun run build && cd ..   # builds the UI into src/orthopan/static
```

[Bun](https://bun.sh) is only needed to build the UI. A wheel built with `uv build` already
contains it, so installing the wheel needs no Bun.

## Quickstart

```sh
uv run orthopan import --orthofinder path/to/OrthoFinder/Results_XXX -o out
uv run orthopan serve out
# open http://127.0.0.1:8000
```

## Using the viewer

- **Matrix:** genomes are rows (in genome tree order) and orthogroups are columns. The mouse wheel zooms orthogroups, Shift+wheel zooms genomes, and dragging pans.
- **Order:** sorts the orthogroups by tree position (default), frequency or any track. Select a clade in the tree and choose "sort by contrast" to bring that clade's genes to the front.
- **Click a cell:** the drawer shows the orthogroup's members, with the clicked genome's genes at the top. **Open alignment view** shows the gene tree, gene names and MSA side by side.
- **Keys:** ←/→ step to the previous or next orthogroup; Esc closes the current view.
- **Other controls:** search box, presence/copy-number toggle, layers menu for tracks and annotations, and downloads (FASTA, aligned FASTA, Newick) for each orthogroup.

The URL keeps the current view, so you can bookmark or share it.

## Other inputs

`orthopan import` also takes a plain orthologs TSV. This can be either a per-genome matrix or a
two-column `orthogroup<TAB>genes` file. Optional flags:

| Flag | Meaning |
| --- | --- |
| `--tree NEWICK` | genome tree (default: OrthoFinder's tree, or UPGMA on gene content) |
| `--annotation TSV` | per-orthogroup descriptions or categories (repeatable) |
| `--genome-metadata TSV` | per-genome metadata, first column = genome name |
| `--track TSV` | per-orthogroup numeric track (repeatable) |
| `--fasta PATH` | protein FASTA files or directories, used for sequences and alignments |
| `--strip-prefix REGEX` | prefix removed from gene IDs for display (default `^gnl\|[^|]+\|`) |

Example:

```sh
uv run orthopan import orthologs.tsv -o out --tree tree.nwk \
    --annotation descriptions.tsv --genome-metadata metadata.tsv
```

Track files can be made with the preprocessors, and a genome tree with `orthopan-tree`:

```sh
uv run python -m orthopan.tracks.basic orthologs.tsv -o basic.tsv
uv run python -m orthopan.tracks.orthofinder RESULTS_DIR -o orthofinder.tsv
uv run python -m orthopan.tracks.msa RESULTS_DIR/MultipleSequenceAlignments -o msa.tsv
uv run orthopan-tree orthologs.tsv -o tree.nwk
```

Orthogroups without an alignment are aligned on demand if
[FAMSA](https://github.com/refresh-bio/FAMSA) is available: `orthopan serve out --famsa PATH`.

## Development

```sh
uv run pytest
cd web && bun test
bun run dev     # UI dev server, proxies /api to 127.0.0.1:8000
bun run mock    # mock API on :8000, no backend needed
uv build        # wheel including the built UI
```

To get a large test dataset, run `uv run python scripts/simulate_pangenome.py --help`. It writes
a seeded, tree-structured pangenome (for example 10,000 genomes × 50,000 orthogroups).
