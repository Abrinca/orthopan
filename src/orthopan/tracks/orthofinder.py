"""Per-orthogroup tracks from an OrthoFinder results directory.

Usage::

    uv run python -m orthopan.tracks.orthofinder RESULTS_DIR [-o OUT]

Only top-level result files are read; ``WorkingDirectory/`` is never touched.

Inputs (all paths relative to ``RESULTS_DIR``):

- ``Orthogroups/Orthogroups.tsv`` defines the OG universe and row order.
- ``Gene_Duplication_Events/Duplications.tsv`` is the per-event duplication
  table (columns ``Orthogroup, Species Tree Node, Gene Tree Node, Support,
  Type, Genes 1, Genes 2``).
- ``Resolved_Gene_Trees/Resolved_Gene_Trees.txt`` holds one resolved Newick
  per OG, ``OGxxxxxxx: <newick>;``.

Missing optional files simply leave their columns empty.

Duplications for OGs absent from ``Duplications.tsv``: OrthoFinder only writes
an event row for an OG that has a gene tree, so such an OG is a true **zero**
when it appears in ``Resolved_Gene_Trees.txt``, and **empty** when it has no
tree (unknown, not zero). If the tree file is missing entirely, every OG
absent from ``Duplications.tsv`` is left empty.

``duplications_supported`` counts events with ``Support >= 0.5``, OrthoFinder's
own threshold. Tree metrics use a minimal parse of the Newick: branch lengths
are the ``:number`` tokens (scientific notation included) and the leaf count is
``commas + 1``, which holds for any rooted Newick whose labels contain no
comma. ``tree_length_per_leaf`` is a cheap divergence proxy.
"""

import argparse
import re
import sys
from pathlib import Path

from orthopan.io import open_text, read_matrix
from orthopan.tracks import write_tracks

COLUMNS = ["duplications", "duplications_supported", "tree_length", "tree_length_per_leaf"]

SUPPORT_THRESHOLD = 0.5

# Signed decimal/scientific branch lengths, including ".5" and "-0.1".
_BRANCH = re.compile(r":(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)")


def read_orthogroups(path):
    """Return the OG IDs, in file order, from Orthogroups.tsv."""
    _, rows = read_matrix(path)
    return [og for og, _ in rows]


def read_duplications(path):
    """Return ({og: events}, {og: events with Support >= 0.5}, skipped row count)."""
    labels, rows = read_matrix(path)
    support = labels.index("Support") if "Support" in labels else 2
    counts, supported, skipped = {}, {}, 0
    for og, cells in rows:
        if not og or len(cells) <= support:
            skipped += 1
            continue
        try:
            value = float(cells[support])
        except ValueError:  # empty Support cell from a short, padded row
            skipped += 1
            continue
        counts[og] = counts.get(og, 0) + 1
        if value >= SUPPORT_THRESHOLD:
            supported[og] = supported.get(og, 0) + 1
    return counts, supported, skipped


def read_gene_trees(path):
    """Return {og: (tree length, length per leaf)} for Resolved_Gene_Trees.txt."""
    trees = {}
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            og, _, newick = line.partition(":")
            newick = newick.strip()
            leaves = newick.count(",") + 1
            length = sum(float(m) for m in _BRANCH.findall(newick))
            trees[og.strip()] = (length, length / leaves if leaves else None)
    return trees


def build_tracks(results_dir):
    """Return {orthogroup: {track: value}} for every OG in Orthogroups.tsv."""
    results_dir = Path(results_dir)
    universe = read_orthogroups(results_dir / "Orthogroups" / "Orthogroups.tsv")

    dup_path = results_dir / "Gene_Duplication_Events" / "Duplications.tsv"
    tree_path = results_dir / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt"
    if dup_path.is_file():
        counts, supported, skipped = read_duplications(dup_path)
        if skipped:
            print(f"warning: skipped {skipped} incomplete rows in {dup_path}", file=sys.stderr)
        dup = (counts, supported)
    else:
        dup = None
    trees = read_gene_trees(tree_path) if tree_path.is_file() else {}

    rows = {}
    for og in universe:
        row = {}
        if dup is not None:
            if og in dup[0]:
                row["duplications"] = float(dup[0][og])
                row["duplications_supported"] = float(dup[1].get(og, 0))
            elif og in trees:
                row["duplications"] = 0.0
                row["duplications_supported"] = 0.0
        if og in trees:
            length, per_leaf = trees[og]
            row["tree_length"] = length
            if per_leaf is not None:
                row["tree_length_per_leaf"] = per_leaf
        rows[og] = row
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", help="OrthoFinder results directory")
    parser.add_argument("-o", "--out", help="output TSV (default: stdout)")
    args = parser.parse_args(argv)
    write_tracks(build_tracks(args.results_dir), COLUMNS, args.out)


if __name__ == "__main__":
    main()
