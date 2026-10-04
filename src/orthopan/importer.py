"""`orthopan import`: orthologs TSV (+ tree, tracks) -> processed pangenome dir.

This is the first pipeline step (docs/design-decisions.md). It reads a Format B
orthologs matrix with :func:`orthopan.io.read_matrix`, optionally validates it
against a genome Newick tree (or builds a UPGMA tree with
:func:`orthopan.tree.upgma_newick`), joins per-orthogroup track TSVs, and writes
the files a future Deck.gl/DuckDB viewer queries:

- ``genomes.parquet``      one row per genome, tree DFS row order
- ``orthogroups.parquet``  one row per orthogroup, orderings + tracks
- ``membership.parquet``   sparse present cells (og, genome, copies, genes)
- ``presence.npy``         bit-packed presence, rows genome_idx, cols og_idx
- ``tree.nwk``             the tree used, labels as given
- ``manifest.json``        format version, inputs, counts, track columns
"""

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl

from orthopan.io import open_text, read_matrix, split_genes
from orthopan.tree import upgma_newick

FORMAT_VERSION = 1
DEFAULT_TRACKS = ("copy_number_mean", "copy_number_max")
# Every column the importer writes to orthogroups.parquet; track columns may
# not collide with any of them.
_ORTHOGROUP_COLUMNS = frozenset({
    "og_idx", "orthogroup", "n_genomes", "n_genes",
    "mrca_node", "mrca_label", "purity", "order_tree", "order_frequency",
    *DEFAULT_TRACKS,
})
_MEMBERSHIP_BATCH = 10_000  # rows per membership Parquet part
_MEMBERSHIP_ROW_GROUP = 10_000  # rows per row group in the final Parquet


class _Node:
    __slots__ = ("label", "branch", "children")

    def __init__(self):
        self.label = None
        self.branch = None
        self.children = []


def _tokens(text):
    """Yield ``(kind, value)`` Newick tokens, kind ``"punct"`` or ``"text"``.

    Comments (``[...]``, nestable) outside quotes are skipped, so a label or
    branch length is never contaminated by them. A quoted label stays a single
    text token even when it contains punctuation such as ``( ) , : ;``.
    """
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif c in "(),;:":
            yield "punct", c
            i += 1
        elif c == "[":
            depth, i = 1, i + 1
            while i < n and depth:
                if text[i] == "[":
                    depth += 1
                elif text[i] == "]":
                    depth -= 1
                i += 1
            if depth:
                raise ValueError("unterminated comment in Newick")
        elif c == "]":
            raise ValueError("unmatched ']' in Newick")
        elif c == "'":
            j, buf, closed = i + 1, [], False
            while j < n:
                if text[j] == "'":
                    if j + 1 < n and text[j + 1] == "'":
                        buf.append("'")
                        j += 2
                        continue
                    closed = True
                    break
                buf.append(text[j])
                j += 1
            if not closed:
                raise ValueError("unterminated quote in Newick")
            yield "text", "".join(buf)
            i = j + 1
        else:
            j = i
            while j < n and text[j] not in "(),:; \t\r\n[]'":
                j += 1
            yield "text", text[i:j]
            i = j


def _parse_newick(text) -> _Node:
    """Parse Newick into a node tree without recursion (deep trees are fine).

    Malformed input raises ValueError: unbalanced parentheses, an unterminated
    quote or comment, a missing/repeated/non-numeric branch length, a missing
    comma, a second root, and any content after the terminating ``;``.

    States: ``subtree`` (need ``(`` or a leaf label), ``length`` (need a branch
    length), ``after_label`` / ``after_length`` / ``after_close`` (a node is
    complete; only a label after a close, then ``:``, ``,``, ``)`` or ``;``).
    """
    root = None
    stack = []
    last = None
    state = "subtree"
    tokens = _tokens(text)

    def add(node):
        nonlocal root
        if stack:
            stack[-1].children.append(node)
        elif root is None:
            root = node
        else:
            raise ValueError("more than one root in Newick")

    def close():
        nonlocal last
        if not stack:
            raise ValueError("unbalanced ')' in Newick")
        last = stack.pop()

    for kind, tok in tokens:
        punct = kind == "punct"
        if state == "subtree":
            if not punct:
                last = _Node()
                last.label = tok
                add(last)
                state = "after_label"
            elif tok == "(":
                node = _Node()
                add(node)
                stack.append(node)
            else:
                raise ValueError(f"unexpected token {tok!r} in Newick")
        elif state == "length":
            if punct:
                raise ValueError("expected a branch length after ':' in Newick")
            try:
                last.branch = float(tok)
            except ValueError:
                raise ValueError(f"invalid branch length: {tok!r}") from None
            state = "after_length"
        elif punct and tok == ";":
            break
        elif punct and tok == ":":
            if state == "after_length":
                raise ValueError("repeated branch length in Newick")
            state = "length"
        elif punct and tok == ",":
            if not stack:
                raise ValueError("comma outside a group in Newick")
            state = "subtree"
        elif punct and tok == ")":
            close()
            state = "after_close"
        elif punct:
            raise ValueError(f"unexpected token {tok!r} in Newick")
        elif state == "after_close":
            last.label = tok
            state = "after_label"
        else:
            raise ValueError(f"unexpected label {tok!r} in Newick")

    if state == "length":
        raise ValueError("missing branch length after ':' in Newick")
    if stack:
        raise ValueError("unbalanced '(' in Newick")
    if root is None:
        raise ValueError("empty Newick")
    if next(tokens, None) is not None:
        raise ValueError("trailing content after ';' in Newick")
    return root


class _Tree:
    """Pre-order-indexed tree with leaf ranges and O(log n) MRCA queries."""

    def __init__(self, root):
        order, parent, stack = [], {}, [root]
        while stack:
            node = stack.pop()
            order.append(node)
            for child in reversed(node.children):
                parent[id(child)] = node
                stack.append(child)
        n = len(order)
        index = {id(node): i for i, node in enumerate(order)}
        self.n = n
        self.label = [node.label for node in order]
        self.children = [[index[id(c)] for c in node.children] for node in order]
        self.parent = [index[id(parent[id(node)])] if id(node) in parent else -1 for node in order]

        self.lo = [0] * n
        self.hi = [0] * n
        self.leaf_pos = [0] * n
        leaf = 0
        self.leaf_order = []
        for i, children in enumerate(self.children):
            if not children:
                self.leaf_pos[i] = leaf
                self.lo[i] = self.hi[i] = leaf
                self.leaf_order.append(i)
                leaf += 1
        self.n_leaves = leaf
        for i in range(n - 1, -1, -1):
            children = self.children[i]
            if children:
                self.lo[i] = self.lo[children[0]]
                self.hi[i] = self.hi[children[-1]]

        up0 = [self.parent[i] if self.parent[i] >= 0 else i for i in range(n)]
        levels = max(1, n.bit_length())
        self.up = [up0]
        for _ in range(1, levels):
            prev = self.up[-1]
            self.up.append([prev[prev[i]] for i in range(n)])

    def mrca(self, min_leaf, max_leaf):
        """Pre-order index of the deepest node covering leaves [min_leaf, max_leaf]."""
        node = self.leaf_order[min_leaf]
        if self.hi[node] >= max_leaf:
            return node
        for up in reversed(self.up):
            ancestor = up[node]
            if self.hi[ancestor] < max_leaf:
                node = ancestor
        return self.up[0][node]


def _read_tracks(paths):
    """Return (column names in file order, list of {og: tuple(values)})."""
    columns, files, seen = [], [], set(_ORTHOGROUP_COLUMNS)
    for path in paths:
        with open_text(path) as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            names = header[1:]
            for name in names:
                if name in seen:
                    raise ValueError(f"duplicate track column name: {name!r}")
                seen.add(name)
            columns.extend(names)
            values = {}
            for line in fh:
                line = line.rstrip("\r\n")
                if not line:
                    continue
                cells = line.split("\t")
                row = []
                for cell in cells[1:1 + len(names)]:
                    if cell == "":
                        row.append(None)
                        continue
                    value = float(cell)
                    row.append(value if math.isfinite(value) else None)
                row += [None] * (len(names) - len(row))
                values[cells[0]] = tuple(row)
            files.append((names, values))
    return columns, files


def _ranks(n, key):
    order = sorted(range(n), key=key)
    ranks = [0] * n
    for rank, i in enumerate(order):
        ranks[i] = rank
    return ranks


def run(args) -> None:
    orthologs = Path(args.orthologs)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- tree ---------------------------------------------------------------
    if args.tree:
        tree_text = Path(args.tree).read_text()
        tree_source = "input"
    else:
        tree_text = upgma_newick(orthologs)
        tree_source = "upgma"
    tree = _Tree(_parse_newick(tree_text))

    # --- matrix labels + tree validation -----------------------------------
    labels, rows = read_matrix(orthologs)
    # `with rows` closes the matrix file even if validation or a row raises
    # before the iterator is exhausted.
    with rows:
        n_genomes = len(labels)
        dupe_labels = sorted(name for name, count in Counter(labels).items() if count > 1)
        if dupe_labels:
            raise ValueError(f"duplicate genome labels in matrix: {dupe_labels}")
        leaf_labels = [tree.label[i] for i in tree.leaf_order]
        dupe_leaves = sorted(name for name, count in Counter(leaf_labels).items() if count > 1)
        if dupe_leaves:
            raise ValueError(f"duplicate leaf labels in tree: {dupe_leaves}")
        leaf_names = set(leaf_labels)
        if leaf_names != set(labels):
            tree_only = sorted(leaf_names - set(labels))
            matrix_only = sorted(set(labels) - leaf_names)
            raise ValueError(
                "tree leaves do not match genome labels; "
                f"tree only: {tree_only}; matrix only: {matrix_only}"
            )
        leaf_pre = {tree.label[i]: i for i in tree.leaf_order}
        genome_leafpos = [tree.leaf_pos[leaf_pre[name]] for name in labels]

        # --- tracks ---------------------------------------------------------
        track_columns, track_files = _read_tracks(args.track)
        track_lookup = [(values, j) for names, values in track_files for j in range(len(names))]

        # --- stream rows ----------------------------------------------------
        og_ids = []
        og_n_genomes = []
        og_n_genes = []
        og_cn_mean = []
        og_cn_max = []
        og_mrca = []
        og_purity = []
        track_values = [[] for _ in track_columns]
        genome_n_og = [0] * n_genomes
        genome_n_genes = [0] * n_genomes
        matrix_ogs = set()

        presence = [bytearray() for _ in range(n_genomes)]
        presence_len = 0
        mem_og, mem_genome, mem_copies, mem_genes = [], [], [], []
        parts = []

        def flush_membership():
            if not mem_og:
                return
            lengths = np.fromiter((len(genes) for genes in mem_genes), np.int64)
            rows = np.repeat(np.arange(len(mem_og), dtype=np.int64), lengths)
            flat = pl.DataFrame(
                {"row": rows, "genes": pl.Series("genes", [g for genes in mem_genes for g in genes], pl.Utf8)}
            ).group_by("row", maintain_order=True).agg(pl.col("genes"))
            frame = (
                pl.DataFrame(
                    {
                        "row": pl.Series("row", np.arange(len(mem_og))),
                        "og_idx": pl.Series("og_idx", mem_og, dtype=pl.UInt32),
                        "genome_idx": pl.Series("genome_idx", mem_genome, dtype=pl.UInt32),
                        "copies": pl.Series("copies", mem_copies, dtype=pl.UInt32),
                    }
                )
                .join(flat, on="row", how="left")
                .sort("row")
                .drop("row")
            )
            part = outdir / f"_membership.{len(parts)}.parquet"
            frame.write_parquet(part)
            parts.append(part)
            del mem_og[:], mem_genome[:], mem_copies[:], mem_genes[:]

        for og, cells in rows:
            og_idx = len(og_ids)
            matrix_ogs.add(og)
            og_ids.append(og)

            byte = og_idx >> 3
            if byte + 1 > presence_len:
                new_len = ((byte >> 12) + 1) << 12
                pad = bytes(new_len - presence_len)
                for buf in presence:
                    buf.extend(pad)
                presence_len = new_len
            mask = 1 << (7 - (og_idx & 7))

            present = total = max_copies = 0
            min_leaf = max_leaf = None
            for genome_idx, cell in enumerate(cells):
                if not cell:
                    continue
                genes = tuple(split_genes(cell))
                if not genes:
                    continue
                copies = len(genes)
                present += 1
                total += copies
                if copies > max_copies:
                    max_copies = copies
                presence[genome_idx][byte] |= mask
                genome_n_og[genome_idx] += 1
                genome_n_genes[genome_idx] += copies
                pos = genome_leafpos[genome_idx]
                if min_leaf is None or pos < min_leaf:
                    min_leaf = pos
                if max_leaf is None or pos > max_leaf:
                    max_leaf = pos
                mem_og.append(og_idx)
                mem_genome.append(genome_idx)
                mem_copies.append(copies)
                mem_genes.append(genes)

            og_n_genomes.append(present)
            og_n_genes.append(total)
            if present:
                og_cn_mean.append(total / present)
                og_cn_max.append(float(max_copies))
                node = tree.mrca(min_leaf, max_leaf)
                og_mrca.append(node)
                og_purity.append(present / (tree.hi[node] - tree.lo[node] + 1))
            else:
                og_cn_mean.append(None)
                og_cn_max.append(None)
                og_mrca.append(None)
                og_purity.append(None)

            for out, (values, j) in zip(track_values, track_lookup):
                row = values.get(og)
                out.append(row[j] if row is not None else None)

            if len(mem_og) >= _MEMBERSHIP_BATCH:
                flush_membership()

        flush_membership()
        n_og = len(og_ids)

    # --- unknown track OGs warning -----------------------------------------
    if track_files:
        track_ogs = set()
        for _, values in track_files:
            track_ogs.update(values)
        unknown = track_ogs - matrix_ogs
        if unknown:
            print(
                f"warning: {len(unknown)} track orthogroup(s) not in the matrix (ignored)",
                file=sys.stderr,
            )

    # --- orderings ----------------------------------------------------------
    order_tree = _ranks(
        n_og,
        lambda i: (1, 0, 0.0, 0, i) if og_mrca[i] is None
        else (0, og_mrca[i], -og_purity[i], -og_n_genomes[i], i),
    )
    order_frequency = _ranks(n_og, lambda i: (-og_n_genomes[i], i))

    # --- write outputs ------------------------------------------------------
    (outdir / "tree.nwk").write_text(tree_text if tree_text.endswith("\n") else tree_text + "\n")

    nbytes = (n_og + 7) // 8
    raw = b"".join(bytes(buf[:nbytes]) for buf in presence)
    np.save(outdir / "presence.npy", np.frombuffer(raw, np.uint8).reshape(n_genomes, nbytes))
    del presence, raw  # release the bit buffers now that presence.npy exists

    orthogroups_schema = {
        "og_idx": pl.UInt32,
        "orthogroup": pl.Utf8,
        "n_genomes": pl.UInt32,
        "n_genes": pl.UInt32,
        "copy_number_mean": pl.Float64,
        "copy_number_max": pl.Float64,
        "mrca_node": pl.UInt32,
        "mrca_label": pl.Utf8,
        "purity": pl.Float64,
        "order_tree": pl.UInt32,
        "order_frequency": pl.UInt32,
    }
    for name in track_columns:
        orthogroups_schema[name] = pl.Float64
    orthogroups_data = {
        "og_idx": list(range(n_og)),
        "orthogroup": og_ids,
        "n_genomes": og_n_genomes,
        "n_genes": og_n_genes,
        "copy_number_mean": og_cn_mean,
        "copy_number_max": og_cn_max,
        "mrca_node": og_mrca,
        "mrca_label": [tree.label[m] if m is not None else None for m in og_mrca],
        "purity": og_purity,
        "order_tree": order_tree,
        "order_frequency": order_frequency,
    }
    for name, values in zip(track_columns, track_values):
        orthogroups_data[name] = values
    pl.DataFrame(orthogroups_data, schema=orthogroups_schema).write_parquet(
        outdir / "orthogroups.parquet"
    )

    pl.DataFrame(
        {
            "genome_idx": list(range(n_genomes)),
            "name": labels,
            "row_order": [genome_leafpos[i] for i in range(n_genomes)],
            "n_orthogroups": genome_n_og,
            "n_genes": genome_n_genes,
        },
        schema={
            "genome_idx": pl.UInt32,
            "name": pl.Utf8,
            "row_order": pl.UInt32,
            "n_orthogroups": pl.UInt32,
            "n_genes": pl.UInt32,
        },
    ).write_parquet(outdir / "genomes.parquet")

    membership_schema = {
        "og_idx": pl.UInt32,
        "genome_idx": pl.UInt32,
        "copies": pl.UInt32,
        "genes": pl.List(pl.Utf8),
    }
    membership = outdir / "membership.parquet"
    if len(parts) == 1:
        parts[0].replace(membership)
    elif parts:
        pl.scan_parquet([str(p) for p in parts]).sink_parquet(
            membership, row_group_size=_MEMBERSHIP_ROW_GROUP
        )
        for part in parts:
            part.unlink()
    else:
        pl.DataFrame(schema=membership_schema).write_parquet(membership)

    manifest = {
        "format_version": FORMAT_VERSION,
        "orthologs": str(orthologs),
        "tree": str(args.tree) if args.tree else None,
        "tree_source": tree_source,
        "tracks": [str(t) for t in args.track],
        "track_columns": track_columns,
        "n_genomes": n_genomes,
        "n_orthogroups": n_og,
        "n_membership": sum(og_n_genomes),
        "n_genes": sum(og_n_genes),
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def add_parser(subparsers):
    p = subparsers.add_parser("import", help="build a processed pangenome directory")
    p.add_argument("orthologs", help="orthologs TSV/TSV.GZ in per-genome matrix format")
    p.add_argument("-o", "--out", required=True, metavar="OUTDIR", help="output directory")
    p.add_argument("--tree", metavar="NEWICK", help="genome tree; default: UPGMA on Jaccard")
    p.add_argument("--track", action="append", default=[], metavar="TSV",
                   help="per-orthogroup track TSV (repeatable)")
    p.set_defaults(func=run)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_parser(parser.add_subparsers(dest="command", required=True))
    ns = parser.parse_args()
    ns.func(ns)
