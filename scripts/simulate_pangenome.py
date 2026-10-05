#!/usr/bin/env python3
"""Simulate a pangenome at scale for performance measurement.

Deterministic (seeded). Writes a rooted genome tree, a Format B orthologs TSV
(gzip), optional orthogroup annotations and genome metadata.

The tree is Yule-like (random pair merging of lineages), which is reasonably
balanced, and every edge has a positive branch length so a phylogram viewer can
lay it out. Presence evolves along the tree: each orthogroup is gained in a
clade (found by walking up from a random leaf) and can be lost again in
subclades, so the matrix shows clade blocks at several tree levels, not one
nested staircase. A small core and a few paralogs are kept.

The TSV is streamed orthogroup by orthogroup; no genome x orthogroup object
matrix is ever held in memory.
"""

import argparse
import gzip
import random
import sys
from pathlib import Path

import numpy as np

COG = list("JKLMNOQUVWXYZABCDEFGHIPST")  # ~24 COG-like letters


class Node:
    __slots__ = ("children", "branch", "label", "lo", "hi")

    def __init__(self):
        self.children = []
        self.branch = None
        self.label = ""
        self.lo = 0
        self.hi = 0


def _yule_tree(n, rng):
    """Random rooted binary tree (random pair merging). Returns (root, leaves).

    Every edge gets a positive, bounded branch length (lognormal, median 1), so
    no single branch dominates the layout.
    """
    leaves = [Node() for _ in range(n)]
    active = list(leaves)
    while len(active) > 1:
        a = active.pop(rng.randrange(len(active)))
        b = active.pop(rng.randrange(len(active)))
        parent = Node()
        a.branch = _branch(rng)
        b.branch = _branch(rng)
        parent.children = [a, b]
        active.append(parent)
    root = active[0]
    root.branch = None
    return root, leaves


def _branch(rng):
    value = rng.lognormvariate(0.0, 0.4)
    return min(6.0, max(0.02, value))


def _serialize(root):
    """Pre-order Newick, with branch lengths, and each leaf's pre-order position."""
    out = []
    leaf_pos = {}  # id(node) -> position
    pos = 0
    stack = [(root, None)]
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            out.append(item)
            continue
        node, ph = item
        suffix = "" if node.branch is None else f":{node.branch:.6g}"
        if not node.children:
            out.append(node.label)
            out.append(suffix)
            leaf_pos[id(node)] = pos
            pos += 1
        else:
            out.append("(")
            stack.append(")" + suffix)
            stack.append((node.children[1], None))
            stack.append(",")
            stack.append((node.children[0], None))
    return "".join(out), leaf_pos


def _tree_index(root, leaf_pos):
    """Return (order, parent, lo, hi, span, leaf_at) in pre-order node indices.

    ``leaf_at[pos]`` is the node index of the leaf at tree position ``pos``.
    """
    order, stack = [], [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(node.children)
    n_nodes = len(order)
    index = {id(node): i for i, node in enumerate(order)}
    parent = [-1] * n_nodes
    lo = [0] * n_nodes
    hi = [0] * n_nodes
    n_leaves = 0
    for i, node in enumerate(order):
        if node.children:
            for child in node.children:
                parent[index[id(child)]] = i
        else:
            lo[i] = hi[i] = leaf_pos[id(node)]
            n_leaves += 1
    for i in reversed(range(n_nodes)):
        if order[i].children:
            lo[i] = lo[index[id(order[i].children[0])]]
            hi[i] = hi[index[id(order[i].children[-1])]]
    span = [hi[i] - lo[i] + 1 for i in range(n_nodes)]
    leaf_at = [0] * n_leaves
    for i, node in enumerate(order):
        if not node.children:
            leaf_at[lo[i]] = i
    return order, parent, lo, hi, span, leaf_at, n_leaves


def _interval_add(intervals, a, b):
    """Add [a, b] to a sorted list of disjoint inclusive intervals."""
    out = []
    placed = False
    for x, y in intervals:
        if b < x - 1:
            if not placed:
                out.append((a, b))
                placed = True
            out.append((x, y))
        elif y < a - 1:
            out.append((x, y))
        else:  # overlap / touch -> merge
            a, b = min(a, x), max(b, y)
    if not placed:
        out.append((a, b))
    return out


def _interval_remove(intervals, a, b):
    """Subtract [a, b] from a sorted list of disjoint inclusive intervals."""
    out = []
    for x, y in intervals:
        if b < x or y < a:  # no overlap
            out.append((x, y))
            continue
        if x < a:
            out.append((x, a - 1))
        if y > b:
            out.append((b + 1, y))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genomes", type=int, default=10_000)
    p.add_argument("--orthogroups", type=int, default=50_000)
    p.add_argument("--out", required=True, type=Path, help="output directory")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--core-frac", type=float, default=0.03)
    p.add_argument("--paralog-frac", type=float, default=0.03)
    a = p.parse_args(argv)

    n, m = a.genomes, a.orthogroups
    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(a.seed)
    nprng = np.random.default_rng(a.seed)

    # --- tree -----------------------------------------------------------------
    root, leaves = _yule_tree(n, rng)
    for i, leaf in enumerate(leaves):
        leaf.label = f"genome_{i:05d}"
    text, leaf_pos = _serialize(root)
    nodes, parent, lo, hi, span, leaf_at, n_leaves = _tree_index(root, leaf_pos)
    (out / "tree.nwk").write_text(text + ";\n")
    assert n_leaves == n

    # genome index == pre-order leaf position, so header order == tree order.
    pos_to_label = {leaf_pos[id(leaf)]: leaf.label for leaf in leaves}
    labels = [pos_to_label[i] for i in range(n)]

    # --- presence-count model -------------------------------------------------
    # ~core-frac present in (nearly) every genome; the rest heavy-tailed, so the
    # matrix has a long tail of rare accessory genes and singletons.
    core = nprng.random(m) < a.core_frac
    counts = np.ones(m, dtype=np.int64)
    acc = ~core
    raw = nprng.lognormal(mean=np.log(31.0), sigma=2.0, size=int(acc.sum()))
    counts[acc] = np.clip(raw.astype(np.int64), 1, int(0.9 * n))
    core_present = nprng.integers(int(0.985 * n), n + 1, size=int(core.sum()))
    counts[core] = core_present
    n_core = int(core.sum())

    # --- gene counters --------------------------------------------------------
    gene_counter = [0] * n

    def next_id(genome_idx):
        c = gene_counter[genome_idx]
        gene_counter[genome_idx] = c + 1
        return f"{labels[genome_idx]}_{c:06d}"

    total_genes = 0
    total_cells = 0
    desc_path = out / "annotations.descriptions.tsv"
    cog_path = out / "annotations.cog.tsv"
    fh = gzip.open(out / "orthologs.tsv.gz", "wt", compresslevel=1)
    fh.write("Orthogroup\t" + "\t".join(labels) + "\n")
    fd = open(desc_path, "w")
    fc = open(cog_path, "w")
    fc.write("orthogroup\tcog_category\n")

    for og in range(m):
        c = int(counts[og])
        if core[og]:
            intervals = [(0, n - 1)]
        else:
            intervals = []
            budget = 0
            # Gain in up to a few clades found by walking up from random leaves.
            tries = 0
            while budget < c and tries < 8:
                tries += 1
                node = leaf_at[rng.randrange(n)]
                while parent[node] != -1 and span[parent[node]] <= c:
                    node = parent[node]
                intervals = _interval_add(intervals, lo[node], hi[node])
                budget = sum(y - x + 1 for x, y in intervals)
            # Trim overshoot from the last clade, keeping counts near target.
            while budget > c and intervals:
                x, y = intervals[-1]
                take = min(budget - c, y - x + 1)
                if take == y - x + 1:
                    intervals = intervals[:-1]
                else:
                    intervals[-1] = (x, y - take)
                budget -= take
            # Loss in a random subclade, so blocks are not solid rectangles.
            if rng.random() < 0.5 and budget > 2:
                i = rng.randrange(len(intervals))
                x, y = intervals[i]
                node = leaf_at[rng.randrange(x, y + 1)]
                limit = max(1, int((y - x + 1) * rng.uniform(0.1, 0.5)))
                while (parent[node] != -1 and span[parent[node]] <= limit
                       and lo[parent[node]] >= x and hi[parent[node]] <= y):
                    node = parent[node]
                intervals = _interval_remove(intervals, lo[node], hi[node])
            if not intervals:  # never emit an all-empty row
                single = leaf_at[rng.randrange(n)]
                intervals = [(lo[single], hi[single])]

        present = np.concatenate(
            [np.arange(x, y + 1, dtype=np.int64) for x, y in intervals]
        )
        cells = [""] * n
        extra = nprng.random(present.size) < a.paralog_frac
        for k, pos in enumerate(present):
            g = int(pos)
            ids = [next_id(g)]
            if extra[k]:
                ids.append(next_id(g))
            cells[g] = ", ".join(ids)
            total_genes += len(ids)
        total_cells += present.size
        og_id = f"OG{og:07d}"
        fh.write(og_id + "\t" + "\t".join(cells) + "\n")
        if rng.random() < 0.5:
            fd.write(f"{og_id}\thypothetical protein\n")
        else:
            fd.write(f"{og_id}\t{rng.choice(['transporter', 'kinase', 'transferase', 'regulator', 'hydrolase'])}\n")
        fc.write(f"{og_id}\t{COG[og % len(COG)]}\n")

    fh.close()
    fd.close()
    fc.close()

    # --- genome metadata: clade (categorical) + gc (numeric) -------------------
    # Clade is a coarse split of tree order, so the metadata column matches it.
    fh = open(out / "genomes.meta.tsv", "w")
    fh.write("genome\tclade\tgc\n")
    half = n // 2
    clades = ["B"] * n
    for i in range(half):
        clades[i] = "A"
    for start in range(0, n, max(1, n // 8)):
        width = max(1, n // 32)
        clades[start:start + width] = ["C"] * len(clades[start:start + width])
    for i, name in enumerate(labels):
        gc = 0.30 + 0.12 * (i / max(1, n - 1)) + nprng.normal(0, 0.005)
        fh.write(f"{name}\t{clades[i]}\t{gc:.4f}\n")
    fh.close()

    size = (out / "orthologs.tsv.gz").stat().st_size
    print(
        f"genomes={n} orthogroups={m} total_genes={total_genes} "
        f"present_cells={total_cells} core_ogs={n_core} "
        f"tsv_gz_bytes={size}",
        file=sys.stderr,
    )
    # quick self-check: streamed row count matches, and (cheap on small runs)
    # the non-empty cells actually written match the modelled present count.
    assert len(labels) == n
    if m <= 2000:
        with gzip.open(out / "orthologs.tsv.gz", "rt") as fh:
            fh.readline()
            written = sum(
                1 for line in fh
                for cell in line.rstrip("\n").split("\t")[1:] if cell
            )
        assert written == total_cells, (written, total_cells)


if __name__ == "__main__":
    main()
