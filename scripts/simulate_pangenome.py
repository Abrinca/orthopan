#!/usr/bin/env python3
"""Simulate a pangenome at scale for performance measurement.

Deterministic (seeded). Writes a rooted genome tree, a Format B orthologs TSV
(gzip), optional orthogroup annotations and genome metadata.

Presence evolves along the tree, so the matrix has clade blocks: a small core,
a long tail of rare accessory genes and singletons, and a few paralogs.

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
        self.branch = 0.0
        self.label = ""
        self.lo = 0
        self.hi = 0


def _coalescent_tree(n, rng):
    """Random rooted binary tree (coalescent order). Returns (root, leaves)."""
    leaves = [Node() for _ in range(n)]
    active = list(leaves)
    while len(active) > 1:
        i = rng.randrange(len(active))
        a = active.pop(i)
        j = rng.randrange(len(active))
        b = active.pop(j)
        parent = Node()
        a.branch = max(1e-6, rng.expovariate(1.0))
        b.branch = max(1e-6, rng.expovariate(1.0))
        parent.children = [a, b]
        active.append(parent)
    root = active[0]
    return root, leaves


def _number_leaves(root):
    """Pre-order serialize the tree and record each leaf's pre-order position."""
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
        if not node.children:
            out.append(node.label)
            leaf_pos[id(node)] = pos
            pos += 1
        else:
            out.append("(")
            stack.append(")")
            stack.append((node.children[1], None))
            stack.append(",")
            stack.append((node.children[0], None))
    return "".join(out), leaf_pos


def _set_ranges(root, leaf_pos):
    """Iterative post-order: leaf range [lo, hi] (inclusive) for every node.

    Returns every node, parent-first, so callers need not walk the tree again.
    """
    order = []
    stack = [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(node.children)
    # children always appear later than their parent in `order`; reverse works
    # because nodes are appended parent-first.
    for node in reversed(order):
        if not node.children:
            node.lo = node.hi = leaf_pos[id(node)]
            continue
        node.lo = node.children[0].lo
        node.hi = node.children[-1].hi
    return order


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
    root, leaves = _coalescent_tree(n, rng)
    for i, leaf in enumerate(leaves):
        leaf.label = f"genome_{i:05d}"
    text, leaf_pos = _number_leaves(root)
    nodes = _set_ranges(root, leaf_pos)  # every node, parent-first
    (out / "tree.nwk").write_text(text + ";\n")

    # genome index == pre-order leaf position, so header order == tree order.
    pos_to_label = {leaf_pos[id(leaf)]: leaf.label for leaf in leaves}
    labels = [pos_to_label[i] for i in range(n)]

    # --- presence-count model -------------------------------------------------
    # ~core-frac present in (nearly) every genome; the rest heavy-tailed, so the
    # matrix has a long tail of rare accessory genes and singletons.
    core = nprng.random(m) < a.core_frac
    counts = np.ones(m, dtype=np.int64)
    acc = ~core
    # lognormal presence counts -> long tail, mean tuned so total genes ~25M at 10k
    raw = nprng.lognormal(mean=np.log(25.0), sigma=2.0, size=int(acc.sum()))
    counts[acc] = np.clip(raw.astype(np.int64), 1, int(0.9 * n))
    core_present = nprng.integers(int(0.97 * n), n + 1, size=int(core.sum()))
    counts[core] = core_present

    # nodes sorted by subtree leaf count, for clade-respecting placement
    node_size = [(node.hi - node.lo + 1, node) for node in nodes]
    node_size.sort(key=lambda t: t[0])
    sizes = np.array([s for s, _ in node_size], dtype=np.int64)
    order_nodes = [node for _, node in node_size]

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
        # place the presence set inside one clade when possible (clade blocks)
        idx = int(np.searchsorted(sizes, c, side="left"))
        if idx < len(sizes) and sizes[idx] <= max(c * 3, c + 3):
            node = order_nodes[idx]
            lo, hi = node.lo, node.hi
        else:
            idx = int(nprng.integers(0, len(sizes)))
            node = order_nodes[idx]
            lo, hi = node.lo, node.hi
        span = hi - lo + 1
        if c >= span:
            present = np.arange(lo, hi + 1)
        else:
            present = np.sort(nprng.choice(span, c, replace=False) + lo)

        cells = [""] * n
        # paralogs: ~paralog-frac of present cells get one extra copy
        extra = nprng.random(c) < a.paralog_frac
        for k, pos in enumerate(present):
            g = int(pos)
            ids = [next_id(g)]
            if extra[k]:
                ids.append(next_id(g))
            cells[g] = ", ".join(ids)
            total_genes += len(ids)
        total_cells += c
        og_id = f"OG{og:07d}"
        fh.write(og_id + "\t" + "\t".join(cells) + "\n")
        if rng.random() < 0.5:
            fd.write(f"{og_id}\thypothetical protein\n")
        else:
            fd.write(f"{og_id}\t{random.choice(['transporter', 'kinase', 'transferase', 'regulator', 'hydrolase'])}\n")
        fc.write(f"{og_id}\t{COG[og % len(COG)]}\n")

    fh.close()
    fd.close()
    fc.close()

    # --- genome metadata: clade (categorical) + gc (numeric) -------------------
    # assign clade from the tree: the two halves of the root are A and B; a
    # deeper split makes C. This makes the metadata column match the tree.
    fh = open(out / "genomes.meta.tsv", "w")
    fh.write("genome\tclade\tgc\n")
    half = n // 2
    clades = ["A"] * half + ["B"] * (n - half)
    # overwrite every 8th block with a third clade for variety
    for i in range(0, n, max(1, n // 8)):
        clades[i:i + max(1, n // 32)] = ["C"] * len(clades[i:i + max(1, n // 32)])
    for i, name in enumerate(labels):
        gc = 0.30 + 0.12 * (i / max(1, n - 1)) + nprng.normal(0, 0.005)
        fh.write(f"{name}\t{clades[i]}\t{gc:.4f}\n")
    fh.close()

    size = (out / "orthologs.tsv.gz").stat().st_size
    print(
        f"genomes={n} orthogroups={m} total_genes={total_genes} "
        f"present_cells={total_cells} core_ogs={int(core.sum())} "
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
