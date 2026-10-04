"""Build a genome dendrogram (Newick) from an orthologs presence/absence matrix."""

import argparse
import gzip
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

from orthopan.io import DEFAULT_GENE_PREFIX, open_text, read_matrix

_SPECIAL = set("(),:;'[]")


def _is_format_a(path):
    """Headerless two-column input -> Format A (genome from the gene ID)."""
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            head = line.split("\t")
            return head[0] not in ("Orthogroup", "HOG") and len(head) == 2
    return False


def _matrix(path, strip_prefix=DEFAULT_GENE_PREFIX):
    """Return (genome labels, presence matrix with genomes as rows)."""
    if _is_format_a(path):
        labels, rows = read_matrix(path, strip_prefix=strip_prefix)
        n = len(labels)
        arr = np.array([[bool(c) for c in cells] for _, cells in rows], bool)
        if arr.size == 0:
            return labels, np.zeros((n, 0), bool)
        return labels, arr.T
    with open(path, "rb") as fh:
        magic = fh.read(2) == b"\x1f\x8b"
    labels, n, buf = [], None, bytearray()
    with (gzip.open if magic else open)(path, "rb") as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.rstrip(b"\r\n")
            if not line:
                continue
            if n is None and line.startswith((b"Orthogroup\t", b"HOG\t")):
                labels = [c.decode() for c in line.split(b"\t")[1:]]
                n = len(labels)
                continue
            cuts = np.flatnonzero(np.frombuffer(line, np.uint8) == 9)
            row = (np.diff(np.concatenate(([-1], cuts, [len(line)]))) > 1)[1:]
            if n is None:
                n = len(row)
                labels = labels or [f"genome_{i}" for i in range(1, n + 1)]
            elif len(row) > n:
                raise ValueError(f"line {line_no}: {len(row)} genome columns, expected {n}")
            buf += row.astype(np.uint8).tobytes()
            buf += b"\x00" * (n - len(row))  # pad dropped trailing empty cells
    if not n:
        raise ValueError("input has no genome columns")
    return labels, np.frombuffer(buf, dtype=bool).reshape(-1, n).T


def _jaccard(matrix):
    """Condensed Jaccard distances (rows); equals scipy pdist, but fast."""
    x = matrix.astype(np.float32)
    inter = x @ x.T
    union = x.sum(1)[:, None] + x.sum(1)[None, :] - inter
    dist = 1.0 - np.divide(inter, union, out=np.ones_like(inter), where=union > 0)
    np.fill_diagonal(dist, 0.0)
    return squareform(dist, checks=False)


def _quote(name):
    if name and not any(ch.isspace() or ch in _SPECIAL for ch in name):
        return name
    return "'" + name.replace("'", "''") + "'"


def _newick(root, labels):
    """Iterative Newick writer (avoids RecursionError on deep trees)."""
    out, stack = [], [(root, None)]
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            out.append(item)
            continue
        node, ph = item
        if node.is_leaf():
            out.append(_quote(labels[node.id]))
            out.append("" if ph is None else f":{(ph - node.dist) / 2:.6g}")
        else:
            out.append("(")
            stack.append(")" if ph is None else f"):{(ph - node.dist) / 2:.6g}")
            stack.append((node.right, node.dist))
            stack.append(",")
            stack.append((node.left, node.dist))
    return "".join(out)


def upgma_newick(path, strip_prefix=DEFAULT_GENE_PREFIX) -> str:
    """Return Newick with UPGMA (average linkage) on Jaccard distances."""
    labels, matrix = _matrix(path, strip_prefix)
    if len(labels) < 2:
        return _quote(labels[0]) + ";"
    root = to_tree(linkage(_jaccard(matrix), method="average"))
    return _newick(root, labels) + ";"


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Genome dendrogram (UPGMA/Jaccard) -> Newick")
    p.add_argument("input", help="orthologs TSV/TSV.GZ in per-sample matrix format")
    p.add_argument("-o", "--output", help="output Newick path (default: stdout)")
    a = p.parse_args(argv)
    nwk = upgma_newick(a.input)
    if a.output:
        Path(a.output).write_text(nwk + "\n")
    else:
        print(nwk)


if __name__ == "__main__":
    main()
