"""Shared readers for orthologs TSV inputs (see docs/ortholog_file_formats_documentation.md)."""

import gzip
from collections.abc import Iterator
from itertools import chain


def open_text(path):
    """Open plain or gzip text (gzip detected by magic bytes)."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_matrix(path) -> tuple[list[str], Iterator[tuple[str, list[str]]]]:
    """Read a Format B (per-genome matrix) TSV.

    Returns (genome labels, iterator of (orthogroup, cells)). Cells are padded to
    the genome count; overwide rows raise ValueError with the line number.
    Headerless input gets labels genome_1..genome_n and keeps its first row.
    """
    fh = open_text(path)
    lines = (ln.rstrip("\r\n") for ln in fh)
    numbered = ((i, ln) for i, ln in enumerate(lines, 1) if ln)
    first = next(numbered, None)
    if first is None:
        fh.close()
        raise ValueError("input has no genome columns")
    head = first[1].split("\t")
    if head[0] in ("Orthogroup", "HOG"):
        labels, rows = head[1:], numbered
    else:
        labels, rows = [f"genome_{i}" for i in range(1, len(head))], chain([first], numbered)
    n = len(labels)
    if not n:
        fh.close()
        raise ValueError("input has no genome columns")

    def gen():
        with fh:
            for line_no, line in rows:
                og, *cells = line.split("\t")
                if len(cells) > n:
                    raise ValueError(f"line {line_no}: {len(cells)} genome columns, expected {n}")
                yield og, cells + [""] * (n - len(cells))

    return labels, gen()


def split_genes(cell: str) -> list[str]:
    """Gene IDs in a matrix cell ("g1, g2"); empty cell -> []."""
    return [g.strip() for g in cell.split(",") if g.strip()]
