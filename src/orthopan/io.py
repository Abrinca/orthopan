"""Shared readers for orthologs TSV inputs (see docs/ortholog_file_formats_documentation.md)."""

import gzip
from collections.abc import Iterator
from itertools import chain


def open_text(path):
    """Open plain or gzip text (gzip detected by magic bytes)."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


class _MatrixRows:
    """Iterator over matrix body rows that owns the input file.

    Iterating to exhaustion closes the file; ``with``/``close()`` close it
    deterministically when the caller stops early (for example after raising).
    """

    def __init__(self, fh, rows, n):
        self._fh = fh
        self._rows = rows
        self._n = n
        self._gen = self._generate()

    def _generate(self):
        try:
            for line_no, line in self._rows:
                og, *cells = line.split("\t")
                if len(cells) > self._n:
                    raise ValueError(
                        f"line {line_no}: {len(cells)} genome columns, expected {self._n}"
                    )
                yield og, cells + [""] * (self._n - len(cells))
        finally:
            self._fh.close()

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._gen)

    def close(self):
        self._fh.close()
        self._gen.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


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

    return labels, _MatrixRows(fh, rows, n)


def split_genes(cell: str) -> list[str]:
    """Gene IDs in a matrix cell ("g1, g2"); empty cell -> []."""
    return [g.strip() for g in cell.split(",") if g.strip()]
