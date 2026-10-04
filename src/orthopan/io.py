"""Shared readers for orthologs TSV inputs (see docs/ortholog_file_formats_documentation.md)."""

import gzip
import re
from collections.abc import Iterator
from itertools import chain

# Default gene-ID prefix (docs/ortholog_file_formats_documentation.md §3D).
DEFAULT_GENE_PREFIX = r"^gnl\|[^|]+\|"


def open_text(path):
    """Open plain or gzip text (gzip detected by magic bytes)."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def compile_prefix(pattern):
    """Compile a strip-prefix regex; ``None``/``""`` disables stripping."""
    return re.compile(pattern) if pattern else None


def strip_gene_prefix(gene: str, pattern) -> str:
    """Remove the configured gene-ID prefix (a compiled pattern or ``None``)."""
    return pattern.sub("", gene) if pattern else gene


def genome_from_gene(gene: str, pattern) -> str:
    """Genome label of a gene: strip prefix, then drop the trailing ``_<n>``."""
    return strip_gene_prefix(gene, pattern).rsplit("_", 1)[0]


class _OwnedRows:
    """Iterator over body rows that owns the input file.

    Iterating to exhaustion closes the file; ``with``/``close()`` close it
    deterministically when the caller stops early (for example after raising).
    """

    def __init__(self, fh):
        self._fh = fh
        self._gen = self._generate()

    def _generate(self):
        raise NotImplementedError

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


class _MatrixRows(_OwnedRows):
    """Format B rows: ``og, [cells]`` padded to the genome count."""

    def __init__(self, fh, rows, n):
        self._rows = rows
        self._n = n
        super().__init__(fh)

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


class _FormatARows(_OwnedRows):
    """Format A rows: ``og, [per-genome gene lists]``.

    The genome label is derived from each gene ID; genes are bucketed by genome
    in first-appearance label order (fixed up front by :func:`read_matrix`).
    Cells keep the gene IDs as written (prefix included).
    """

    def __init__(self, fh, rows, n, index, pattern):
        self._rows = rows
        self._n = n
        self._index = index
        self._pattern = pattern
        super().__init__(fh)

    def _generate(self):
        try:
            for line_no, line in self._rows:
                og, _, field = line.partition("\t")
                buckets = [[] for _ in range(self._n)]
                for gene in split_genes(field):
                    label = genome_from_gene(gene, self._pattern)
                    idx = self._index.get(label)
                    if idx is None:
                        raise ValueError(
                            f"line {line_no}: unknown genome {label!r} in gene {gene!r}"
                        )
                    buckets[idx].append(gene)
                yield og, [", ".join(b) for b in buckets]
        finally:
            self._fh.close()


def _format_a_labels(path, pattern):
    """First pass: collect genome labels in first-appearance order."""
    labels, index = [], {}
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            _, _, field = line.partition("\t")
            for gene in split_genes(field):
                label = genome_from_gene(gene, pattern)
                if label not in index:
                    index[label] = len(labels)
                    labels.append(label)
    return labels, index


def read_matrix(path, strip_prefix: str | None = DEFAULT_GENE_PREFIX):
    """Read an orthologs TSV.

    Returns (genome labels, iterator of (orthogroup, cells)). Cells are padded to
    the genome count. Two layouts are supported:

    - Format B (per-genome matrix): header ``Orthogroup<TAB>g1<TAB>...`` or a
      headerless matrix with more than two columns. Headerless input gets labels
      genome_1..genome_n and keeps its first row.
    - Format A (flat two-column): headerless ``HOG<TAB>gene1, gene2, ...``. The
      genome comes from each gene ID: strip ``strip_prefix``, then take
      everything before the last underscore. A first pass collects the genome
      labels; genes are distributed into per-genome cells.

    Overwide Format B rows raise ValueError with the line number.
    """
    pattern = compile_prefix(strip_prefix)
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
        n = len(labels)
        if not n:
            fh.close()
            raise ValueError("input has no genome columns")
        return labels, _MatrixRows(fh, rows, n)
    if len(head) == 2:
        # Format A: exactly two columns and no header.
        fh.close()
        labels, index = _format_a_labels(path, pattern)
        if not labels:
            raise ValueError("input has no genome columns")
        fh = open_text(path)
        numbered = ((i, ln) for i, ln in enumerate(
            (ln.rstrip("\r\n") for ln in fh), 1) if ln)
        return labels, _FormatARows(fh, numbered, len(labels), index, pattern)
    labels, rows = [f"genome_{i}" for i in range(1, len(head))], chain([first], numbered)
    n = len(labels)
    if not n:
        fh.close()
        raise ValueError("input has no genome columns")
    return labels, _MatrixRows(fh, rows, n)


def split_genes(cell: str) -> list[str]:
    """Gene IDs in a matrix cell ("g1, g2"); empty cell -> []."""
    return [g.strip() for g in cell.split(",") if g.strip()]
