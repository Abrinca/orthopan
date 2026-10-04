"""Byte-offset index for protein FASTA files (docs/design-decisions.md v2).

Sequences are never copied: the index stores ``(file, offset, length)`` per gene
and the slice is read on demand. Matching is on the first header word, keyed
both as written and with the gene-ID prefix stripped.

Gzip cannot be random-accessed, so ``.gz`` inputs are skipped with a warning;
decompress them if they must be indexed. OrthoFinder's
``Orthogroup_Sequences/<OG>.fa`` is read directly when available and needs no
index.

Duplicate IDs across files follow ``tracks/basic.py``'s rule: an ID present in
several files is ambiguous only when its sequence lengths differ; then it is
resolved by matching the FASTA file stem to the genome label. Same-length
duplicates keep the first occurrence.
"""

import sys
from pathlib import Path

from orthopan.io import strip_gene_prefix
from orthopan.tracks.basic import _stem, fasta_files

__all__ = [
    "fasta_files", "index_files", "resolve", "read_record", "is_gzip", "stem",
]


def stem(path) -> str:
    """File/label name without a trailing fasta extension (and optional .gz)."""
    return _stem(path)


def is_gzip(path) -> bool:
    """True for a ``.gz`` suffix (random access is impossible there)."""
    return Path(path).name.lower().endswith(".gz")


def index_files(files, prefix):
    """Return ``{key: [(file_idx, offset, length, seq_len), ...]}``.

    ``files`` are plain (non-gz) paths. ``prefix`` is a compiled gene-ID prefix
    (or None). Each record is keyed by its first header word and, when the
    prefix actually matches, by the stripped word too.
    """
    index: dict[str, list[tuple[int, int, int, int]]] = {}

    def add(key, record):
        if key:
            index.setdefault(key, []).append(record)

    for fi, path in enumerate(files):
        with open(path, "rb") as fh:
            start = None
            word = ""
            seq_len = 0
            while True:
                pos = fh.tell()
                line = fh.readline()
                if not line:
                    if start is not None:
                        _add(add, fi, start, pos - start, seq_len, word, prefix)
                    break
                if line[:1] == b">":
                    if start is not None:
                        _add(add, fi, start, pos - start, seq_len, word, prefix)
                    start = pos
                    text = line[1:].strip()
                    word = text.split(None, 1)[0].decode() if text else ""
                    seq_len = 0
                elif start is not None:
                    seq_len += sum(1 for c in line if c > 32)  # non-whitespace ASCII
    return index


def _add(add, fi, start, length, seq_len, word, prefix):
    record = (fi, start, length, seq_len)
    add(word, record)
    stripped = strip_gene_prefix(word, prefix)
    if stripped != word:
        add(stripped, record)


def resolve(index, raw, display, genome, files, stems):
    """Resolve one membership gene to a record, or None.

    Duplicates across files are resolved by preferring a unique genome/file-stem
    match; otherwise they are ambiguous unless all candidates carry identical
    sequences (then the first wins). ``stems`` maps a file stem to a genome.
    """
    for key in (raw, display):
        records = index.get(key)
        if not records:
            continue
        if len(records) == 1:
            return records[0]
        matches = [r for r in records if stems.get(stem(files[r[0]])) == genome]
        if len(matches) == 1:
            return matches[0]
        if len({r[3] for r in records}) == 1:
            sequences = [read_record(files, r)[1] for r in records]
            if len(set(sequences)) == 1:
                return records[0]
        return None
    return None


def read_record(files, record):
    """Return ``(header_word, sequence)``.

    ``record`` is either a :func:`index_files` 4-tuple or a
    ``(file, offset, length)`` triple from ``sequences.parquet``.
    """
    fi, offset, length = record[0], record[1], record[2]
    with open(files[fi], "rb") as fh:
        fh.seek(offset)
        data = fh.read(length)
    header, _, body = data.partition(b"\n")
    word = header[1:].strip().split(None, 1)[0].decode()
    return word, b"".join(body.split()).decode()


def warn_gzip(paths):
    """Return the plain files of ``paths`` (dirs expanded), warning for gzip."""
    plain = []
    for path in paths:
        if is_gzip(path):
            print(f"warning: skipping gzip FASTA (cannot be indexed): {path}",
                  file=sys.stderr)
        else:
            plain.append(str(path))
    return plain
