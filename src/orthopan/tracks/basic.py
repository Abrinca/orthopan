"""Basic per-orthogroup tracks: paralog fraction and protein-length variation.

Copy number is **not** here: ``orthopan import`` always computes
``copy_number_mean``/``copy_number_max``, and a track with those names would
collide with the built-in columns. Works with any orthology tool in Format A or
B; parsing is delegated to :mod:`orthopan.io`. ``length_cv`` needs protein FASTA
(``--fasta``) to look up gene lengths; gene ID is the first whitespace-delimited
word of the header. Only lengths are held, not sequences.
"""

import argparse
import math
import sys
from pathlib import Path

from orthopan.io import open_text, read_matrix, split_genes
from orthopan.tracks import write_tracks

COLUMNS = ["paralog_fraction", "length_cv"]

_FASTA_SUFFIXES = (".fa", ".faa", ".fasta")


def _is_fasta(name: str) -> bool:
    name = name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return name.endswith(_FASTA_SUFFIXES)


def fasta_files(paths):
    """Expand files and directories (*.fa, *.faa, *.fasta, optionally .gz)."""
    for path in paths:
        path = Path(path)
        if path.is_dir():
            yield from sorted(p for p in path.iterdir() if p.is_file() and _is_fasta(p.name))
        else:
            yield path


def _stem(path) -> str:
    """File/label name without a trailing fasta extension (and optional .gz)."""
    name = Path(path).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    low = name.lower()
    for ext in _FASTA_SUFFIXES:
        if low.endswith(ext):
            return name[: -len(ext)]
    return name


def _read_fasta(path):
    """Yield ``(gene, length)`` for one protein FASTA (length excludes ``*``/space)."""
    gene, n = None, 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if gene is not None:
                    yield gene, n
                parts = line[1:].split()
                gene = parts[0] if parts else None
                n = 0
            elif gene is not None:
                n += sum(1 for ch in line if not ch.isspace() and ch != "*")
    if gene is not None:
        yield gene, n


def load_lengths(paths, labels=()) -> dict:
    """Map gene ID -> length from protein FASTAs.

    Unique gene IDs are keyed alone (common, low-memory path). An ID that occurs
    with different lengths across proteomes is ambiguous; it is keyed by
    ``(genome, gene)``, where ``genome`` is the matrix column label whose name
    matches the FASTA file stem (extension stripped). Ambiguous IDs that cannot
    be mapped to a genome are dropped and counted in a stderr warning.
    """
    files = list(fasta_files(paths))
    lengths: dict[str, int] = {}
    ambiguous: set[str] = set()
    for path in files:
        for gene, n in _read_fasta(path):
            if gene not in lengths:
                lengths[gene] = n
            elif lengths[gene] != n:
                ambiguous.add(gene)
    if not ambiguous:
        return lengths
    for gene in ambiguous:
        del lengths[gene]

    stems = {}
    for label in labels:
        stems.setdefault(_stem(label), label)
    unresolved = 0
    for path in files:
        genome = stems.get(_stem(path))
        for gene, n in _read_fasta(path):
            if gene not in ambiguous:
                continue
            if genome is None:
                unresolved += 1
            else:
                lengths[(genome, gene)] = n
    if unresolved:
        print(f"warning: {unresolved} ambiguous gene length(s) could not be mapped to a genome",
              file=sys.stderr)
    return lengths


def _length_cv(vals) -> float | None:
    """Population CV of member-gene lengths, or None if fewer than 2 are known."""
    vals = [v for v in vals if v is not None]
    if len(vals) < 2:
        return None
    mean = sum(vals) / len(vals)
    if mean == 0:
        return None
    var = sum((v - mean) ** 2 for v in vals) / len(vals)
    return math.sqrt(var) / mean


def basic_tracks(path, fasta=()) -> dict[str, dict[str, float | None]]:
    """Return {orthogroup: {track: value}} for the two basic tracks."""
    labels, rows = read_matrix(path)
    lengths = load_lengths(fasta, labels) if fasta else {}
    result: dict[str, dict[str, float | None]] = {}
    for og, cells in rows:
        copies = [cell.count(",") + 1 for cell in cells if cell]
        if not copies:
            result[og] = dict.fromkeys(COLUMNS)
            continue
        present = len(copies)
        values: dict[str, float | None] = {
            "paralog_fraction": sum(c > 1 for c in copies) / present,
            "length_cv": None,
        }
        if lengths:
            vals = [lengths.get((labels[i], g), lengths.get(g))
                    for i, cell in enumerate(cells)
                    for g in split_genes(cell)]
            values["length_cv"] = _length_cv(vals)
        result[og] = values
    return result


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Paralog-fraction and length-CV tracks for orthogroups")
    p.add_argument("input", help="orthologs TSV/TSV.GZ (Format A or B)")
    p.add_argument("--fasta", nargs="+", default=[], metavar="PATH",
                   help="protein FASTA files or directories (for length_cv)")
    p.add_argument("-o", "--output", help="output TSV path (default: stdout)")
    a = p.parse_args(argv)
    write_tracks(basic_tracks(a.input, a.fasta), COLUMNS, a.output)


if __name__ == "__main__":
    main()
