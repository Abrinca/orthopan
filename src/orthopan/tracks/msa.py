"""Per-orthogroup tracks from OrthoFinder MSAs.

Reads one aligned FASTA per orthogroup (``<OG>.fa``) and writes four tracks:

``geometric_homogeneity``
    anvi'o geometric homogeneity in full mode (the arithmetic mean of the
    site-level ``by_residue`` and gene-level ``by_gene`` terms), computed from
    the closed form in ``docs/anvio-internals.md`` §1.  All-gap columns are
    dropped from the alignment first, because HOG pruning creates them and
    anvi'o would otherwise count them against every gene.  Needs two or more
    columns.
``functional_homogeneity``
    anvi'o functional homogeneity: amino-acid groups and 3/2/0 scoring, gaps
    skipped, a single ratio over the whole alignment (``docs/anvio-internals.md``
    §2).  ``-`` is the only gap (anvi'o semantics).  Works with a single column;
    columns with fewer than two non-gap sequences contribute nothing and the
    result is missing when no pair of non-gap residues exists anywhere.
``gap_fraction``
    Fraction of gap cells in the cleaned alignment (all-gap columns dropped).
``trimmed``
    1.0 if any confirmed-matched sequence's non-gap residue count in the MSA is
    smaller than its unaligned length in ``--seqs/<OG>.fa``; 0.0 only when every
    aligned sequence was matched and none is shorter; missing otherwise.
    Duplicate unaligned gene ids with conflicting lengths are unresolved.  Empty
    when ``--seqs`` is not given.

Unknown characters
------------------
anvi'o raises ``KeyError`` for residues outside its 24-letter table (``*``,
``U``, ``O``, lowercase, ``.``) and the caller then reports ``-1``.  This
module normalises instead: lowercase letters are upper-cased and every
character that is not one of anvi'o's table entries is mapped to ``X``, which
anvi'o scores 0 against everything (and 0 for ``X``/``X``).  A track value is
therefore never negative; uncomputable cases are written as an empty cell.

Paralog weighting
-----------------
anvi'o pools every gene in a cluster, so a genome with ``k`` paralogs
contributes ``k`` rows and ``O(k·n)`` pairs.  ``--weighting genome`` (the
default) instead gives each genome equal total weight: a genome with ``k``
copies contributes ``1/k`` per copy.  ``--weighting pooled`` reproduces anvi'o
exactly.

Genome assignment is exact, never by prefix.  ``--orthogroups`` is read with
``orthopan.io.read_matrix`` and streamed: each row's genes per genome column
form that OG's expected MSA headers, trying ``<genome>_<gene>`` and
OrthoFinder's normalised form with ``.`` replaced by ``_``.  A header that
matches more than one candidate (for example genome ``A`` with gene ``B_g``
against genome ``A_B`` with gene ``g``, or ``A.B`` normalising onto ``A_B``) or
none makes that OG fall back to pooled weighting, with one summary warning.
Membership is built one OG at a time and submitted with bounded in-flight work,
so the whole dataset is never held in memory; in pooled mode without a
``--seqs``-driven need it is not built at all.
"""

import argparse
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import (
    ALL_COMPLETED,
    FIRST_COMPLETED,
    ProcessPoolExecutor,
    wait,
)
from pathlib import Path

import numpy as np

from orthopan.io import read_matrix, split_genes
from orthopan.tracks import write_tracks

TRACKS = ["geometric_homogeneity", "functional_homogeneity", "gap_fraction", "trimmed"]

_GAP = 26  # code for '-'; letters A-Z are 0..25
_X = ord("X") - 65

# anvi'o conserved_amino_acid_groups (constants.py:438-481).  conserved(a) is
# the set of residues b for which is_amino_acid_functionally_conserved(a, b).
_CONSERVED = {
    "A": "LVIMCHA", "V": "LVIMCHA", "M": "LVIMCHA", "C": "LVIMCHA",
    "L": "LVIMCHAJ", "I": "LVIMCHAJ",
    "F": "FWY", "W": "FWY",
    "K": "KRH", "R": "KRH",
    "H": "HYLVIMCAKR",
    "Y": "HYFW",
    "S": "ST", "T": "ST",
    "N": "QNB", "Q": "QNZ", "D": "DEB", "E": "DEZ",
    "B": "BND", "Z": "ZQE", "J": "JLI",
}


def _score_matrix():
    """26x26 pair-score table: S[x, y] is anvi'o's score for an unordered pair."""
    s = np.zeros((26, 26), dtype=np.float64)
    for a, partners in _CONSERVED.items():
        ai = ord(a) - 65
        for b in partners:
            s[ai, ord(b) - 65] = 2.0
    for code in range(26):
        c = chr(65 + code)
        if c in "XJBZ":
            s[code, code] = 2.0 if c in "BZJ" else 0.0  # X self-conserves with nothing
        else:
            s[code, code] = 3.0
    return s


_SCORE = _score_matrix()
_DIAG = np.diag(_SCORE).copy()
_NONSTANDARD = np.full(256, _X, dtype=np.uint8)
for _c in range(65, 91):  # A-Z
    if chr(_c) not in "OU":
        _NONSTANDARD[_c] = _c - 65
        _NONSTANDARD[_c + 32] = _c - 65  # lowercase
_NONSTANDARD[ord("-")] = _GAP


def read_fasta(path):
    """Return ``[(header, sequence), ...]`` with every line stripped.

    Stripping handles OrthoFinder's trailing space and 80-column wrapping.
    """
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def msa_metrics(seqs, weights=None):
    """Tracks for one alignment.

    ``seqs`` are equal-length aligned strings; ``weights`` is a per-sequence
    weight vector (``None`` means anvi'o's pooled, all-ones weighting).
    """
    n = len(seqs)
    if n == 0 or any(len(s) != len(seqs[0]) for s in seqs):
        return {t: None for t in TRACKS}
    codes = _NONSTANDARD[np.frombuffer("".join(seqs).encode("latin-1", "replace"), np.uint8)]
    codes = codes.reshape(n, -1)
    gap = codes == _GAP
    keep = ~gap.all(0)  # drop all-gap columns
    codes, gap = codes[:, keep], gap[:, keep]
    L = codes.shape[1]

    result = {t: None for t in TRACKS}
    if L == 0 or n < 1:
        result["gap_fraction"] = 0.0
        return result
    result["gap_fraction"] = float(gap.mean())
    if n < 2:
        return result

    wr = np.ones(n, dtype=np.float64) if weights is None else np.asarray(weights, dtype=np.float64)
    G = float(wr.sum())
    Sw2 = float((wr * wr).sum())

    if L >= 2:
        k = gap.sum(1).astype(np.float64)
        by_residue = float((wr * (k * (k - 1) + (L - k) * (L - k - 1))).sum()) / (G * L * (L - 1))
        W = gap.T @ wr  # weighted gap count per column
        by_gene = float((W * W + (G - W) * (G - W) - Sw2).sum()) / (L * (G * G - Sw2))
        result["geometric_homogeneity"] = (by_residue + by_gene) / 2.0

    wr2 = wr * wr
    T = np.zeros((L, 26), dtype=np.float64)
    Uc = np.zeros((L, 26), dtype=np.float64)
    for x in np.unique(codes):
        if x == _GAP:
            continue
        m = codes == x
        T[:, x] = m.T @ wr
        Uc[:, x] = m.T @ wr2
    score = 0.5 * (np.einsum("lx,xy,ly->l", T, _SCORE, T) - (Uc * _DIAG).sum(1))

    # Non-gap totals straight from the residue counts.  Subtracting weighted gap
    # totals cancels catastrophically for many small weights (paralogs), so a
    # full column of one-per-genome fragments used to yield a negative maximum.
    weighted_non_gap = T.sum(1)
    squared_non_gap = Uc.sum(1)
    eligible = (~gap).sum(0) >= 2  # a pair needs two non-gap sequences
    score = np.where(eligible, score, 0.0)
    max_score = np.where(eligible, 1.5 * (weighted_non_gap ** 2 - squared_non_gap), 0.0)
    total = float(max_score.sum())
    if total > 0:
        # Clamp only as a guard against float noise; the maths is already in [0, 1].
        result["functional_homogeneity"] = min(1.0, max(0.0, float(score.sum()) / total))
    return result


def _expected_headers(members):
    """{MSA header: {(genome, gene)}} for one OG's Orthogroups.tsv membership."""
    expected = defaultdict(set)
    for genome, genes in (members or {}).items():
        normalised = genome.replace(".", "_")
        for gene in genes:
            expected[f"{genome}_{gene}"].add((genome, gene))
            expected[f"{normalised}_{gene}"].add((genome, gene))
    return expected


def _match_headers(headers, expected):
    """Return ({header: (genome, gene) | None}, n_ambiguous, n_unmatched)."""
    matched, n_ambiguous, n_unmatched = {}, 0, 0
    for header in headers:
        candidates = expected.get(header)
        if candidates is None:
            n_unmatched += 1
            matched[header] = None
        elif len(candidates) == 1:
            matched[header] = next(iter(candidates))
        else:
            n_ambiguous += 1
            matched[header] = None
    return matched, n_ambiguous, n_unmatched


def _unaligned_lengths(path):
    """{gene id: unaligned length}, with conflicting duplicate ids set to None."""
    lengths = {}
    for name, seq in read_fasta(path):
        length = len(seq.rstrip("*"))
        if name not in lengths:
            lengths[name] = length
        elif lengths[name] != length:
            lengths[name] = None
    return lengths


def _unique_suffix_gene(header, unaligned):
    """Unambiguously suffix-match a header against the unaligned gene ids."""
    candidates = {header} if header in unaligned else set()
    for i, char in enumerate(header):
        if char == "_" and header[i + 1:] in unaligned:
            candidates.add(header[i + 1:])
    if len(candidates) != 1:
        return None
    gene = candidates.pop()
    return gene if unaligned[gene] is not None else None


def _trimmed(headers, seqs, matched, unaligned):
    any_shorter = False
    all_matched = True
    for header, seq in zip(headers, seqs):
        match = matched.get(header)
        gene = match[1] if match else None
        if gene is None or gene not in unaligned:
            gene = _unique_suffix_gene(header, unaligned)
        if gene is None or unaligned.get(gene) is None:
            all_matched = False
            continue
        if len(seq) - seq.count("-") < unaligned[gene]:
            any_shorter = True
    if any_shorter:
        return 1.0
    return 0.0 if all_matched else None


_CFG = {}


def _init(cfg):
    _CFG.update(cfg)


def _process(og, members):
    cfg = _CFG
    result = {t: None for t in TRACKS}
    msa_path = Path(cfg["msa_dir"]) / f"{og}.fa"
    if not msa_path.is_file():
        return og, result, False
    pairs = read_fasta(msa_path)
    if not pairs:
        return og, result, False
    headers = [h for h, _ in pairs]
    seqs = [s for _, s in pairs]

    expected = _expected_headers(members)
    matched, n_ambiguous, n_unmatched = _match_headers(headers, expected)
    weights = None
    fallback = False
    if cfg["weighting"] == "genome":
        if expected and not n_ambiguous and not n_unmatched:
            counts = Counter(matched[h][0] for h in headers)
            weights = [1.0 / counts[matched[h][0]] for h in headers]
        else:
            fallback = True

    result.update(msa_metrics(seqs, weights))

    seqs_dir = cfg["seqs_dir"]
    if seqs_dir:
        seq_path = Path(seqs_dir) / f"{og}.fa"
        if seq_path.is_file():
            result["trimmed"] = _trimmed(headers, seqs, matched, _unaligned_lengths(seq_path))
    return og, result, fallback


def _process_batch(batch):
    return [_process(og, members) for og, members in batch]


def _members_from_cells(labels, cells):
    return {labels[i]: split_genes(cell) for i, cell in enumerate(cells) if cell}


def build_tracks(msa_dir, seqs=None, orthogroups=None, threads=None, weighting="genome"):
    """Compute the MSA tracks for every orthogroup; return ``{og: {track: value}}``."""
    msa_dir = Path(msa_dir)
    if not msa_dir.is_dir():
        raise ValueError(f"{msa_dir} is not a directory")
    if threads is None:
        threads = os.cpu_count() or 1
    msa_ogs = {f.stem for f in msa_dir.glob("*.fa")}

    if weighting == "genome" and not orthogroups:
        print("warning: no genome mapping (--orthogroups not given); falling back to pooled", file=sys.stderr)
        weighting = "pooled"
    # Membership is only needed to weight genomes or to exactly match headers for trimmed.
    need_members = weighting == "genome" or seqs is not None

    cfg = {
        "msa_dir": str(msa_dir),
        "seqs_dir": str(seqs) if seqs else None,
        "weighting": weighting,
    }

    t0 = time.time()
    results = {}
    n_fallback = 0

    if threads and threads > 1 and (msa_ogs or orthogroups):
        # Platform default start method (forkserver/spawn): fork after native threads
        # can deadlock. _init passes config explicitly; batches bound in-flight memory.
        with ProcessPoolExecutor(
            max_workers=threads, initializer=_init, initargs=(cfg,)
        ) as ex:
            pending, processed, batch = [], set(), []
            limit, chunk = 2 * threads, 32

            def record(item):
                nonlocal n_fallback
                og, values, fallback = item
                results[og] = values
                n_fallback += fallback

            def flush(*, force=False):
                if not batch:
                    return
                pending.append(ex.submit(_process_batch, list(batch)))
                batch.clear()
                if force or len(pending) >= limit:
                    done, still = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        for item in future.result():
                            record(item)
                    pending[:] = list(still)

            if orthogroups:
                labels, rows = read_matrix(orthogroups)
                for og, cells in rows:
                    processed.add(og)
                    batch.append((og, _members_from_cells(labels, cells) if need_members else {}))
                    if len(batch) >= chunk:
                        flush()
            for og in sorted(msa_ogs - processed):
                batch.append((og, {}))
                if len(batch) >= chunk:
                    flush()
            flush(force=True)
            done, _ = wait(pending, return_when=ALL_COMPLETED)
            for future in done:
                for item in future.result():
                    record(item)
    else:
        _init(cfg)
        processed = set()
        if orthogroups:
            labels, rows = read_matrix(orthogroups)
            for og, cells in rows:
                processed.add(og)
                members = _members_from_cells(labels, cells) if need_members else {}
                _, results[og], fallback = _process(og, members)
                n_fallback += fallback
        for og in sorted(msa_ogs - processed):
            _, results[og], fallback = _process(og, {})
            n_fallback += fallback

    if n_fallback:
        print(
            f"warning: {n_fallback} orthogroups fell back to pooled weighting "
            "(ambiguous or unmatched MSA headers)",
            file=sys.stderr,
        )

    elapsed = time.time() - t0
    n_msa = len(msa_ogs)
    print(
        f"{len(results)} orthogroups ({n_msa} with MSAs) in {elapsed:.2f}s "
        f"({n_msa / elapsed:.0f} MSAs/s, {threads} threads)",
        file=sys.stderr,
    )
    return results


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("msa_dir", help="MultipleSequenceAlignments directory")
    p.add_argument("--seqs", help="Orthogroup_Sequences directory (for the trimmed flag)")
    p.add_argument("-o", "--out", help="output track TSV (default: stdout)")
    p.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    p.add_argument("--weighting", choices=["genome", "pooled"], default="genome")
    p.add_argument("--orthogroups", help="Orthogroups.tsv (genome mapping for --weighting genome)")
    args = p.parse_args(argv)
    if not Path(args.msa_dir).is_dir():
        p.error(f"{args.msa_dir} is not a directory")
    results = build_tracks(args.msa_dir, args.seqs, args.orthogroups, args.threads, args.weighting)
    write_tracks({og: results[og] for og in sorted(results)}, TRACKS, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
