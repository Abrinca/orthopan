"""Tests for orthopan.tracks.msa."""

import importlib.util
import sys
import types
from pathlib import Path

import numpy as np
import pytest

from orthopan.tracks import msa as M

# ---------------------------------------------------------------------------
# anvi'o reference: load homogeneityindex.py with lightweight stubs and a naive
# transcription of is_amino_acid_functionally_conserved for the functional index.
# ---------------------------------------------------------------------------

_GROUPS = {
    "Nonpolar": "LVIMCHA",
    "Aromatic": "FWY",
    "Bases": "KRH",
    "Neutral Amines": "QN",
    "Acids": "DE",
    "Polar and Nonpolar": "HY",
    "Mostly nonpolar": "ST",
    "B": "BND",
    "Z": "ZQE",
    "J": "JLI",
    "None": "",
}
_GROUPS["N"] = _GROUPS["Neutral Amines"] + "B"
_GROUPS["D"] = _GROUPS["Acids"] + "B"
_GROUPS["Q"] = _GROUPS["Neutral Amines"] + "Z"
_GROUPS["E"] = _GROUPS["Acids"] + "Z"
_GROUPS["LI"] = _GROUPS["Nonpolar"] + "J"

_PROP = {}
for _k in "AVMC":
    _PROP[_k] = "Nonpolar"
for _k in "FW":
    _PROP[_k] = "Aromatic"
for _k in "KR":
    _PROP[_k] = "Bases"
for _k in "HY":
    _PROP[_k] = "Polar and Nonpolar"
for _k in "ST":
    _PROP[_k] = "Mostly nonpolar"
for _k in "GPX":
    _PROP[_k] = "None"
for _k in "BZJNDQE":
    _PROP[_k] = _k
for _k in "LI":
    _PROP[_k] = "LI"


def _conserved(a, b):
    group = _PROP[a]
    if b in _GROUPS[group]:
        return True
    if group == "Polar and Nonpolar":
        if a == "H" and (b in _GROUPS["Nonpolar"] or b in _GROUPS["Bases"]):
            return True
        if a == "Y" and b in _GROUPS["Aromatic"]:
            return True
    return False


def _load_anvio():
    an = types.ModuleType("anvio")
    an.__version__ = "test"
    term = types.ModuleType("anvio.terminal")
    term.Run = term.Progress = lambda *a, **k: None
    ut = types.ModuleType("anvio.utils")
    ut.is_amino_acid_functionally_conserved = _conserved
    sys.modules.update({"anvio": an, "anvio.terminal": term, "anvio.utils": ut})
    an.terminal, an.utils = term, ut
    path = Path(__file__).resolve().parents[1] / "resources/anvio/anvio/homogeneityindex.py"
    spec = importlib.util.spec_from_file_location("anvio_homogeneityindex", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.HomogeneityCalculator()


_ANVIO = _load_anvio()


# ---------------------------------------------------------------------------
# pooled equals anvi'o
# ---------------------------------------------------------------------------

HANDMADE = [
    ["AAAA", "AA--"],
    ["ACDEFG", "ACDEFG"],
    ["ACDEFG", "ACDE--"],
    ["LIVMCH", "LIVMCJ"],
    ["BNDQEZ", "BNDQEZ"],
    ["XXAA--", "XXAABB"],
    ["QQQNNN", "QQQN--"],
]


@pytest.mark.parametrize("seqs", HANDMADE)
def test_pooled_geometric_matches_anvio(seqs):
    ours = M.msa_metrics(seqs, None)["geometric_homogeneity"]
    assert ours == pytest.approx(_ANVIO.compute_geometric_index(seqs), abs=1e-12)


@pytest.mark.parametrize("seqs", HANDMADE)
def test_pooled_functional_matches_anvio(seqs):
    ours = M.msa_metrics(seqs, None)["functional_homogeneity"]
    assert ours == pytest.approx(_ANVIO.compute_functional_index(seqs), abs=1e-12)


def test_pooled_matches_anvio_random():
    rng = np.random.default_rng(0)
    alphabet = "ACDEFGHIKLMNPQRSTVWYBZJX--"
    for _ in range(150):
        n = int(rng.integers(2, 8))
        L = int(rng.integers(2, 25))
        seqs = ["".join(rng.choice(list(alphabet), L)) for _ in range(n)]
        if any(all(s[c] == "-" for s in seqs) for c in range(L)):
            continue  # anvi'o keeps all-gap columns; ours drops them
        ours = M.msa_metrics(seqs, None)
        assert ours["geometric_homogeneity"] == pytest.approx(
            _ANVIO.compute_geometric_index(seqs), abs=1e-12
        )
        assert ours["functional_homogeneity"] == pytest.approx(
            _ANVIO.compute_functional_index(seqs), abs=1e-12
        )


# ---------------------------------------------------------------------------
# all-gap columns, single sequence, unknown characters, one column
# ---------------------------------------------------------------------------


def test_all_gap_columns_dropped():
    base = M.msa_metrics(["AAAA", "AA--"], None)
    with_gap = M.msa_metrics(["AAAA-", "AA---"], None)  # last column is all gaps
    assert with_gap == pytest.approx(base)


def test_single_sequence_empty():
    m = M.msa_metrics(["AAAA"], None)
    assert m["geometric_homogeneity"] is None
    assert m["functional_homogeneity"] is None
    assert m["gap_fraction"] == 0.0


def test_unknown_characters_never_negative():
    for seqs in (["AAAA", "AA*U"], ["aaaa", "AAAA"], ["M*UO", "M*UO"], ["A.B.", "A-B-"]):
        m = M.msa_metrics(seqs, None)
        for key in ("geometric_homogeneity", "functional_homogeneity"):
            assert m[key] is None or 0.0 <= m[key] <= 1.0
    # lowercase is upper-cased, so a/a is an identity, not a KeyError
    assert M.msa_metrics(["aaaa", "AAAA"], None)["functional_homogeneity"] == pytest.approx(1.0)


def test_one_column_identical():
    m = M.msa_metrics(["A", "A"], None)
    assert m["geometric_homogeneity"] is None  # needs two columns
    assert m["functional_homogeneity"] == pytest.approx(1.0)
    assert m["gap_fraction"] == 0.0


def test_one_column_conserved():
    # A is "conserved" with V (both Nonpolar) -> 2 out of 3.
    assert M.msa_metrics(["A", "V"], None)["functional_homogeneity"] == pytest.approx(2 / 3)


def test_one_column_non_overlapping_is_missing():
    assert M.msa_metrics(["A", "-"], None)["functional_homogeneity"] is None
    assert M.msa_metrics(["A-", "-A"], None)["functional_homogeneity"] is None


# ---------------------------------------------------------------------------
# paralog weighting
# ---------------------------------------------------------------------------


def test_genome_weighting_differs_with_paralogs():
    # Genome A has two identical truncated paralogs; genome B one full copy.
    # Pooling lets A's two copies dominate and pulls the score down; genome
    # weighting makes A count once.
    seqs = ["AA--", "AA--", "AAAA"]
    pooled = M.msa_metrics(seqs, None)["geometric_homogeneity"]
    weighted = M.msa_metrics(seqs, [0.5, 0.5, 1.0])["geometric_homogeneity"]
    assert pooled == pytest.approx(0.6111111111111112)
    assert weighted > pooled
    assert weighted == pytest.approx(0.6333333333333333)


def test_pooled_equals_unit_weights():
    seqs = ["ACDEFG", "ACDE--", "ACD---"]
    assert M.msa_metrics(seqs, None) == pytest.approx(M.msa_metrics(seqs, [1.0, 1.0, 1.0]))


def test_weighted_non_overlapping_fragments_missing_not_negative():
    # 14 paralogs, each with a single residue in its own column, weighted 1/14.
    # Subtracting weighted gap totals used to cancel and give -0.3889.
    n = 14
    seqs = ["".join("-" * i + "A" + "-" * (n - 1 - i)) for i in range(n)]
    w = [1.0 / n] * n
    for weights in (w, None):
        m = M.msa_metrics(seqs, weights)
        assert m["functional_homogeneity"] is None


def test_weighted_shared_column_is_one():
    n = 14
    seqs = ["A" + "".join("-" * i + "A" + "-" * (n - 1 - i)) for i in range(n)]
    m = M.msa_metrics(seqs, [1.0 / n] * n)
    assert m["functional_homogeneity"] == pytest.approx(1.0)  # clamped, shared column only


# ---------------------------------------------------------------------------
# end-to-end helpers
# ---------------------------------------------------------------------------


def _write_fasta(path, entries):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in entries))


def _run(msa_dir, *, seqs=None, orthogroups=None, weighting="genome", threads=1):
    out = Path(msa_dir).parent / "out.tsv"
    args = [str(msa_dir)]
    if seqs:
        args += ["--seqs", str(seqs)]
    if orthogroups:
        args += ["--orthogroups", str(orthogroups)]
    args += ["--weighting", weighting, "-o", str(out), "--threads", str(threads)]
    M.main(args)
    rows = {}
    for line in out.read_text().splitlines()[1:]:
        cells = line.split("\t")
        rows[cells[0]] = [float(c) if c != "" else None for c in cells[1:]]
    return rows


def _geo(row):
    return row[0]


# ---------------------------------------------------------------------------
# trimmed detection
# ---------------------------------------------------------------------------


def test_trimmed_detection(tmp_path):
    msa_dir = tmp_path / "msa"
    seqs_dir = tmp_path / "seqs"
    msa_dir.mkdir()
    seqs_dir.mkdir()
    # OG1: first sequence is 5 residues in the MSA but 10 unaligned -> trimmed.
    _write_fasta(seqs_dir / "OG1.fa", [("gnl|g1", "AAAAAAAAAA"), ("gnl|g2", "CCCCCCCCCC")])
    _write_fasta(
        msa_dir / "OG1.fa",
        [("FAM1_1_gnl|g1", "AAAAA-----"), ("FAM1_1_gnl|g2", "CCCCCCCCCC")],
    )
    # OG2: same lengths -> not trimmed; trailing space must be ignored.
    _write_fasta(seqs_dir / "OG2.fa", [("gnl|g1", "AAAAAAAAAA")])
    (msa_dir / "OG2.fa").write_text(">FAM1_1_gnl|g1\nAAAAAAAAAA ")

    rows = _run(msa_dir, seqs=seqs_dir)
    assert rows["OG1"][3] == 1.0
    assert rows["OG2"][3] == 0.0


def test_trimmed_unique_suffix(tmp_path):
    msa_dir = tmp_path / "msa"
    seqs_dir = tmp_path / "seqs"
    msa_dir.mkdir()
    seqs_dir.mkdir()
    _write_fasta(seqs_dir / "OG1.fa", [("g1", "AAAA")])
    _write_fasta(msa_dir / "OG1.fa", [("prefix_g1", "AAAA")])
    _write_fasta(seqs_dir / "OG2.fa", [("g1", "AAAA")])
    _write_fasta(msa_dir / "OG2.fa", [("prefix_g1", "A---")])  # 1 residue < 4
    rows = _run(msa_dir, seqs=seqs_dir)
    assert rows["OG1"][3] == 0.0
    assert rows["OG2"][3] == 1.0


def test_trimmed_conflicting_duplicate_ids(tmp_path):
    msa_dir = tmp_path / "msa"
    seqs_dir = tmp_path / "seqs"
    msa_dir.mkdir()
    seqs_dir.mkdir()
    # Same id "g1" appears with two different lengths -> unresolved.
    _write_fasta(seqs_dir / "OG1.fa", [("g1", "A" * 10), ("g1", "A" * 20)])
    _write_fasta(msa_dir / "OG1.fa", [("g1", "A" * 10 + "-" * 10), ("g1", "A" * 20)])
    rows = _run(msa_dir, seqs=seqs_dir)
    assert rows["OG1"][3] is None  # unresolved, not a false "trimmed"
    # A third, unambiguous sequence that is shorter proves trimming anyway.
    _write_fasta(
        seqs_dir / "OG2.fa", [("g1", "A" * 10), ("g1", "A" * 20), ("g2", "A" * 10)]
    )
    _write_fasta(
        msa_dir / "OG2.fa",
        [("g1", "A" * 10 + "-" * 10), ("g1", "A" * 20), ("g2", "A" * 5 + "-" * 15)],
    )
    rows = _run(msa_dir, seqs=seqs_dir)
    assert rows["OG2"][3] == 1.0


def test_trimmed_unresolved_id_not_replaced_by_other_suffix(tmp_path):
    msa_dir = tmp_path / "msa"
    seqs_dir = tmp_path / "seqs"
    msa_dir.mkdir()
    seqs_dir.mkdir()
    # "B_g" is unresolved (two lengths); the known "g" must not stand in for it.
    _write_fasta(seqs_dir / "OG1.fa", [("B_g", "AA"), ("B_g", "AAAA"), ("g", "A" * 10)])
    _write_fasta(msa_dir / "OG1.fa", [("A_B_g", "AA")])
    rows = _run(msa_dir, seqs=seqs_dir)
    assert rows["OG1"][3] is None


def test_trimmed_empty_and_unrelated_and_ambiguous(tmp_path):
    msa_dir = tmp_path / "msa"
    seqs_dir = tmp_path / "seqs"
    msa_dir.mkdir()
    seqs_dir.mkdir()
    _write_fasta(msa_dir / "EMPTY.fa", [("h1", "AAAA")])
    (seqs_dir / "EMPTY.fa").write_text("")  # empty FASTA
    _write_fasta(msa_dir / "UNREL.fa", [("h1", "AAAA")])
    _write_fasta(seqs_dir / "UNREL.fa", [("otherXYZ", "AAAA")])  # no match
    _write_fasta(msa_dir / "AMBIG.fa", [("x_a_g", "AAAA")])
    _write_fasta(seqs_dir / "AMBIG.fa", [("a_g", "AAAA"), ("g", "AAAA")])  # two suffixes
    rows = _run(msa_dir, seqs=seqs_dir)
    assert rows["EMPTY"][3] is None
    assert rows["UNREL"][3] is None
    assert rows["AMBIG"][3] is None


# ---------------------------------------------------------------------------
# genome mapping / weighting
# ---------------------------------------------------------------------------

PARALOG_SEQS = ["AA--", "AA--", "AAAA"]
PARALOG_GEO = 0.6333333333333333
POOLED_GEO = 0.6111111111111112


def test_dotted_genome_weighting_is_numeric(tmp_path):
    # Genome "FAM.a" is written normalised as FAM_a in the MSA headers; two of
    # its copies are truncated paralogs, genome FAM.b has one full copy.
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    _write_fasta(
        msa_dir / "OG1.fa",
        [("FAM_a_g1", "AA--"), ("FAM_a_g2", "AA--"), ("FAM_b_g3", "AAAA")],
    )
    og = tmp_path / "Orthogroups.tsv"
    og.write_text(
        "Orthogroup\tFAM.a\tFAM.b\n"
        "OG1\tg1, g2\tg3\n"
    )
    rows = _run(msa_dir, orthogroups=og, weighting="genome")
    assert _geo(rows["OG1"]) == pytest.approx(PARALOG_GEO)
    assert _geo(rows["OG1"]) == pytest.approx(
        M.msa_metrics(PARALOG_SEQS, [0.5, 0.5, 1.0])["geometric_homogeneity"]
    )
    pooled = _run(msa_dir, orthogroups=og, weighting="pooled")
    assert _geo(pooled["OG1"]) == pytest.approx(POOLED_GEO)


def test_genome_names_with_dash_dot_underscore(tmp_path):
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    _write_fasta(
        msa_dir / "OG1.fa",
        [("FAM_a-b_c_g1", "AA--"), ("FAM_a-b_c_g2", "AA--"), ("FAM_b_g3", "AAAA")],
    )
    og = tmp_path / "Orthogroups.tsv"
    og.write_text(
        "Orthogroup\tFAM.a-b_c\tFAM.b\n"
        "OG1\tg1, g2\tg3\n"
    )
    rows = _run(msa_dir, orthogroups=og, weighting="genome")
    assert _geo(rows["OG1"]) == pytest.approx(PARALOG_GEO)


def test_prefix_overlap_fam1_fam10(tmp_path):
    # FAM10 must not be confused with FAM1; FAM1 has a paralog, FAM10 does not.
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    _write_fasta(
        msa_dir / "OG1.fa",
        [("FAM1_g1", "AA--"), ("FAM1_g2", "AA--"), ("FAM10_g1", "AAAA")],
    )
    og = tmp_path / "Orthogroups.tsv"
    og.write_text(
        "Orthogroup\tFAM1\tFAM10\n"
        "OG1\tg1, g2\tg1\n"
    )
    rows = _run(msa_dir, orthogroups=og, weighting="genome")
    assert _geo(rows["OG1"]) == pytest.approx(PARALOG_GEO)
    pooled = _run(msa_dir, orthogroups=og, weighting="pooled")
    assert _geo(pooled["OG1"]) == pytest.approx(POOLED_GEO)


def test_ambiguous_a_and_a_b_falls_back(tmp_path, capsys):
    # Header "A_B_g" is genome A with gene B_g, or genome A_B with gene g.
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    _write_fasta(msa_dir / "OG1.fa", [("A_B_g", "AA--"), ("A_B_x", "AA--")])
    og = tmp_path / "Orthogroups.tsv"
    og.write_text(
        "Orthogroup\tA\tA_B\n"
        "OG1\tB_g, B_x\tg, x\n"
    )
    genome = _run(msa_dir, orthogroups=og, weighting="genome")
    genome_err = capsys.readouterr().err
    pooled = _run(msa_dir, orthogroups=og, weighting="pooled")
    assert genome == pooled
    assert "fell back to pooled" in genome_err


def test_normalisation_collision_falls_back(tmp_path, capsys):
    # A.B and A_B both normalise to A_B, so the header is ambiguous.
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    _write_fasta(msa_dir / "OG1.fa", [("A_B_g", "AA--"), ("A_B_h", "AA--")])
    og = tmp_path / "Orthogroups.tsv"
    og.write_text(
        "Orthogroup\tA.B\tA_B\n"
        "OG1\tg, h\tg, h\n"
    )
    genome = _run(msa_dir, orthogroups=og, weighting="genome")
    genome_err = capsys.readouterr().err
    pooled = _run(msa_dir, orthogroups=og, weighting="pooled")
    assert genome == pooled
    assert "fell back to pooled" in genome_err


def test_headerless_orthogroups_falls_back(tmp_path, capsys):
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    _write_fasta(msa_dir / "OG1.fa", [("FAM_a_g1", "AA--"), ("FAM_a_g2", "AA--"), ("FAM_b_g3", "AAAA")])
    og = tmp_path / "Orthogroups.tsv"
    # No "Orthogroup" header -> read_matrix labels genomes genome_1, genome_2.
    og.write_text("OG1\tg1, g2\tg3\n")
    genome = _run(msa_dir, orthogroups=og, weighting="genome")
    genome_err = capsys.readouterr().err
    pooled = _run(msa_dir, orthogroups=og, weighting="pooled")
    assert _geo(genome["OG1"]) == pytest.approx(_geo(pooled["OG1"]))
    assert "fell back to pooled" in genome_err


def test_orthogroups_extends_universe(tmp_path):
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    _write_fasta(msa_dir / "OG1.fa", [("FAM_a_g1", "AAAA")])
    og = tmp_path / "Orthogroups.tsv"
    og.write_text("Orthogroup\tFAM.a\nOG1\tg1\nOG2\tg4\n")
    rows = _run(msa_dir, orthogroups=og, weighting="pooled")
    assert set(rows) == {"OG1", "OG2"}
    assert rows["OG2"] == [None, None, None, None]


def test_parallel_matches_sequential(tmp_path):
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    for i in range(40):
        _write_fasta(msa_dir / f"OG{i:03d}.fa", [("FAM_a_g1", "AA--"), ("FAM_b_g1", "AAAA")])
    og = tmp_path / "Orthogroups.tsv"
    og.write_text(
        "Orthogroup\tFAM.a\tFAM.b\n" + "".join(f"OG{i:03d}\tg1\tg1\n" for i in range(40))
    )
    sequential = _run(msa_dir, orthogroups=og, weighting="pooled", threads=1)
    parallel = _run(msa_dir, orthogroups=og, weighting="pooled", threads=2)
    assert parallel == sequential


# ---------------------------------------------------------------------------
# CLI errors
# ---------------------------------------------------------------------------


def test_missing_msa_dir_errors(tmp_path):
    with pytest.raises(SystemExit) as exc:
        M.main([str(tmp_path / "nope"), "--threads", "1"])
    assert exc.value.code != 0


def test_empty_msa_dir_is_fine(tmp_path):
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    og = tmp_path / "Orthogroups.tsv"
    og.write_text("Orthogroup\tFAM.a\nOG1\tg1\n")
    rows = _run(msa_dir, orthogroups=og, weighting="pooled")
    assert rows == {"OG1": [None, None, None, None]}
