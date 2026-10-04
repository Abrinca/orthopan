from pathlib import Path

import pytest

from orthopan.tracks.orthofinder import COLUMNS, build_tracks, main

ORTHOGROUPS = (
    "Orthogroup\tA\tB\n"
    "OG1\tx\tx\n"
    "OG2\tx\tx\n"
    "OG3\tx\t\n"
    "OG4\tx\t\n"
    "OG5\t\tx\n"
)

DUPLICATIONS = (
    "Orthogroup\tSpecies Tree Node\tGene Tree Node\tSupport\tType\tGenes 1\tGenes 2\n"
    "OG1\tn1\tn1\t0.49\tTerminal\ta\tb\n"
    "OG1\tn2\tn2\t0.5\tTerminal\ta\tb\n"
    "OG1\tn3\tn3\t1.0\tTerminal\ta\tb\n"
    "OG2\tn1\tn1\t0.3076923076923077\tNon-Terminal\ta\tb\n"
)

GENE_TREES = (
    "OG1: (gnl|extdb|A.1_0001:0.5,gnl|extdb|A.1_0002:1.5e-2)n12:0.25;\n"
    "OG2: ((gnl|extdb|A.1_0001:1,gnl|extdb|A.1_0002:1)n12,gnl|extdb|B.1_0003:2);\n"
    "OG3: (gnl|extdb|A.1_0001:2.5,gnl|extdb|B.1_0003:2.5);\n"
)


def _results(tmp_path, *, duplications=True, trees=True):
    root = tmp_path / "Results"
    (root / "Orthogroups").mkdir(parents=True)
    (root / "Orthogroups" / "Orthogroups.tsv").write_text(ORTHOGROUPS)
    if duplications:
        (root / "Gene_Duplication_Events").mkdir()
        (root / "Gene_Duplication_Events" / "Duplications.tsv").write_text(DUPLICATIONS)
    if trees:
        (root / "Resolved_Gene_Trees").mkdir()
        (root / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text(GENE_TREES)
    return root


def test_duplication_counts_and_support_threshold(tmp_path):
    rows = build_tracks(_results(tmp_path))
    assert rows["OG1"] == {
        "duplications": 3.0,
        "duplications_supported": 2.0,  # 0.49 excluded, 0.5 included
        "tree_length": pytest.approx(0.5 + 0.015 + 0.25),
        "tree_length_per_leaf": pytest.approx((0.5 + 0.015 + 0.25) / 2),
    }
    assert rows["OG2"]["duplications"] == 1.0
    assert rows["OG2"]["duplications_supported"] == 0.0  # its only event is below 0.5


def test_absent_from_duplications_zero_or_empty(tmp_path):
    rows = build_tracks(_results(tmp_path))
    # Has a gene tree -> true zero.
    assert rows["OG3"] == {"duplications": 0.0, "duplications_supported": 0.0, **{
        "tree_length": 5.0, "tree_length_per_leaf": 2.5}}
    # No gene tree -> empty (unknown).
    assert rows["OG4"] == {}
    assert rows["OG5"] == {}


def test_tree_metrics_with_odd_labels(tmp_path):
    rows = build_tracks(_results(tmp_path))
    # Labels contain | and ., internal label n12 is :-free, branch is scientific.
    assert rows["OG1"]["tree_length"] == pytest.approx(0.765)
    assert rows["OG1"]["tree_length_per_leaf"] == pytest.approx(0.3825)
    # :-free internal label, 3 leaves, all integer branch lengths.
    assert rows["OG2"]["tree_length"] == pytest.approx(4.0)
    assert rows["OG2"]["tree_length_per_leaf"] == pytest.approx(4.0 / 3)


def test_missing_files_leave_columns_empty(tmp_path):
    rows = build_tracks(_results(tmp_path, duplications=False, trees=False))
    assert list(rows) == ["OG1", "OG2", "OG3", "OG4", "OG5"]  # one row per OG
    assert all(row == {} for row in rows.values())


def test_writes_header_and_rows(tmp_path):
    root = _results(tmp_path)
    out = tmp_path / "tracks.tsv"
    main([str(root), "-o", str(out)])
    lines = Path(out).read_text().splitlines()
    assert lines[0] == "\t".join(["orthogroup", *COLUMNS])
    assert lines[1].split("\t")[:3] == ["OG1", "3", "2"]
    assert lines[4].split("\t") == ["OG4", "", "", "", ""]


def test_branch_length_signed_and_leading_dot(tmp_path):
    root = tmp_path / "Results"
    (root / "Orthogroups").mkdir(parents=True)
    (root / "Orthogroups" / "Orthogroups.tsv").write_text("Orthogroup\tA\nOG1\tx\n")
    (root / "Resolved_Gene_Trees").mkdir()
    (root / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt").write_text("OG1: (a:.5,b:-0.1);\n")
    rows = build_tracks(root)
    assert rows["OG1"]["tree_length"] == pytest.approx(0.4)
    assert rows["OG1"]["tree_length_per_leaf"] == pytest.approx(0.2)


def test_short_duplication_row_skipped_with_warning(tmp_path, capsys):
    root = tmp_path / "Results"
    (root / "Orthogroups").mkdir(parents=True)
    (root / "Orthogroups" / "Orthogroups.tsv").write_text("Orthogroup\tA\nOG1\tx\nOG2\tx\n")
    (root / "Gene_Duplication_Events").mkdir()
    (root / "Gene_Duplication_Events" / "Duplications.tsv").write_text(
        "Orthogroup\tSpecies Tree Node\tGene Tree Node\tSupport\tType\tGenes 1\tGenes 2\n"
        "OG1\tn1\tn1\t1.0\tTerminal\ta\tb\n"
        "OG2\tN0\n"  # short row: padded Support cell is empty
    )
    rows = build_tracks(root)
    assert rows["OG1"]["duplications"] == 1.0
    assert rows["OG2"] == {}  # skipped, not a crash
    assert "skipped 1 incomplete rows" in capsys.readouterr().err


