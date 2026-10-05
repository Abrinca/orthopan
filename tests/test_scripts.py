"""Tests for the maintenance scripts under `scripts/`."""

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(
        f"_script_{name}", REPO / "scripts" / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


make_test_metadata = _load("make_test_metadata")


def test_read_genomes_headerless_format_b(tmp_path):
    # The first line is data, so the genomes are head[1:]: labels are
    # genome_1..genome_{len(head)-1}, not one extra.
    path = tmp_path / "matrix.tsv"
    path.write_text("OG1\tx\ty\tz\nOG2\t\tx\t\n")
    labels, ogs = make_test_metadata.read_genomes_and_ogs(path)
    assert labels == ["genome_1", "genome_2", "genome_3"]
    assert ogs == ["OG1", "OG2"]


def test_read_genomes_header_format_b(tmp_path):
    path = tmp_path / "matrix.tsv"
    path.write_text("Orthogroup\tA\tB\nOG1\tx\t\n")
    assert make_test_metadata.read_genomes_and_ogs(path) == (["A", "B"], ["OG1"])


def test_read_genomes_format_a(tmp_path):
    path = tmp_path / "matrix.tsv"
    path.write_text("OG1\tgnl|x|A_1, gnl|x|B_2\n")
    labels, ogs = make_test_metadata.read_genomes_and_ogs(path)
    assert labels == ["A", "B"]
    assert ogs == ["OG1"]
