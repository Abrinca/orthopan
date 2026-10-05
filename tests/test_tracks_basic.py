import gzip
import statistics

import pytest

from orthopan.tracks.basic import basic_tracks, load_lengths

MATRIX = (
    "Orthogroup\tA\tB\tC\n"
    "OG1\tg1\tg1, g2\tg1\n"
    "OG2\tg1\t\t\n"
    "OG3\t\t\t\n"
)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_paralog_fraction(tmp_path):
    rows = basic_tracks(_write(tmp_path, "og.tsv", MATRIX))
    og1, og2, og3 = rows["OG1"], rows["OG2"], rows["OG3"]
    assert og1["paralog_fraction"] == pytest.approx(1 / 3)
    assert og2["paralog_fraction"] == 0
    assert og3 == {"paralog_fraction": None, "length_cv": None}


def test_paralog_fraction_gnl_prefix(tmp_path):
    # Real OrthoFinder IDs carry a `gnl|extdb|` prefix; cell count is by commas only.
    text = "Orthogroup\tA\tB\nOG1\tgnl|extdb|x_1, gnl|extdb|x_2\tgnl|extdb|y_1\n"
    row = basic_tracks(_write(tmp_path, "p.tsv", text))["OG1"]
    assert row["paralog_fraction"] == pytest.approx(0.5)


def test_length_cv_filled_and_gaps(tmp_path):
    fasta = _write(tmp_path, "p.faa", ">g1\nAAAAAAAAAA\n>g2\n" + "A" * 20 + "\n>g3\n" + "A" * 30 + "\n")
    text = (
        "Orthogroup\tA\tB\n"
        "OG1\tg1\tg2\n"
        "OG2\tg1\t\n"
        "OG3\tg1, g2, g3\tg1\n"
        "OG4\tghost\tg1\n"
    )
    rows = basic_tracks(_write(tmp_path, "og.tsv", text), [fasta])
    assert rows["OG1"]["length_cv"] == pytest.approx(statistics.pstdev([10, 20]) / 15)
    assert rows["OG2"]["length_cv"] is None  # <2 genes
    assert rows["OG4"]["length_cv"] is None  # only one gene found in the FASTA
    assert rows["OG3"]["length_cv"] == pytest.approx(statistics.pstdev([10, 20, 30, 10]) / 17.5)


def test_length_cv_absent_without_fasta(tmp_path):
    rows = basic_tracks(_write(tmp_path, "og.tsv", MATRIX))
    assert all(r["length_cv"] is None for r in rows.values())


def test_load_lengths_strips_stop_and_whitespace(tmp_path):
    fasta = _write(tmp_path, "s.faa", ">gx description here\nACDEF*\nGH\n")
    assert load_lengths([fasta]) == {"gx": 7}


def test_load_lengths_from_directory(tmp_path):
    d = tmp_path / "faa"
    d.mkdir()
    (d / "a.faa").write_text(">g1\nAAA\n")
    with gzip.open(d / "b.fa.gz", "wt") as fh:
        fh.write(">g2\n" + "A" * 5 + "\n")
    (d / "notes.txt").write_text("ignored\n")
    assert load_lengths([d]) == {"g1": 3, "g2": 5}


def test_short_rows_padded(tmp_path):
    # Width comes from the header; a short row is just fewer present genomes.
    text = "Orthogroup\tA\tB\tC\nOG1\tg1\nOG2\t\tg2\n"
    rows = basic_tracks(_write(tmp_path, "s.tsv", text))
    assert rows["OG1"] == {"paralog_fraction": 0, "length_cv": None}
    assert rows["OG2"]["paralog_fraction"] == 0


def test_headerless(tmp_path):
    text = "OG1\tg1, g2\tg1\nOG2\t\tg1\n"
    rows = basic_tracks(_write(tmp_path, "nh.tsv", text))
    assert rows["OG1"]["paralog_fraction"] == pytest.approx(0.5)
    assert rows["OG2"]["paralog_fraction"] == 0


def test_gz_input_matches_plain(tmp_path):
    plain = basic_tracks(_write(tmp_path, "og.tsv", MATRIX))
    gz = tmp_path / "og.tsv.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(MATRIX)
    assert basic_tracks(gz) == plain


def test_repeated_gene_id_across_proteomes(tmp_path):
    # Same ID, different lengths in A and B -> keyed by (genome, gene), not overwritten.
    a = _write(tmp_path, "A.faa", ">g1\n" + "A" * 10 + "\n")
    b = _write(tmp_path, "B.faa", ">g1\n" + "A" * 20 + "\n")
    text = "Orthogroup\tA\tB\nOG1\tg1\tg1\n"
    row = basic_tracks(_write(tmp_path, "og.tsv", text), [a, b])["OG1"]
    assert row["length_cv"] == pytest.approx(statistics.pstdev([10, 20]) / 15)


def test_ambiguous_id_unresolvable_is_dropped(tmp_path, capsys):
    # Ambiguous ID but the file stems don't match the matrix labels -> unresolved.
    a = _write(tmp_path, "A.faa", ">g1\n" + "A" * 10 + "\n")
    b = _write(tmp_path, "B.faa", ">g1\n" + "A" * 20 + "\n")
    text = "Orthogroup\tX\tY\nOG1\tg1\tg1\n"
    row = basic_tracks(_write(tmp_path, "og.tsv", text), [a, b])["OG1"]
    assert row["length_cv"] is None
    assert "could not be mapped" in capsys.readouterr().err


def test_overwide_row_raises(tmp_path):
    text = "Orthogroup\tA\tB\nOG1\tg1\nOG2\tg1\tg2\tg3\n"
    with pytest.raises(ValueError, match="line 3"):
        basic_tracks(_write(tmp_path, "bad.tsv", text))


def test_zero_genome_columns_raises(tmp_path):
    with pytest.raises(ValueError, match="no genome columns"):
        basic_tracks(_write(tmp_path, "z.tsv", "OG1\nOG2\n"))


def test_write_tracks_missing_inf_and_ordinary(tmp_path):
    from orthopan.tracks import write_tracks

    rows = {"OG1": {"a": None, "b": float("nan"), "c": float("inf"),
                    "d": float("-inf"), "e": 0.0, "f": 1.5}}
    out = tmp_path / "t.tsv"
    write_tracks(rows, ["a", "b", "c", "d", "e", "f"], out)
    assert out.read_text().splitlines() == [
        "orthogroup\ta\tb\tc\td\te\tf",
        "OG1\t\t\t\t\t0\t1.5",
    ]
