import gzip

import pytest

from orthopan.io import read_matrix


def _write(tmp_path, name, data):
    path = tmp_path / name
    if isinstance(data, bytes):
        path.write_bytes(data)
    else:
        path.write_text(data)
    return path


def _read(path):
    labels, rows = read_matrix(path)
    return labels, list(rows)


def test_header(tmp_path):
    text = "Orthogroup\tA\tB\nOG1\tg1\tg2\nOG2\tg1\t\n"
    labels, rows = _read(_write(tmp_path, "h.tsv", text))
    assert labels == ["A", "B"]
    assert rows == [("OG1", ["g1", "g2"]), ("OG2", ["g1", ""])]


def test_headerless_keeps_first_row(tmp_path):
    text = "OG1\tg1\tg2\nOG2\t\tg2\n"
    labels, rows = _read(_write(tmp_path, "nh.tsv", text))
    assert labels == ["genome_1", "genome_2"]
    assert rows == [("OG1", ["g1", "g2"]), ("OG2", ["", "g2"])]


def test_gz_detected_by_magic_without_gz_suffix(tmp_path):
    text = "Orthogroup\tA\tB\nOG1\tg1\tg2\n"
    raw = gzip.compress(text.encode())
    labels, rows = _read(_write(tmp_path, "matrix.bin", raw))
    assert labels == ["A", "B"]
    assert rows == [("OG1", ["g1", "g2"])]


def test_short_rows_padded(tmp_path):
    text = "Orthogroup\tA\tB\tC\nOG1\tg1\nOG2\t\tg2\n"
    labels, rows = _read(_write(tmp_path, "s.tsv", text))
    assert labels == ["A", "B", "C"]
    assert rows == [("OG1", ["g1", "", ""]), ("OG2", ["", "g2", ""])]


def test_overwide_row_raises_line_number(tmp_path):
    text = "Orthogroup\tA\tB\nOG1\tg1\nOG2\tg1\tg2\tg3\n"
    with pytest.raises(ValueError, match="line 3"):
        _read(_write(tmp_path, "wide.tsv", text))


def test_empty_input_raises(tmp_path):
    with pytest.raises(ValueError, match="no genome columns"):
        read_matrix(_write(tmp_path, "e.tsv", ""))


def test_single_id_no_columns_raises(tmp_path):
    with pytest.raises(ValueError, match="no genome columns"):
        read_matrix(_write(tmp_path, "z.tsv", "OG1\nOG2\n"))
