import gzip

import pytest

from orthopan.io import read_matrix, split_genes


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


def test_close_before_iterating_closes_file(tmp_path):
    _, rows = read_matrix(_write(tmp_path, "c.tsv", "Orthogroup\tA\tB\nOG1\tg1\tg2\n"))
    assert not rows._fh.closed
    rows.close()
    assert rows._fh.closed


def test_iterating_to_exhaustion_closes_file(tmp_path):
    _, rows = read_matrix(_write(tmp_path, "e.tsv", "Orthogroup\tA\tB\nOG1\tg1\tg2\n"))
    assert list(rows) == [("OG1", ["g1", "g2"])]
    assert rows._fh.closed


def test_with_closes_file_when_body_raises(tmp_path):
    _, rows = read_matrix(_write(tmp_path, "w.tsv", "Orthogroup\tA\tB\nOG1\tg1\tg2\n"))
    with pytest.raises(RuntimeError):
        with rows:
            assert not rows._fh.closed
            raise RuntimeError("boom")
    assert rows._fh.closed


# --- Format A: flat two-column (genome from gene ID) ------------------------

FORMAT_A = (
    "N0.HOG0000000\tgnl|extdb|FAM1-i1-1.1_001153, gnl|extdb|FAM2-i1-1.1_001100\n"
    "N0.HOG0000001\tgnl|extdb|FAM1-i1-1.1_001457, gnl|extdb|FAM1-i1-1.1_001458\n"
    "N0.HOG0000002\t\n"
)


def _read_a(path, **kw):
    labels, rows = read_matrix(path, **kw)
    return labels, list(rows)


def test_format_a_genome_from_gene(tmp_path):
    labels, rows = _read_a(_write(tmp_path, "a.tsv", FORMAT_A))
    assert labels == ["FAM1-i1-1.1", "FAM2-i1-1.1"]
    assert rows[0] == (
        "N0.HOG0000000",
        ["gnl|extdb|FAM1-i1-1.1_001153", "gnl|extdb|FAM2-i1-1.1_001100"],
    )
    # both genes of row 2 come from the same genome, joined back into one cell
    assert rows[1] == (
        "N0.HOG0000001",
        ["gnl|extdb|FAM1-i1-1.1_001457, gnl|extdb|FAM1-i1-1.1_001458", ""],
    )
    assert rows[2] == ("N0.HOG0000002", ["", ""])


def test_format_a_gz(tmp_path):
    raw = gzip.compress(FORMAT_A.encode())
    labels, rows = _read_a(_write(tmp_path, "a.bin", raw))
    assert labels == ["FAM1-i1-1.1", "FAM2-i1-1.1"]
    assert split_genes(rows[1][1][0]) == [
        "gnl|extdb|FAM1-i1-1.1_001457", "gnl|extdb|FAM1-i1-1.1_001458"]


def test_format_a_prefix_disabled(tmp_path):
    labels, rows = _read_a(_write(tmp_path, "a.tsv", FORMAT_A), strip_prefix="")
    assert labels == ["gnl|extdb|FAM1-i1-1.1", "gnl|extdb|FAM2-i1-1.1"]


def test_format_a_no_underscore_gene(tmp_path):
    text = "HOG1\tgnl|extdb|genomeX, gnl|extdb|genomeY_1\n"
    labels, rows = _read_a(_write(tmp_path, "a.tsv", text))
    assert labels == ["genomeX", "genomeY"]
    assert rows[0][1] == ["gnl|extdb|genomeX", "gnl|extdb|genomeY_1"]


def test_format_a_two_column_with_header_is_matrix(tmp_path):
    # A header makes it Format B even with one genome column.
    _, rows = _read_a(_write(tmp_path, "h.tsv", "Orthogroup\tA\nOG1\tgene\n"))
    assert rows == [("OG1", ["gene"])]


def test_format_a_close_semantics(tmp_path):
    _, rows = read_matrix(_write(tmp_path, "a.tsv", FORMAT_A))
    assert not rows._fh.closed
    rows.close()
    assert rows._fh.closed
