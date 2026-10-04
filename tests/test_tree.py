import gzip
import itertools
import re

import pytest

from orthopan.tree import _matrix, upgma_newick

TSV = (
    "Orthogroup\tA\tB\tC\tD\n"
    "OG1\tx\tx\t\t\n"
    "OG2\tx\tx\t\t\n"
    "OG3\tx\tx\t\t\n"
    "OG4\tx\t\t\t\n"
    "OG5\t\tx\t\t\n"
    "OG6\t\t\tx\tx\n"
    "OG7\t\t\tx\tx\n"
    "OG8\t\t\t\tx\n"
)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def _tip_paths(newick):
    """Map frozenset({a, b}) -> sum of branch lengths between two tips."""
    toks = re.findall(r"[(),;]|:[^(),;]+|[^(),;:]+", newick)
    pos = 0
    pairs = {}

    def clade():
        nonlocal pos
        if toks[pos] == "(":
            pos += 1
            left, bl = clade()
            pos += 1  # ','
            right, br = clade()
            pos += 1  # ')'
            left = {k: v + bl for k, v in left.items()}
            right = {k: v + br for k, v in right.items()}
            for a, da in left.items():
                for b, db in right.items():
                    pairs[frozenset((a, b))] = da + db
            merged = {**left, **right}
        else:
            merged = {toks[pos]: 0.0}
            pos += 1
        branch = 0.0
        if pos < len(toks) and toks[pos].startswith(":"):
            branch = float(toks[pos][1:])
            pos += 1
        return merged, branch

    clade()
    return pairs


def test_upgma_topology(tmp_path):
    nwk = upgma_newick(_write(tmp_path, "a.tsv", TSV))

    assert nwk.endswith(";")
    assert nwk.count("(") == nwk.count(")") == 3
    topology = re.sub(r":[\d.eE+-]+", "", nwk)
    assert "(A,B)" in topology or "(B,A)" in topology
    assert "(C,D)" in topology or "(D,C)" in topology


def test_gz_input(tmp_path):
    plain = _write(tmp_path, "a.tsv", TSV)
    gz = tmp_path / "a.tsv.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(TSV)
    assert upgma_newick(gz) == upgma_newick(plain)


def test_branch_lengths_disjoint(tmp_path):
    # Disjoint single-gene genomes -> Jaccard 1 -> each tip branch 1/2.
    nwk = upgma_newick(_write(tmp_path, "d.tsv", "Orthogroup\tA\tB\nOG1\tx\t\nOG2\t\tx\n"))
    paths = _tip_paths(nwk)
    assert paths == {frozenset(("A", "B")): pytest.approx(1.0)}
    assert re.findall(r":([\d.]+)", nwk) == ["0.5", "0.5"]


def test_tip_to_tip_equals_jaccard(tmp_path):
    # Ultrametric case: (A,B) and (C,D) closer than the cross distances.
    tsv = (
        "Orthogroup\tA\tB\tC\tD\n"
        "OG1\tx\tx\t\t\n"
        "OG2\tx\tx\t\t\n"
        "OG3\tx\t\t\t\n"
        "OG4\t\tx\t\t\n"
        "OG5\t\t\tx\tx\n"
        "OG6\t\t\tx\tx\n"
        "OG7\t\t\tx\t\n"
        "OG8\t\t\t\tx\n"
    )
    sets = {"A": {1, 2, 3}, "B": {1, 2, 4}, "C": {5, 6, 7}, "D": {5, 6, 8}}
    expected = {
        frozenset((a, b)): 1 - len(sets[a] & sets[b]) / len(sets[a] | sets[b])
        for a, b in itertools.combinations(sorted(sets), 2)
    }
    paths = _tip_paths(upgma_newick(_write(tmp_path, "u.tsv", tsv)))
    assert paths.keys() == expected.keys()
    for key, value in expected.items():
        assert paths[key] == pytest.approx(value)


def test_deep_tree_is_iterative(tmp_path):
    names = [f"g{i}" for i in range(1100)]
    header = "Orthogroup\t" + "\t".join(names)
    present = "OG1\t" + "\t".join(["x"] * 1100)
    absent = "OG2\t" + "\t".join([""] * 1100)
    path = _write(tmp_path, "deep.tsv", f"{header}\n{present}\n{absent}\n")
    nwk = upgma_newick(path)
    assert nwk.endswith(";")
    assert nwk.count("(") == nwk.count(")") == 1099


def test_quote_special_chars(tmp_path):
    tsv = (
        "Orthogroup\tA[1]\tB c\tD's\n"
        "OG1\tx\tx\t\n"
        "OG2\tx\t\tx\n"
        "OG3\t\tx\tx\n"
    )
    nwk = upgma_newick(_write(tmp_path, "q.tsv", tsv))
    assert "'A[1]'" in nwk
    assert "'B c'" in nwk
    assert "'D''s'" in nwk


def test_single_genome(tmp_path):
    assert upgma_newick(_write(tmp_path, "one.tsv", "Orthogroup\tsolo\nOG1\tx\nOG2\tx\n")) == "solo;"
    one = _write(tmp_path, "one2.tsv", "Orthogroup\ta b\nOG1\tx\nOG2\tx\n")
    assert upgma_newick(one) == "'a b';"


def test_pads_short_rows(tmp_path):
    # Width comes from the header; short first and later rows are padded with absent.
    tsv = "Orthogroup\tA\tB\tC\nOG1\tx\nOG2\t\tx\nOG3\tx\n"
    labels, m = _matrix(_write(tmp_path, "p.tsv", tsv))
    assert labels == ["A", "B", "C"]
    assert m.tolist() == [[True, False, True], [False, True, False], [False, False, False]]


def test_reproducer_short_rows(tmp_path):
    tsv = "Orthogroup\tA\tB\tC\nOG1\tx\t\t\nOG2\tx\nOG3\tx\nOG4\tx\n"
    labels, m = _matrix(_write(tmp_path, "r.tsv", tsv))
    assert labels == ["A", "B", "C"]
    assert m.tolist() == [[True, True, True, True], [False] * 4, [False] * 4]


def test_rejects_too_long_row(tmp_path):
    tsv = "Orthogroup\tA\tB\nOG1\tx\t\nOG2\t\tx\nOG3\tx\tx\tx\n"
    with pytest.raises(ValueError, match="line 4"):
        _matrix(_write(tmp_path, "l.tsv", tsv))


def test_header_only_allowed(tmp_path):
    path = _write(tmp_path, "h.tsv", "Orthogroup\tA\tB\n")
    labels, m = _matrix(path)
    assert labels == ["A", "B"]
    assert m.shape == (2, 0)
    assert upgma_newick(path).endswith(";")


def test_empty_input_raises(tmp_path):
    with pytest.raises(ValueError, match="no genome columns"):
        _matrix(_write(tmp_path, "e.tsv", ""))


def test_no_genome_columns_raises(tmp_path):
    with pytest.raises(ValueError, match="no genome columns"):
        _matrix(_write(tmp_path, "ng.tsv", "OG1\nOG2\n"))
