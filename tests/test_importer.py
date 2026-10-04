import json
import random

import numpy as np
import polars as pl
import pytest

from orthopan import main
from orthopan.importer import _parse_newick, _Tree

TREE = "((A:1,B:1)N3:0.5,(C:1,D:1)N4:0.5)N0;"

# header order deliberately differs from tree DFS order (C,A,D,B)
MATRIX = (
    "Orthogroup\tC\tA\tD\tB\n"
    "OG_core\tx\tx\tx\tx\n"
    "OG_ab\t\tx\t\tx\n"
    "OG_cd\tx\t\tx\t\n"
    "OG_ac\tx\tx\t\t\n"
    "OG_a\t\tga1, ga2, ga3\t\t\n"
    "OG_c\tx\t\t\t\n"
    "OG_empty\t\t\t\t\n"
)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def _import(tmp_path, matrix=MATRIX, tree=TREE, tracks=()):
    argv = [str(_write(tmp_path, "matrix.tsv", matrix)), "-o", str(tmp_path / "out")]
    if tree is not None:
        argv += ["--tree", str(_write(tmp_path, "tree.nwk", tree))]
    for i, text in enumerate(tracks):
        argv += ["--track", str(_write(tmp_path, f"track{i}.tsv", text))]
    main(["import", *argv])
    return tmp_path / "out"


def test_genome_rows_follow_tree_order(tmp_path):
    out = _import(tmp_path)
    genomes = pl.read_parquet(out / "genomes.parquet").sort("genome_idx")
    # input column order C,A,D,B -> leaf positions 2,0,3,1
    assert genomes["name"].to_list() == ["C", "A", "D", "B"]
    assert genomes["row_order"].to_list() == [2, 0, 3, 1]
    assert genomes["n_orthogroups"].to_list() == [4, 4, 2, 2]  # C,A,D,B
    assert genomes["n_genes"].to_list() == [4, 6, 2, 2]


def test_mrca_purity_and_order_tree(tmp_path):
    out = _import(tmp_path)
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    assert og["orthogroup"].to_list() == [
        "OG_core", "OG_ab", "OG_cd", "OG_ac", "OG_a", "OG_c", "OG_empty"]
    assert og["n_genomes"].to_list() == [4, 2, 2, 2, 1, 1, 0]
    assert og["mrca_node"].to_list() == [0, 1, 4, 0, 2, 5, None]
    assert og["mrca_label"].to_list() == ["N0", "N3", "N4", "N0", "A", "C", None]
    assert og["purity"].to_list() == [1.0, 1.0, 1.0, 0.5, 1.0, 1.0, None]
    # order_tree: core(root), ac(root,purity .5), ab(N3), a(A), cd(N4), c(C), empty
    assert og["order_tree"].to_list() == [0, 2, 4, 1, 3, 5, 6]
    assert og["copy_number_mean"].to_list() == [1.0, 1.0, 1.0, 1.0, 3.0, 1.0, None]


def test_order_frequency(tmp_path):
    out = _import(tmp_path)
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    # frequency: core(4), ab(2), cd(2), a(1), c(1), empty
    assert og["order_frequency"].to_list() == [0, 1, 2, 3, 4, 5, 6]


def test_presence_bits_roundtrip(tmp_path):
    out = _import(tmp_path)
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    n_og = og.height
    arr = np.load(out / "presence.npy")
    assert arr.shape == (4, (n_og + 7) // 8)
    unpacked = np.unpackbits(arr, axis=1)[:, :n_og].astype(bool)
    membership = pl.read_parquet(out / "membership.parquet")
    dense = np.zeros((4, n_og), bool)
    for row in membership.iter_rows(named=True):
        dense[row["genome_idx"], row["og_idx"]] = True
    assert np.array_equal(unpacked, dense)


def test_membership_sparse_cells(tmp_path):
    out = _import(tmp_path)
    membership = pl.read_parquet(out / "membership.parquet")
    assert membership.height == 12  # 4+2+2+2+1+1 present cells
    a_row = membership.filter((pl.col("og_idx") == 4)).row(0, named=True)
    assert a_row["genome_idx"] == 1 and a_row["copies"] == 3
    assert a_row["genes"] == ["ga1", "ga2", "ga3"]


def test_track_join_and_unknown_warning(tmp_path, capsys):
    track = "orthogroup\tscore\nOG_core\t1.5\nOG_ab\t\nOG_nothere\t9.0\n"
    out = _import(tmp_path, tracks=[track])
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    assert og["score"].to_list() == [1.5, None, None, None, None, None, None]
    assert json.loads((out / "manifest.json").read_text())["track_columns"] == ["score"]
    assert "1 track orthogroup(s) not in the matrix" in capsys.readouterr().err


def test_duplicate_track_column_errors(tmp_path):
    t1 = "orthogroup\tscore\nOG_core\t1\n"
    t2 = "orthogroup\tscore\nOG_core\t2\n"
    with pytest.raises(ValueError, match="duplicate track column name"):
        _import(tmp_path, tracks=[t1, t2])


def test_tree_label_mismatch_errors(tmp_path):
    with pytest.raises(ValueError, match="tree leaves do not match genome labels"):
        _import(tmp_path, tree="((A:1,B:1)N3,(C:1,E:1)N4)N0;")


def test_upgma_fallback_when_tree_omitted(tmp_path):
    out = _import(tmp_path, tree=None)
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["tree_source"] == "upgma"
    assert (out / "tree.nwk").read_text().strip().endswith(";")
    genomes = pl.read_parquet(out / "genomes.parquet")
    assert set(genomes["row_order"].to_list()) == {0, 1, 2, 3}
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    assert og["mrca_node"].to_list()[-1] is None


def test_deep_caterpillar_tree_is_iterative(tmp_path):
    n = 2000
    newick = "g0"
    for i in range(1, n):
        newick = f"({newick},g{i})"
    tree = _Tree(_parse_newick(newick + ";"))
    assert tree.n_leaves == n
    assert tree.hi[0] == n - 1
    # MRCA of the two extreme leaves is the root.
    assert tree.mrca(0, n - 1) == 0


def test_cli_through_main(tmp_path):
    out = _import(tmp_path)
    for name in ("genomes.parquet", "orthogroups.parquet", "membership.parquet",
                 "presence.npy", "tree.nwk", "manifest.json"):
        assert (out / name).is_file()


def test_newick_quoted_labels_and_branches():
    tree = _Tree(_parse_newick("('a b':1e-2,'c,d''s':0.5)N0;"))
    assert tree.n_leaves == 2
    assert sorted(tree.label[i] for i in tree.leaf_order) == ["a b", "c,d's"]


# --- track columns (collisions, non-finite values) --------------------------

@pytest.mark.parametrize(
    "name",
    [
        "og_idx", "orthogroup", "n_genomes", "n_genes", "copy_number_mean",
        "copy_number_max", "mrca_node", "mrca_label", "purity", "order_tree",
        "order_frequency",
    ],
)
def test_track_column_collision_with_builtin_errors(tmp_path, name):
    track = f"orthogroup\t{name}\nOG_core\t1\n"
    with pytest.raises(ValueError, match="duplicate track column name"):
        _import(tmp_path, tracks=[track])


def test_nonfinite_track_values_become_null(tmp_path):
    track = "orthogroup\tscore\nOG_core\tnan\nOG_ab\tinf\nOG_cd\t-inf\nOG_ac\t2.5\n"
    out = _import(tmp_path, tracks=[track])
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    assert og["score"].to_list() == [None, None, None, 2.5, None, None, None]


# --- duplicate labels -------------------------------------------------------

def test_duplicate_matrix_label_errors(tmp_path):
    matrix = "Orthogroup\tA\tA\tB\nOG1\tx\tx\t\n"
    with pytest.raises(ValueError, match="duplicate genome labels in matrix"):
        _import(tmp_path, matrix=matrix)


def test_duplicate_tree_leaf_errors(tmp_path):
    with pytest.raises(ValueError, match="duplicate leaf labels in tree"):
        _import(tmp_path, tree="(A,A,B,C,D)N0;")


def test_each_genome_has_exactly_one_leaf(tmp_path):
    # A tree that names every genome once imports fine; a duplicate leaf fails.
    out = _import(tmp_path, tree="((A:1,B:1)N3,(C:1,D:1)N4)N0;")
    og = pl.read_parquet(out / "orthogroups.parquet")
    assert og.height == 7


# --- Newick comments and quoted punctuation ---------------------------------

def test_newick_comments_are_skipped():
    tree = _Tree(_parse_newick("(A[comment]:1,B:2[after])R;"))
    assert sorted(tree.label[i] for i in tree.leaf_order) == ["A", "B"]


def test_newick_nested_comments_are_skipped():
    tree = _Tree(_parse_newick("(A[x[y]z]:1,B:2)R;"))
    assert sorted(tree.label[i] for i in tree.leaf_order) == ["A", "B"]


def test_newick_quoted_punctuation_is_a_single_label():
    tree = _Tree(_parse_newick("('(A,B)':1,':;'[note]:2)R;"))
    assert sorted(tree.label[i] for i in tree.leaf_order) == ["(A,B)", ":;"]


def test_upgma_fallback_quotes_punctuation_labels(tmp_path):
    matrix = "Orthogroup\tA,B\tC(D)\nOG1\tx\t\n"
    out = _import(tmp_path, matrix=matrix, tree=None)
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["tree_source"] == "upgma"
    tree = _Tree(_parse_newick((out / "tree.nwk").read_text()))
    assert sorted(tree.label[i] for i in tree.leaf_order) == ["A,B", "C(D)"]


# --- malformed Newick -------------------------------------------------------

@pytest.mark.parametrize(
    "text",
    [
        "(A,B;",            # unbalanced '('
        "(A,B))R;",         # unbalanced ')'
        "'A,B;",            # unterminated quote
        "(A[comment,B);",   # unterminated comment
        "(A:,B)R;",         # no branch length after ':'
        "(A:B)R;",          # non-numeric branch length
        "(A,B)R; extra",    # trailing content after ';'
        "(A,B)R;;",         # a second ';' is trailing content
        "(A B)R;",          # missing comma between leaves
        "(A:1:2,B)R;",      # repeated branch length
        "(A,B)R:",          # missing length at EOF
        "A,(B,C)R;",        # a second root
    ],
)
def test_malformed_newick_raises(text):
    with pytest.raises(ValueError):
        _parse_newick(text)


def test_valid_newick_without_trailing_semicolon():
    tree = _Tree(_parse_newick("(A:1,B:2)R"))
    assert tree.n_leaves == 2


# --- trees: polytomies, ties, MRCA ------------------------------------------

def test_polytomy_mrca(tmp_path):
    tree = "((A:1,B:1,C:1)N3:0.5,D:1)N0;"
    matrix = "Orthogroup\tA\tB\tC\tD\nOG_abc\tx\tx\tx\t\n"
    out = _import(tmp_path, matrix=matrix, tree=tree)
    og = pl.read_parquet(out / "orthogroups.parquet")
    assert og["mrca_label"].to_list() == ["N3"]
    assert og["n_genomes"].to_list() == [3]
    assert og["purity"].to_list() == [1.0]


def test_same_clade_rank_ties_break_by_index(tmp_path):
    tree = "((A:1,B:1)N3:0.5,(C:1,D:1)N4:0.5)N0;"
    matrix = (
        "Orthogroup\tA\tB\tC\tD\n"
        "OG_a\tx\t\t\t\n"
        "OG_ab\tx\tx\t\t\n"
        "OG_abd\tx\tx\t\tx\n"
        "OG_abc\tx\tx\tx\t\n"
    )
    out = _import(tmp_path, matrix=matrix, tree=tree)
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    # OG_abd and OG_abc share MRCA N0 and purity 3/4; index order breaks the tie.
    assert og["mrca_label"].to_list() == ["A", "N3", "N0", "N0"]
    assert og["purity"].to_list() == [1.0, 1.0, 0.75, 0.75]
    assert og["order_tree"].to_list() == [3, 2, 0, 1]


def _random_newick(rng, n):
    nodes = [chr(ord("a") + i) for i in range(n)]
    while len(nodes) > 1:
        i, j = sorted(rng.sample(range(len(nodes)), 2))
        nodes[i] = f"({nodes[i]},{nodes[j]})"
        del nodes[j]
    return nodes[0] + ";"


def test_mrca_matches_bruteforce_on_random_trees():
    rng = random.Random(1234)
    for _ in range(40):
        newick = _random_newick(rng, rng.randint(2, 8))
        tree = _Tree(_parse_newick(newick))
        for a in range(tree.n_leaves):
            for b in range(a, tree.n_leaves):
                expected = min(
                    (i for i in range(tree.n)
                     if tree.lo[i] <= a and tree.hi[i] >= b),
                    key=lambda i: tree.hi[i] - tree.lo[i],
                )
                assert tree.mrca(a, b) == expected, (newick, a, b)


# --- membership batching ----------------------------------------------------

def test_membership_multiple_batches_consolidate(tmp_path, monkeypatch):
    rows = "".join(f"OG{i:02d}\tx\t\t\tx\n" for i in range(17))
    matrix = "Orthogroup\tC\tA\tD\tB\n" + rows
    monkeypatch.setattr("orthopan.importer._MEMBERSHIP_BATCH", 2)
    out = _import(tmp_path, matrix=matrix)
    membership = pl.read_parquet(out / "membership.parquet")
    assert membership.height == 34  # 17 OGs x 2 genomes
    assert membership["og_idx"].n_unique() == 17
    assert not list(out.glob("_membership.*.parquet"))  # parts removed
    arr = np.unpackbits(np.load(out / "presence.npy"), axis=1)[:, :17].astype(bool)
    dense = np.zeros((4, 17), bool)
    for row in membership.iter_rows(named=True):
        dense[row["genome_idx"], row["og_idx"]] = True
    assert np.array_equal(arr, dense)
