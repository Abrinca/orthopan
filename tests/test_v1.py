"""Tests for the v1 additions: annotations, metadata, patterns, copies, search."""

import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from fastapi.testclient import TestClient

from orthopan import main
from orthopan.server import create_app

TREE = "((A:1,B:1)N3:0.5,(C:1,D:1)N4:0.5)N0;"

# multi-copy cells, and two OGs (OG_f/OG_g) sharing a presence pattern
MATRIX = (
    "Orthogroup\tC\tA\tD\tB\n"
    "OG_a\t\tx,y\t\t\n"
    "OG_b\t\tx\tx,y,z\t\n"
    "OG_c\tx,y\tx\tx\t\n"
    "OG_d\t\t\t\tx\n"
    "OG_e\tx\tx\tx\tx\n"
    "OG_f\tx\tx\t\t\n"
    "OG_g\tx\tx\t\t\n"
)
DESCRIPTIONS = (
    "OG_a\tamino acid transporter ATP-binding protein\n"
    "OG_b\tabc transporter permease\n"
    "OG_c\tkinase\n"
    "OG_d\thypothetical protein\n"
    "OG_e\ttransporter\n"
    "OG_f\tkinase\n"
    "OG_g\tkinase\n"
)
COG = (
    "Orthogroup\tcog_category\n"
    "OG_a\tJ\nOG_b\tK\nOG_c\tJ\nOG_d\tK\nOG_e\tL\nOG_f\tJ\nOG_g\tJ\n"
)
GENOME_METADATA = "genome\tclade\tsize\nA\tX\t1.5\nB\tX\t2.5\nC\tY\t3.5\nD\tY\t4.5\n"
DEMO = (Path(__file__).resolve().parent.parent
        / "resources/orthofinder-demo/results/Results_Oct04_1")


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def _import(tmp_path, matrix=MATRIX, tree=TREE, extras=(), name="out"):
    argv = [str(_write(tmp_path, "matrix.tsv", matrix)), "-o", str(tmp_path / name)]
    if tree is not None:
        argv += ["--tree", str(_write(tmp_path, "tree.nwk", tree))]
    argv += list(extras)
    main(["import", *argv])
    return tmp_path / name


@pytest.fixture(scope="module")
def out(tmp_path_factory):
    d = tmp_path_factory.mktemp("v1")
    return _import(d, extras=[
        "--annotation", str(_write(d, "desc-descriptions.tsv", DESCRIPTIONS)),
        "--annotation", str(_write(d, "cog.tsv", COG)),
        "--genome-metadata", str(_write(d, "gm.tsv", GENOME_METADATA)),
    ])


@pytest.fixture(scope="module")
def client(out):
    return TestClient(create_app(out))


# --- prefix stripping / raw retention ---------------------------------------

def test_strip_prefix_and_raw_retention(tmp_path):
    matrix = "Orthogroup\tC\tA\tD\tB\nOG1\t\tgnl|extdb|A_1, gnl|extdb|A_2\t\t\n"
    out = _import(tmp_path, matrix=matrix)
    row = pl.read_parquet(out / "membership.parquet").row(0, named=True)
    assert row["genes"] == ["A_1", "A_2"]
    assert row["genes_raw"] == ["gnl|extdb|A_1", "gnl|extdb|A_2"]


def test_strip_prefix_disabled(tmp_path):
    matrix = "Orthogroup\tC\tA\tD\tB\nOG1\t\tgnl|extdb|A_1\t\t\n"
    out = _import(tmp_path, matrix=matrix, extras=["--strip-prefix", ""])
    row = pl.read_parquet(out / "membership.parquet").row(0, named=True)
    assert row["genes"] == ["gnl|extdb|A_1"]
    assert row["genes_raw"] == ["gnl|extdb|A_1"]


def test_format_a_import_without_tree(tmp_path):
    # headerless two-column -> Format A, genome from the gene ID; UPGMA fallback
    matrix = (
        "HOG1\tgnl|extdb|A_1, gnl|extdb|B_2\n"
        "HOG2\tgnl|extdb|A_3\n"
    )
    out = _import(tmp_path, matrix=matrix, tree=None, name="fa")
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["tree_source"] == "upgma"
    genomes = pl.read_parquet(out / "genomes.parquet").sort("genome_idx")
    assert set(genomes["name"].to_list()) == {"A", "B"}
    assert pl.read_parquet(out / "orthogroups.parquet").height == 2


def test_orthogroup_endpoint_gene_objects(client):
    body = client.get("/api/orthogroup/OG_e").json()
    genes = {m["genome"]: m["genes"] for m in body["members"]}
    assert set(genes) == {"A", "B", "C", "D"}
    assert genes["A"] == [{"id": "x", "raw": "x"}]


# --- annotations -------------------------------------------------------------

def test_annotation_specs_and_text_categorical(out):
    manifest = json.loads((out / "manifest.json").read_text())
    specs = {s["name"]: s for s in manifest["annotations"]}
    assert specs["description"]["kind"] == "text"
    assert specs["cog_category"]["kind"] == "categorical"
    assert set(specs["cog_category"]["categories"]) == {"J", "K", "L"}
    # most common first
    assert specs["cog_category"]["categories"][0] == "J"
    ann = pl.read_parquet(out / "annotations.parquet").sort("og_idx")
    assert ann["description"].to_list()[0] == "amino acid transporter ATP-binding protein"


def test_annotation_typing_by_distinct_count(tmp_path):
    n = 70
    matrix = "Orthogroup\tA\tB\n" + "".join(
        f"OG{i:03d}\tx\t\n" for i in range(n))
    rows = ["Orthogroup\tcat64\ttext65"]
    for i in range(n):
        rows.append(f"OG{i:03d}\tc{i % 64}\tv{i}")
    out = _import(tmp_path, matrix=matrix, tree="((A:1,B:1)N0);",
                  extras=["--annotation", str(_write(tmp_path, "ann.tsv", "\n".join(rows) + "\n"))])
    specs = {s["name"]: s for s in json.loads((out / "manifest.json").read_text())["annotations"]}
    assert specs["cat64"]["kind"] == "categorical"
    assert len(specs["cat64"]["categories"]) == 64
    assert specs["text65"]["kind"] == "text"


def test_annotation_headerless_default_name(tmp_path):
    text = "OG_a\thello\nOG_b\tworld\n"
    out = _import(tmp_path, extras=[
        "--annotation", str(_write(tmp_path, "my-descriptions.tsv", text))])
    specs = json.loads((out / "manifest.json").read_text())["annotations"]
    assert specs == [{"name": "description", "kind": "text"}]


def test_annotation_unknown_og_warns(tmp_path, capsys):
    text = "OG_a\thello\nOG_nothere\tworld\n"
    out = _import(tmp_path, extras=["--annotation", str(_write(tmp_path, "ann.tsv", text))])
    assert "not in the matrix" in capsys.readouterr().err


def test_annotation_unknown_ids_do_not_affect_specs(tmp_path, capsys):
    n = 70
    matrix = "Orthogroup\tA\tB\n" + "".join(f"OG{i:03d}\tx\ty\n" for i in range(n))
    rows = ["Orthogroup\tcat"] + [f"OG{i:03d}\tA" for i in range(n)]
    rows += [f"UNKNOWN{i}\tu{i}" for i in range(100)]
    out = _import(
        tmp_path, matrix=matrix, tree="((A:1,B:1)N0);",
        extras=["--annotation", str(_write(tmp_path, "ann.tsv", "\n".join(rows) + "\n"))])
    assert "not in the matrix" in capsys.readouterr().err
    spec = json.loads((out / "manifest.json").read_text())["annotations"][0]
    # unknown IDs would have pushed the distinct count to 101 -> text
    assert spec == {"name": "cat", "kind": "categorical", "categories": ["A"]}


def test_annotation_duplicate_name_clash(tmp_path):
    text = "OG_a\thello\n"
    ann = str(_write(tmp_path, "ann.tsv", text))
    with pytest.raises(ValueError, match="duplicate annotation column name"):
        _import(tmp_path, extras=["--annotation", ann, "--annotation", ann])


def test_annotation_builtin_name_clash(tmp_path):
    text = "Orthogroup\tpurity\nOG_a\t1\n"
    with pytest.raises(ValueError, match="duplicate annotation column name"):
        _import(tmp_path, extras=["--annotation", str(_write(tmp_path, "ann.tsv", text))])


# --- genome metadata ---------------------------------------------------------

def test_genome_metadata_typing(out):
    specs = {s["name"]: s for s in json.loads(
        (out / "manifest.json").read_text())["genome_metadata"]}
    assert specs["clade"]["kind"] == "categorical"
    assert set(specs["clade"]["categories"]) == {"X", "Y"}
    assert specs["size"]["kind"] == "numeric"
    assert (specs["size"]["min"], specs["size"]["max"]) == (1.5, 4.5)
    gm = pl.read_parquet(out / "genome_metadata.parquet").sort("genome_idx")
    genomes = pl.read_parquet(out / "genomes.parquet").sort("genome_idx")
    assert dict(zip(genomes["name"].to_list(), gm["size"].to_list())) == {
        "A": 1.5, "B": 2.5, "C": 3.5, "D": 4.5}


def test_genome_metadata_unknown_and_missing(tmp_path, capsys):
    # Z is unknown and carries a category ("Q") that must not enter the specs.
    gm = "genome\tclade\tgc\nA\tX\t0.4\nC\tY\t0.5\nZ\tQ\t0.9\n"
    out = _import(tmp_path, extras=["--genome-metadata", str(_write(tmp_path, "gm.tsv", gm))])
    err = capsys.readouterr().err
    assert "genome(s) not in the matrix" in err
    assert "1 genome(s) not in the matrix" in err
    spec = json.loads((out / "manifest.json").read_text())["genome_metadata"][0]
    assert set(spec["categories"]) == {"X", "Y"}  # Q from unknown Z excluded
    frame = pl.read_parquet(out / "genome_metadata.parquet").sort("genome_idx")
    genomes = pl.read_parquet(out / "genomes.parquet").sort("genome_idx")
    name_by_idx = dict(zip(genomes["genome_idx"].to_list(), genomes["name"].to_list()))
    by = {name_by_idx[int(g)]: c
          for g, c in zip(frame["genome_idx"].to_list(), frame["clade"].to_list())}
    assert by["B"] is None  # B missing from the file
    assert by["A"] == "X"


def test_genomes_endpoint_metadata(client):
    rows = {r["name"]: r for r in client.get("/api/genomes").json()}
    assert rows["A"]["metadata"] == {"clade": "X", "size": 1.5}
    assert rows["D"]["metadata"] == {"clade": "Y", "size": 4.5}


def test_meta_lists_annotations_and_metadata(client):
    body = client.get("/api/meta").json()
    assert {a["name"] for a in body["annotations"]} == {"description", "cog_category"}
    assert {m["name"] for m in body["genome_metadata"]} == {"clade", "size"}
    assert body["detail"] == {
        "msa": False, "gene_trees": False, "sequences": False, "famsa": False}


# --- pattern grouping --------------------------------------------------------

def test_pattern_grouping_bruteforce(out):
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    ids = og["orthogroup"].to_list()
    pids = og["pattern_id"].to_list()
    same = og["n_same_pattern"].to_list()
    dense = _presence(out)
    for i in range(len(ids)):
        for j in range(len(ids)):
            equal = bool(np.array_equal(dense[:, i], dense[:, j]))
            assert (pids[i] == pids[j]) == equal, (ids[i], ids[j])
    # OG_f and OG_g share a pattern
    assert pids[ids.index("OG_f")] == pids[ids.index("OG_g")]
    assert same[ids.index("OG_f")] == 2
    # the unique empty pattern
    assert same[ids.index("OG_d")] == 1


def test_pattern_endpoint(client):
    body = client.get("/api/orthogroup/OG_f/pattern").json()
    assert body["total"] == 2
    assert {o["orthogroup"] for o in body["orthogroups"]} == {"OG_f", "OG_g"}
    assert "annotations" in body["orthogroups"][0]
    assert client.get("/api/orthogroup/NOPE/pattern").status_code == 404


def test_pattern_order(client):
    r = client.post("/api/orders", json={"by": "pattern", "og": "OG_f"})
    assert r.status_code == 200
    oid = r.json()["order"]
    first = client.get("/api/orthogroups", params=dict(order=oid, og0=0, og1=2)).json()
    assert set(first) == {"OG_f", "OG_g"}


def test_orthogroup_detail_pattern_fields(client):
    body = client.get("/api/orthogroup/OG_f").json()
    assert body["n_same_pattern"] == 2
    assert body["annotations"]["description"] == "kinase"
    assert "pattern_id" in body


# --- copies ------------------------------------------------------------------

def _presence(out):
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    mg = pl.read_parquet(out / "genomes.parquet").sort("genome_idx")
    n_og = og.height
    inv = np.empty(mg.height, np.int64)
    inv[mg["row_order"].to_numpy()] = np.arange(mg.height)
    packed = np.load(out / "presence.npy")
    return np.unpackbits(packed, axis=1, bitorder="big")[:, :n_og][inv]


def _copies(out):
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    mg = pl.read_parquet(out / "genomes.parquet").sort("genome_idx")
    mem = pl.read_parquet(out / "membership.parquet")
    dense = np.zeros((mg.height, og.height), np.int32)
    for r in mem.iter_rows(named=True):
        dense[int(r["genome_idx"]), int(r["og_idx"])] = int(r["copies"])
    inv = np.empty(mg.height, np.int64)
    inv[mg["row_order"].to_numpy()] = np.arange(mg.height)
    return dense[inv]


def _copies_brute(copies_order, og0, og1, g0, g1, w, h):
    col = (np.arange(w + 1) * (og1 - og0)) // w
    row = (np.arange(h + 1) * (g1 - g0)) // h
    out = np.zeros((h, w), np.uint8)
    for i in range(h):
        for j in range(w):
            block = copies_order[g0 + row[i]:g0 + row[i + 1], og0 + col[j]:og0 + col[j + 1]]
            present = int((block > 0).sum())
            if not present:
                continue
            extra = int((block[block > 0] - 1).sum())
            out[i, j] = min(255, int(np.rint((present + extra) / present * 10.0)))
    return out


def _perms(out):
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    means = og["copy_number_mean"].to_numpy()
    missing = np.isnan(means)
    present = np.flatnonzero(~missing)
    track = np.concatenate([present[np.argsort(means[present], kind="stable")],
                            np.flatnonzero(missing)]).astype(np.int64)
    return {
        "tree": np.argsort(og["order_tree"].to_numpy(), kind="stable"),
        "frequency": np.argsort(og["order_frequency"].to_numpy(), kind="stable"),
        "track:copy_number_mean": track,
    }


def test_copies_matrix_bruteforce(out, client):
    copies = _copies(out)
    perms = _perms(out)
    n_genomes = copies.shape[0]
    n_og = copies.shape[1]
    for order, perm in perms.items():
        ordered = copies[:, perm]
        for og0, og1, w in [(0, n_og, n_og), (1, n_og, 3), (0, n_og, 4), (2, 6, 2)]:
            for g0, g1, h in [(0, n_genomes, n_genomes), (1, n_genomes, 2), (0, n_genomes, 3)]:
                r = client.get("/api/matrix", params=dict(
                    order=order, og0=og0, og1=og1, g0=g0, g1=g1, w=w, h=h,
                    value="copies"))
                assert r.status_code == 200
                got = np.frombuffer(r.content, "<u1").reshape(h, w)
                assert np.array_equal(got, _copies_brute(ordered, og0, og1, g0, g1, w, h)), \
                    (order, og0, og1, g0, g1, w, h)


def test_copies_dynamic_order(out, client):
    oid = client.post("/api/orders", json={"by": "contrast", "genomes": ["A", "B"]}).json()["order"]
    perm = client.app.state.pg.order_perm(oid)
    ordered = _copies(out)[:, perm]
    n_og, n_genomes = ordered.shape[1], ordered.shape[0]
    for og0, og1, w in [(0, n_og, 3), (1, n_og - 1, 4), (0, n_og, n_og)]:
        for g0, g1, h in [(0, n_genomes, 4), (1, n_genomes, 2)]:
            r = client.get("/api/matrix", params=dict(
                order=oid, og0=og0, og1=og1, g0=g0, g1=g1, w=w, h=h, value="copies"))
            got = np.frombuffer(r.content, "<u1").reshape(h, w)
            assert np.array_equal(got, _copies_brute(ordered, og0, og1, g0, g1, w, h)), \
                (og0, og1, g0, g1, w, h)


def test_copies_absent_bin_is_zero(out, client):
    n_og = len(client.get("/api/orthogroups", params=dict(og0=0, og1=7)).json())
    ids = client.get("/api/orthogroups", params=dict(order="tree", og0=0, og1=n_og)).json()
    pos = ids.index("OG_d")
    r = client.get("/api/matrix", params=dict(
        order="tree", og0=pos, og1=pos + 1, g0=0, g1=4, w=1, h=4, value="copies"))
    got = np.frombuffer(r.content, "<u1")
    # OG_d is present only in B (tree row 1, one copy)
    assert got.tolist() == [0, 10, 0, 0]


# --- search ------------------------------------------------------------------

def test_search_text_case_insensitive(client):
    r = client.get("/api/search", params=dict(q="TRANSPORTER"))
    body = r.json()
    hits = body["hits"]
    assert body["total"] == 3
    assert {h["orthogroup"] for h in hits} == {"OG_a", "OG_b", "OG_e"}
    assert all(h["field"] == "description" for h in hits)
    assert [h["position"] for h in hits] == sorted(h["position"] for h in hits)


def test_search_orthogroup_id(client):
    hits = client.get("/api/search", params=dict(q="og_a")).json()["hits"]
    assert hits[0]["field"] == "orthogroup"
    assert hits[0]["orthogroup"] == "OG_a"


@pytest.fixture(scope="module")
def many_hits(tmp_path_factory):
    d = tmp_path_factory.mktemp("hits")
    n = 600
    matrix = "Orthogroup\tA\tB\n" + "".join(f"OG{i:03d}\tx\ty\n" for i in range(n))
    rows = ["Orthogroup\tnote"] + [f"OG{i:03d}\tcommonword" for i in range(n)]
    out = _import(
        d, matrix=matrix, tree="((A:1,B:1)N0);",
        extras=["--annotation", str(_write(d, "ann.tsv", "\n".join(rows) + "\n"))])
    return TestClient(create_app(out))


def test_search_limit_cap(many_hits):
    r = many_hits.get("/api/search", params=dict(q="commonword", limit=1000))
    assert r.status_code == 200
    body = r.json()
    assert body["total"] == 600
    assert len(body["hits"]) == 500  # capped below the 600 hits


def test_search_order_tree_stable(client):
    body = client.get("/api/search", params=dict(q="transporter")).json()
    hits = sorted(body["hits"], key=lambda h: h["position"])
    assert body["total"] == 3
    oid = client.post("/api/orders", json={"by": "search", "q": "transporter"}).json()["order"]
    all_ids = client.get("/api/orthogroups", params=dict(order=oid, og0=0, og1=7)).json()
    tree_ids = client.get("/api/orthogroups", params=dict(order="tree", og0=0, og1=7)).json()
    hit_ids = [h["orthogroup"] for h in hits]
    # hits first, in their tree-relative order, then the non-hits in tree order
    assert all_ids[:3] == hit_ids
    assert all_ids[3:] == [og for og in tree_ids if og not in set(hit_ids)]
    tree_pos = {og: i for i, og in enumerate(tree_ids)}
    assert hit_ids == sorted(hit_ids, key=lambda o: tree_pos[o])
    assert [h["position"] for h in hits] == [tree_pos[o] for o in hit_ids]


# --- annotation track --------------------------------------------------------

TIE_MATRIX = (
    "Orthogroup\tA\tB\n"
    "OG0\tx\ty\nOG1\tx\ty\nOG2\tx\ty\nOG3\tx\ty\nOG4\tx\ty\n"
)
TIE_ANN = "Orthogroup\tcat\nOG0\tA\nOG1\tA\nOG2\tB\nOG3\tB\nOG4\t\n"


@pytest.fixture(scope="module")
def tie_client(tmp_path_factory):
    d = tmp_path_factory.mktemp("tie")
    out = _import(
        d, matrix=TIE_MATRIX, tree="((A:1,B:1)N0);",
        extras=["--annotation", str(_write(d, "cat.tsv", TIE_ANN))])
    return TestClient(create_app(out))


def test_annotation_track_tie_breaks_lower_index(tie_client):
    # bin 0: A,A -> A; bin 1: B,B,missing -> B; full: A=A,B=B tie -> A
    r = tie_client.get("/api/annotation-track", params=dict(
        order="tree", og0=0, og1=5, w=2, name="cat"))
    assert r.status_code == 200
    assert np.frombuffer(r.content, "<u2").tolist() == [1, 2]
    full = tie_client.get("/api/annotation-track", params=dict(
        order="tree", og0=0, og1=5, w=1, name="cat"))
    assert np.frombuffer(full.content, "<u2").tolist() == [1]
    # missing category -> 0
    tail = tie_client.get("/api/annotation-track", params=dict(
        order="tree", og0=4, og1=5, w=1, name="cat"))
    assert np.frombuffer(tail.content, "<u2").tolist() == [0]


def test_annotation_track_text_rejected(client):
    r = client.get("/api/annotation-track", params=dict(
        order="tree", og0=0, og1=7, w=1, name="description"))
    assert r.status_code == 400
    assert client.get("/api/annotation-track", params=dict(
        order="tree", og0=0, og1=7, w=1, name="nope")).status_code == 400


# --- nullable annotations / search robustness --------------------------------

def test_search_nullable_annotations(tmp_path):
    matrix = "Orthogroup\tA\tB\nOG_a\tx\ty\nOG_b\tx\ty\nOG_c\tx\ty\n"
    sparse = "Orthogroup\tnote\nOG_a\talpha\nOG_b\t\nOG_c\tbeta\n"
    allnull = "Orthogroup\temptycol\nOG_a\t\nOG_b\t\nOG_c\t\n"
    out = _import(tmp_path, matrix=matrix, tree="((A:1,B:1)N0);", extras=[
        "--annotation", str(_write(tmp_path, "sparse.tsv", sparse)),
        "--annotation", str(_write(tmp_path, "allnull.tsv", allnull)),
    ])
    ann = pl.read_parquet(out / "annotations.parquet")
    assert ann["emptycol"].dtype == pl.Utf8  # not inferred as Null
    c = TestClient(create_app(out))
    # GET search over sparse + all-null columns
    r = c.get("/api/search", params=dict(q="alpha"))
    assert r.status_code == 200
    assert [h["orthogroup"] for h in r.json()["hits"]] == ["OG_a"]
    assert c.get("/api/search", params=dict(q="beta")).status_code == 200
    # POST search order must not 500 on nulls
    r2 = c.post("/api/orders", json={"by": "search", "q": "beta"})
    assert r2.status_code == 200
    oid = r2.json()["order"]
    assert c.get("/api/orthogroups", params=dict(order=oid, og0=0, og1=1)).json() == ["OG_c"]


def test_annotation_empty_header_name(tmp_path):
    text = "Orthogroup\t\tcog\nOG_a\tX\n"
    with pytest.raises(ValueError, match="annotation column name is empty"):
        _import(tmp_path, extras=["--annotation", str(_write(tmp_path, "ann.tsv", text))])


# --- format version guard ----------------------------------------------------

def test_outdated_format_version_fails(tmp_path):
    d = tmp_path / "old"
    d.mkdir()
    (d / "manifest.json").write_text(json.dumps({"format_version": 1}))
    with pytest.raises(ValueError, match="format 1 is outdated; re-run orthopan import"):
        create_app(d)


# --- UPGMA honours --strip-prefix -------------------------------------------

def test_upgma_respects_custom_strip_prefix(tmp_path):
    matrix = "HOG1\tcustom|A_1, custom|B_2\nHOG2\tcustom|A_3\n"
    out = _import(tmp_path, matrix=matrix, tree=None,
                  extras=["--strip-prefix", r"^custom\|"], name="sp1")
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["tree_source"] == "upgma"
    genomes = pl.read_parquet(out / "genomes.parquet")
    assert set(genomes["name"].to_list()) == {"A", "B"}

    out2 = _import(tmp_path, matrix=matrix, tree=None,
                   extras=["--strip-prefix", ""], name="sp2")
    genomes2 = pl.read_parquet(out2 / "genomes.parquet")
    assert set(genomes2["name"].to_list()) == {"custom|A", "custom|B"}


# --- non-finite genome metadata ---------------------------------------------

def test_genome_metadata_nonfinite_is_missing(tmp_path):
    gm = "genome\tgc\nA\tnan\nB\tinf\nC\t0.5\nD\t-0.25\n"
    out = _import(tmp_path, extras=["--genome-metadata", str(_write(tmp_path, "gm.tsv", gm))])
    specs = json.loads((out / "manifest.json").read_text())["genome_metadata"]
    assert specs == [{"name": "gc", "kind": "numeric", "min": -0.25, "max": 0.5}]
    frame = pl.read_parquet(out / "genome_metadata.parquet").sort("genome_idx")
    genomes = pl.read_parquet(out / "genomes.parquet").sort("genome_idx")
    by = dict(zip(genomes["name"].to_list(), frame["gc"].to_list()))
    assert by["A"] is None and by["B"] is None
    assert by["C"] == 0.5


# --- metadata column-name validation ----------------------------------------

def test_metadata_reserved_name_rejected(tmp_path):
    gm = "genome\tname\nA\tX\nB\tY\nC\tX\nD\tY\n"
    with pytest.raises(ValueError, match="clash with built-in"):
        _import(tmp_path, extras=["--genome-metadata", str(_write(tmp_path, "gm.tsv", gm))])


def test_metadata_duplicate_name_rejected(tmp_path):
    gm = "genome\tclade\tclade\nA\tX\tY\n"
    with pytest.raises(ValueError, match="duplicate genome metadata column name"):
        _import(tmp_path, extras=["--genome-metadata", str(_write(tmp_path, "gm.tsv", gm))])


def test_metadata_empty_name_rejected(tmp_path):
    gm = "genome\t\tclade\nA\tX\tY\n"
    with pytest.raises(ValueError, match="genome metadata column name is empty"):
        _import(tmp_path, extras=["--genome-metadata", str(_write(tmp_path, "gm.tsv", gm))])


# --- end-to-end OrthoFinder import ------------------------------------------

@pytest.mark.skipif(not DEMO.is_dir(), reason="demo results not available")
def test_orthofinder_end_to_end(tmp_path):
    out = tmp_path / "demo"
    track = _write(tmp_path, "track.tsv", "orthogroup\tcustom_score\nOG0000000\t1.5\n")
    main([
        "import", "--orthofinder", str(DEMO), "-o", str(out),
        "--tree", str(DEMO / "Species_Tree" / "SpeciesTree_rooted.txt"),
        "--track", str(track),
    ])
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["tree_source"] == "input"  # explicit --tree overrides
    # explicit --track adds to the preprocessor columns
    assert {"duplications", "tree_length_per_leaf", "geometric_homogeneity",
            "custom_score"} <= set(manifest["track_columns"])
    assert manifest["detail"]["msa"]
    assert manifest["detail"]["gene_trees"]
    assert manifest["detail"]["sequences"]
    c = TestClient(create_app(out))
    meta = c.get("/api/meta").json()
    assert meta["detail"] == {
        "msa": True, "gene_trees": True, "sequences": True, "famsa": False}
    assert any(t["name"] == "geometric_homogeneity" for t in meta["tracks"])
