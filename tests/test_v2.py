"""Tests for the v2 single-orthogroup detail view (docs/api.md)."""

import concurrent.futures
import json
import os
import shutil
import subprocess
import threading
import time
from pathlib import Path
from urllib.parse import quote

import pytest
from fastapi import HTTPException
from fastapi.testclient import TestClient

from orthopan import main
from orthopan.importer import _parse_newick
from orthopan.io import read_matrix
from orthopan.server import _content_disposition, _famsa_ok, _to_newick, create_app

REPO = Path(__file__).resolve().parents[1]
DEMO = REPO / "resources/orthofinder-demo/results/Results_Oct04_1"
FAMSA = (REPO / "scratch/.venv/lib/python3.12/site-packages"
         / "orthofinder/bin/famsa")
requires_famsa = pytest.mark.skipif(not FAMSA.is_file(), reason="famsa binary not present")


def _write(path, text, binary=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.write_bytes if binary else path.write_text)(
        text.encode() if binary else text)
    return path


def _fake_orthofinder(root: Path):
    """A minimal OrthoFinder results dir with tricky genome and OG names."""
    _write(root / "Orthogroups/Orthogroups.tsv",
           "Orthogroup\tG.1\tG.2\tA_B\n"
           "OG1\tgnl|extdb|G.1_alpha\tgnl|extdb|G.2_beta\tgnl|extdb|A_B_gamma\n"
           "OG2\tgnl|extdb|G.1_delta\t\t\n"
           "OG\u6f22\tgnl|extdb|G.1_gi\t\t\n"
           "OG\"x\tgnl|extdb|G.1_gj\t\t\n")
    _write(root / "Species_Tree/SpeciesTree_rooted.txt",
           "((G.1:0.1,G.2:0.2)n1:0.3,A_B:0.4)root;")
    _write(root / "MultipleSequenceAlignments/OG1.fa",
           ">A_B_gnl|extdb|A_B_gamma\nMK-\n"
           ">G_1_gnl|extdb|G.1_alpha\nMKT\n"
           ">G_2_gnl|extdb|G.2_beta\nMKS\n")
    _write(root / "Orthogroup_Sequences/OG1.fa",
           ">gnl|extdb|G.1_alpha\nMKT\n"
           ">gnl|extdb|G.2_beta\nMKS\n"
           ">gnl|extdb|A_B_gamma\nMK\n")
    _write(root / "Orthogroup_Sequences/OG\u6f22.fa", ">gnl|extdb|G.1_gi\nMK\n")
    _write(root / 'Orthogroup_Sequences/OG"x.fa', ">gnl|extdb|G.1_gj\nMK\n")
    _write(root / "Resolved_Gene_Trees/Resolved_Gene_Trees.txt",
           "OG1: ((G_1_gnl|extdb|G.1_alpha:0.1,G_2_gnl|extdb|G.2_beta:0.2)n2:0.15,"
           "A_B_gnl|extdb|A_B_gamma:0.4)n1:0.3;\n"
           "OG\u6f22: G_1_gnl|extdb|G.1_gi:0.1;\n"
           "OG\"x: G_1_gnl|extdb|G.1_gj:0.1;\n")
    _write(root / "Gene_Duplication_Events/Duplications.tsv",
           "Orthogroup\tSpecies Tree Node\tGene Tree Node\tSupport\tType\tGenes 1\tGenes 2\n"
           "OG1\tn2\tn2\t0.9\tNon-Terminal\tx\ty\n"
           "OG1\tn1\tn1\t\tNon-Terminal\tx\ty\n")  # empty Support -> flag true, support null
    return root


@pytest.fixture(scope="module")
def fake(tmp_path_factory):
    out = tmp_path_factory.mktemp("v2fake") / "out"
    main(["import", "--orthofinder", str(_fake_orthofinder(tmp_path_factory.mktemp("of"))),
          "-o", str(out)])
    return out


@pytest.fixture(scope="module")
def fake_client(fake):
    return TestClient(create_app(fake))


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    if os.environ.get("ORTHOPAN_TEST_NO_RESOURCES") == "1" or not DEMO.is_dir():
        pytest.skip("resources/orthofinder-demo not available")
    out = tmp_path_factory.mktemp("v2demo") / "out"
    main(["import", "--orthofinder", str(DEMO), "-o", str(out)])
    return out


@pytest.fixture(scope="module")
def demo_client(demo):
    return TestClient(create_app(demo))


def _no_msa_og(demo):
    msa_dir = Path(json.loads((demo / "manifest.json").read_text())["detail"]["msa"])
    msa_stems = {p.stem for p in msa_dir.glob("*.fa")}
    for path in sorted((DEMO / "Orthogroup_Sequences").glob("*.fa")):
        if path.stem not in msa_stems:
            return path.stem
    raise AssertionError("no OG without an MSA")


# --- header mapping, row order, tricky genome names -------------------------

def test_msa_header_mapping_and_row_order(fake_client):
    body = fake_client.get("/api/orthogroup/OG1/msa").json()
    assert body["source"] == "orthofinder"
    assert body["length"] == 3
    keys = [(r["genome"], r["row"], r["gene"], r["raw"]) for r in body["rows"]]
    # sorted by genome row (G.1=0, G.2=1, A_B=2), then gene
    assert keys == [
        ("G.1", 0, "G.1_alpha", "gnl|extdb|G.1_alpha"),
        ("G.2", 1, "G.2_beta", "gnl|extdb|G.2_beta"),
        ("A_B", 2, "A_B_gamma", "gnl|extdb|A_B_gamma"),
    ]


def test_msa_unknown_orthogroup_404(fake_client):
    assert fake_client.get("/api/orthogroup/NOPE/msa").status_code == 404


# --- gene tree ---------------------------------------------------------------

def test_gene_tree_leaf_mapping_and_duplications(fake_client):
    body = fake_client.get("/api/orthogroup/OG1/gene-tree").json()
    assert body["parent"][0] == -1
    leaves = {
        body["label"][i]: (body["gene"][i], body["genome_row"][i])
        for i in range(len(body["label"]))
        if body["gene"][i] is not None
    }
    assert leaves == {
        "G_1_gnl|extdb|G.1_alpha": ("G.1_alpha", 0),
        "G_2_gnl|extdb|G.2_beta": ("G.2_beta", 1),
        "A_B_gnl|extdb|A_B_gamma": ("A_B_gamma", 2),
    }
    dup = {body["label"][i]: body["support"][i]
           for i in range(len(body["label"])) if body["duplication"][i]}
    # n1 has an empty Support cell: still a duplication, support stays null
    assert dup == {"n2": 0.9, "n1": None}


def test_gene_tree_missing_404(fake_client):
    assert fake_client.get("/api/orthogroup/OG2/gene-tree").status_code == 404


def test_gene_tree_real_duplications(demo_client):
    body = demo_client.get("/api/orthogroup/OG0000000/gene-tree").json()
    labels, rows = read_matrix(DEMO / "Gene_Duplication_Events" / "Duplications.tsv")
    node_col = labels.index("Gene Tree Node")
    support_col = labels.index("Support")
    expected = {}
    for og, cells in rows:
        if og == "OG0000000":
            expected[cells[node_col]] = float(cells[support_col])
    assert expected  # the OG really has duplication events
    got = {body["label"][i]: body["support"][i]
           for i in range(len(body["label"])) if body["duplication"][i]}
    for node, support in expected.items():
        assert got[node] == pytest.approx(support)


# --- downloads ---------------------------------------------------------------

def test_downloads_display_and_raw(fake_client):
    raw = fake_client.get("/api/orthogroup/OG1/fasta",
                          params=dict(aligned=0, ids="raw"))
    assert raw.status_code == 200
    assert 'filename="OG1.fa"' in raw.headers["content-disposition"]
    assert raw.text.splitlines()[0] == ">gnl|extdb|G.1_alpha"
    assert ">gnl|extdb|A_B_gamma" in raw.text

    display = fake_client.get("/api/orthogroup/OG1/fasta",
                              params=dict(aligned=0, ids="display"))
    assert display.text.splitlines()[0] == ">G.1_alpha"

    aligned = fake_client.get("/api/orthogroup/OG1/fasta",
                              params=dict(aligned=1, ids="display"))
    assert 'filename="OG1.aligned.fa"' in aligned.headers["content-disposition"]
    assert aligned.text.splitlines()[1] == "MKT"  # G.1 row first

    nwk = fake_client.get("/api/orthogroup/OG1/gene-tree.nwk", params=dict(ids="raw"))
    assert nwk.status_code == 200
    assert 'filename="OG1.nwk"' in nwk.headers["content-disposition"]
    assert "gnl|extdb|G.1_alpha" in nwk.text and nwk.text.endswith(";")

    nwk_disp = fake_client.get("/api/orthogroup/OG1/gene-tree.nwk",
                               params=dict(ids="display"))
    assert "G.1_alpha" in nwk_disp.text and "gnl|extdb" not in nwk_disp.text

    assert fake_client.get(
        "/api/orthogroup/OG1/fasta", params=dict(ids="bogus")).status_code == 400


def test_download_no_sequences_404(fake_client):
    assert fake_client.get("/api/orthogroup/OG2/fasta").status_code == 404


# --- famsa on-demand alignment ----------------------------------------------

@requires_famsa
def test_famsa_aligns_missing_msa_and_caches(demo):
    shutil.rmtree(demo / "cache", ignore_errors=True)
    og = _no_msa_og(demo)
    c = TestClient(create_app(demo, famsa=str(FAMSA)))
    body = c.get(f"/api/orthogroup/{og}/msa").json()
    assert body["source"] == "famsa"
    assert body["rows"]
    # every row maps back to a genome/membership entry, with a stripped display ID
    assert all(r["genome"] for r in body["rows"])
    assert all(r["row"] >= 0 for r in body["rows"])
    assert all("gnl|extdb" not in r["gene"] for r in body["rows"])
    assert [(r["row"], r["gene"]) for r in body["rows"]] == sorted(
        (r["row"], r["gene"]) for r in body["rows"])
    cached = demo / "cache" / "msa" / f"{og}.fa"
    assert cached.is_file()
    headers = {line[1:] for line in cached.read_text().splitlines()
               if line.startswith(">")}
    expected = {f"{r['genome'].replace('.', '_')}_{r['raw']}"
                for r in body["rows"] if r["genome"]}
    assert headers == expected  # cached with canonical genome_raw headers
    # second request is a cache hit and returns the same data
    assert c.get(f"/api/orthogroup/{og}/msa").json() == body


@requires_famsa
def test_famsa_concurrency_503(demo, tmp_path, monkeypatch):
    og = _no_msa_og(demo)
    app = create_app(demo, famsa=str(FAMSA))
    pg = app.state.pg
    monkeypatch.setattr(pg, "cache_dir", tmp_path / "cache")  # force a miss
    assert pg._famsa_sem.acquire(blocking=False)
    assert pg._famsa_sem.acquire(blocking=False)
    c = TestClient(app)
    assert c.get(f"/api/orthogroup/{og}/msa").status_code == 503
    pg._famsa_sem.release()
    pg._famsa_sem.release()
    assert c.get(f"/api/orthogroup/{og}/msa").status_code == 200


def test_famsa_unavailable_404(demo):
    shutil.rmtree(demo / "cache", ignore_errors=True)
    og = _no_msa_og(demo)
    c = TestClient(create_app(demo, famsa="definitely-not-a-famsa-binary"))
    r = c.get(f"/api/orthogroup/{og}/msa")
    assert r.status_code == 404
    assert "famsa is not available" in r.json()["detail"]
    assert c.get("/api/meta").json()["detail"]["famsa"] is False


def test_too_many_sequences_404(demo, monkeypatch):
    import orthopan.server as server
    og = _no_msa_og(demo)
    monkeypatch.setattr(server, "_FAMSA_MAX_SEQS", 0)
    c = TestClient(create_app(demo, famsa=str(FAMSA) if FAMSA.is_file() else "famsa"))
    r = c.get(f"/api/orthogroup/{og}/msa")
    assert r.status_code == 404
    assert "exceed" in r.json()["detail"]


# --- --fasta index -----------------------------------------------------------

def test_fasta_index_duplicate_ids_and_gz_policy(tmp_path, capsys):
    matrix = "Orthogroup\tA\tB\nOG1\tgnl|x|g1\tgnl|x|g1\n"
    (tmp_path / "matrix.tsv").write_text(matrix)
    (tmp_path / "tree.nwk").write_text("((A:1,B:1)N0);")
    fa_dir = tmp_path / "fa"
    _write(fa_dir / "A.faa", ">gnl|x|g1\nAAAAAA\n")
    _write(fa_dir / "B.faa", ">gnl|x|g1\nAAAA\n")  # same ID, different length
    _write(fa_dir / "C.faa.gz", ">gnl|x|other\nAAAA\n")  # gz -> skipped
    out = tmp_path / "out"
    main(["import", str(tmp_path / "matrix.tsv"), "-o", str(out),
          "--tree", str(tmp_path / "tree.nwk"), "--fasta", str(fa_dir)])
    err = capsys.readouterr().err
    assert "skipping gzip FASTA" in err

    manifest = json.loads((out / "manifest.json").read_text())
    assert len(manifest["fasta"]) == 2  # A.faa, B.faa; C skipped
    c = TestClient(create_app(out))
    body = c.get("/api/orthogroup/OG1/fasta", params=dict(ids="display")).text
    lines = body.splitlines()
    # A (row 0) resolved to the longer A.faa sequence, B to the shorter B.faa
    assert lines[0] == ">g1" and lines[1] == "AAAAAA"
    assert lines[2] == ">g1" and lines[3] == "AAAA"
    assert c.get("/api/meta").json()["detail"]["sequences"] is True


def test_fasta_prefix_stripped_headers(tmp_path):
    # FASTA headers already stripped; membership raw still matches via display.
    matrix = "Orthogroup\tA\tB\nOG1\tgnl|x|g1\t\n"
    (tmp_path / "matrix.tsv").write_text(matrix)
    (tmp_path / "tree.nwk").write_text("((A:1,B:1)N0);")
    _write(tmp_path / "fa" / "A.faa", ">g1\nAAAA\n")
    out = tmp_path / "out"
    main(["import", str(tmp_path / "matrix.tsv"), "-o", str(out),
          "--tree", str(tmp_path / "tree.nwk"), "--fasta", str(tmp_path / "fa")])
    c = TestClient(create_app(out))
    body = c.get("/api/orthogroup/OG1/fasta", params=dict(ids="raw")).text
    assert body.splitlines()[0] == ">gnl|x|g1"
    assert body.splitlines()[1] == "AAAA"


def test_reimport_without_fasta_drops_stale_sequences(tmp_path):
    matrix = "Orthogroup\tA\tB\nOG1\tgnl|x|g1\t\n"
    (tmp_path / "matrix.tsv").write_text(matrix)
    (tmp_path / "tree.nwk").write_text("((A:1,B:1)N0);")
    _write(tmp_path / "fa" / "A.faa", ">gnl|x|g1\nAAAA\n")
    out = tmp_path / "out"
    main(["import", str(tmp_path / "matrix.tsv"), "-o", str(out),
          "--tree", str(tmp_path / "tree.nwk"), "--fasta", str(tmp_path / "fa")])
    assert (out / "sequences.parquet").is_file()
    # a reimport into the same directory without --fasta must not leave it behind
    main(["import", str(tmp_path / "matrix.tsv"), "-o", str(out),
          "--tree", str(tmp_path / "tree.nwk")])
    assert not (out / "sequences.parquet").exists()


def test_fasta_equal_length_duplicates(tmp_path):
    def build(name, fa_names, seq_a, seq_b):
        d = tmp_path / name
        _write(d / f"{fa_names[0]}.faa", f">gnl|x|g1\n{seq_a}\n")
        _write(d / f"{fa_names[1]}.faa", f">gnl|x|g1\n{seq_b}\n")
        (d / "matrix.tsv").write_text("Orthogroup\tA\tB\nOG1\tgnl|x|g1\t\n")
        (d / "tree.nwk").write_text("((A:1,B:1)N0);")
        out = d / "out"
        main(["import", str(d / "matrix.tsv"), "-o", str(out),
              "--tree", str(d / "tree.nwk"), "--fasta", str(d)])
        return out

    # unique genome/file-stem match wins even with different equal-length content
    out = build("stem_match", ("A", "B"), "AAAAAA", "CCCCCC")
    body = TestClient(create_app(out)).get("/api/orthogroup/OG1/fasta").text
    assert body.splitlines()[1] == "AAAAAA"

    # no stem match + equal length + different content -> ambiguous, unresolved
    out2 = build("ambiguous", ("X", "Y"), "AAAAAA", "CCCCCC")
    assert not (out2 / "sequences.parquet").exists()
    assert TestClient(create_app(out2)).get("/api/orthogroup/OG1/fasta").status_code == 404

    # no stem match but identical sequences -> first wins
    out3 = build("identical", ("X", "Y"), "AAAAAA", "AAAAAA")
    body3 = TestClient(create_app(out3)).get("/api/orthogroup/OG1/fasta").text
    assert body3.splitlines()[1] == "AAAAAA"


def _script(path, body):
    path.write_text("#!/bin/sh\n" + body + "\n")
    path.chmod(0o755)
    return path


def test_famsa_probe_validates_output(tmp_path):
    good = _script(tmp_path / "good.sh", "printf '>a\\nAAAA\\n>b\\nAAAT\\n' > \"$2\"")
    bad = _script(tmp_path / "bad.sh", "printf 'not a fasta' > \"$2\"")
    empty = _script(tmp_path / "empty.sh", "exit 0")
    assert _famsa_ok(str(good)) is True
    assert _famsa_ok(str(bad)) is False
    assert _famsa_ok(str(empty)) is False


@pytest.mark.parametrize("mode", ["empty", "wrong_count", "unequal", "empty_seq", "duplicate_header"])
def test_famsa_invalid_output_not_cached(demo, tmp_path, monkeypatch, mode):
    og = _no_msa_og(demo)
    app = create_app(demo, famsa=str(FAMSA) if FAMSA.is_file() else "famsa")
    pg = app.state.pg
    monkeypatch.setattr(pg, "cache_dir", tmp_path / f"cache_{mode}")

    def runner(og_name, seqs):
        n = len(seqs)
        if mode == "empty":
            return ""
        if mode == "wrong_count":
            return "".join(f">m{i}\nA\n" for i in range(max(0, n - 1)))
        if mode == "duplicate_header":  # right count, m0 twice, last member missing
            return "".join(f">m{min(i, n - 2)}\nA\n" for i in range(n))
        if mode == "unequal":
            return "".join(f">m{i}\n{'A' * (i + 1)}\n" for i in range(n))
        return "".join(f">m{i}\n\n" for i in range(n))

    monkeypatch.setattr(pg, "_run_famsa", runner)
    r = TestClient(app).get(f"/api/orthogroup/{og}/msa")
    assert r.status_code == 500
    assert not (pg.cache_dir / f"{og}.fa").exists()
    assert not list(pg.cache_dir.glob("*")) if pg.cache_dir.exists() else True


@requires_famsa
def test_famsa_same_og_single_run(demo, tmp_path, monkeypatch):
    og = _no_msa_og(demo)
    app = create_app(demo, famsa=str(FAMSA))
    pg = app.state.pg
    monkeypatch.setattr(pg, "cache_dir", tmp_path / "cache")
    calls = []

    def runner(og_name, seqs):
        calls.append(og_name)
        time.sleep(0.3)
        return "".join(f">m{i}\n{'A' * 8}\n" for i in range(len(seqs)))

    monkeypatch.setattr(pg, "_run_famsa", runner)
    c = TestClient(app)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        statuses = list(ex.map(
            lambda _: c.get(f"/api/orthogroup/{og}/msa").status_code, range(2)))
    assert statuses == [200, 200]
    assert calls == [og]  # the per-OG lock + cache recheck align only once


@requires_famsa
def test_famsa_timeout_and_failure_cleanup(demo, tmp_path, monkeypatch):
    og = _no_msa_og(demo)
    app = create_app(demo, famsa=str(FAMSA))
    pg = app.state.pg
    monkeypatch.setattr(pg, "cache_dir", tmp_path / "cache")
    c = TestClient(app)
    errors = [
        subprocess.TimeoutExpired("famsa", 1),
        subprocess.CalledProcessError(1, "famsa", b"", b"boom"),
    ]
    for exc in errors:
        def boom(*args, **kwargs):
            raise exc
        monkeypatch.setattr(subprocess, "run", boom)
        r = c.get(f"/api/orthogroup/{og}/msa")
        assert r.status_code == 500
        assert not (pg.cache_dir / f"{og}.fa").exists()
        if pg.cache_dir.exists():
            assert list(pg.cache_dir.iterdir()) == []  # temp files cleaned up
        assert pg._famsa_sem._value == 2  # slot returned after the failure
    # semaphore still usable: the next attempt fails again with 500, not 503
    assert c.get(f"/api/orthogroup/{og}/msa").status_code == 500


# --- newick labels, filename encoding ---------------------------------------

def test_newick_round_trip_preserves_internal_labels():
    root = _parse_newick("((A:1,B:2)'n,1':0.5,C:1)'root x';")
    text = _to_newick(root, lambda node: node.label)
    labels = []
    stack = [_parse_newick(text)]
    while stack:
        node = stack.pop()
        labels.append(node.label)
        stack.extend(node.children)
    assert "n,1" in labels and "root x" in labels


def test_newick_endpoint_keeps_internal_labels(fake_client):
    text = fake_client.get("/api/orthogroup/OG1/gene-tree.nwk").text
    assert "n2" in text and "n1" in text
    _parse_newick(text)  # still valid Newick


def test_download_filename_encoding(fake_client):
    for og in ["OG\u6f22", 'OG"x']:
        for path, params in (("fasta", {"aligned": 0, "ids": "raw"}),
                             ("gene-tree.nwk", {"ids": "raw"})):
            r = fake_client.get(f"/api/orthogroup/{quote(og)}/{path}", params=params)
            assert r.status_code == 200, (og, path, r.text)
            cd = r.headers["content-disposition"]
            assert "filename*=UTF-8''" in cd
            assert 'OG"' not in cd  # quote escaped in the ASCII fallback
    with pytest.raises(HTTPException):
        _content_disposition("bad\x01name.fa")
