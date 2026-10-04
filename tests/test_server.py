"""Tests for the FastAPI backend (docs/api.md)."""

import numpy as np
import polars as pl
import pytest
from fastapi.testclient import TestClient

import orthopan.server as server
from orthopan import main
from orthopan.server import MAX_CELLS, create_app

TREE = "((A:1,B:1)N3:0.5,(C:1,D:1)N4:0.5)N0;"

# header order C,A,D,B differs from tree DFS order (A,B,C,D)
MATRIX = (
    "Orthogroup\tC\tA\tD\tB\n"
    "OG_core\tx\tx\tx\tx\n"
    "OG_ab\t\tx\t\tx\n"
    "OG_cd\tx\t\tx\t\n"
    "OG_ac\tx\tx\t\t\n"
    "OG_a\t\tx\t\t\n"
    "OG_c\tx\t\t\t\n"
    "OG_empty\t\t\t\t\n"
)
TRACK = (
    "orthogroup\tscore\n"
    "OG_core\t1\n"
    "OG_ab\t4\n"
    "OG_cd\t3\n"
    "OG_ac\t2\n"
    "OG_a\t5\n"
    "OG_c\t\n"
    "OG_empty\t\n"
)
OG_IDS = ["OG_core", "OG_ab", "OG_cd", "OG_ac", "OG_a", "OG_c", "OG_empty"]
OG_SCORES = np.array([1.0, 4.0, 3.0, 2.0, 5.0, np.nan, np.nan])

WIDE_N = 18  # > 17 OGs so the packed representation spans three bytes


def _wide_matrix(n=WIDE_N):
    labels = ["C", "A", "D", "B"]
    lines = ["Orthogroup\t" + "\t".join(labels)]
    for i in range(n):
        cells = ["x" if ((i * 3 + j * 5) % 7) < 4 else "" for j in range(len(labels))]
        lines.append(f"OG{i:02d}\t" + "\t".join(cells))
    return "\n".join(lines) + "\n"


def _wide_track(n=WIDE_N):
    lines = ["orthogroup\tscore"]
    for i in range(n):
        lines.append(f"OG{i:02d}\t{'' if i % 5 == 0 else (i % 9) + 1}")
    return "\n".join(lines) + "\n"


def _import(dirpath, matrix, tree, tracks=()):
    dirpath.mkdir(parents=True, exist_ok=True)
    (dirpath / "matrix.tsv").write_text(matrix)
    (dirpath / "tree.nwk").write_text(tree)
    argv = [str(dirpath / "matrix.tsv"), "-o", str(dirpath / "out"),
            "--tree", str(dirpath / "tree.nwk")]
    for i, text in enumerate(tracks):
        (dirpath / f"track{i}.tsv").write_text(text)
        argv += ["--track", str(dirpath / f"track{i}.tsv")]
    main(["import", *argv])
    return dirpath / "out"


@pytest.fixture(scope="module")
def tiny(tmp_path_factory):
    return _import(tmp_path_factory.mktemp("tiny"), MATRIX, TREE, [TRACK])


@pytest.fixture(scope="module")
def wide(tmp_path_factory):
    return _import(tmp_path_factory.mktemp("wide"), _wide_matrix(), TREE, [_wide_track()])


@pytest.fixture(scope="module")
def big_og(tmp_path_factory):
    n = server.MAX_LIST + 1
    lines = ["Orthogroup\tA\tB"]
    lines += [f"OG{i:05d}\tx\t" + ("x" if i % 2 else "") for i in range(n)]
    d = tmp_path_factory.mktemp("big_og")
    return _import(d, "\n".join(lines) + "\n", "((A:1,B:1)N0);")


@pytest.fixture(scope="module")
def client(tiny):
    return TestClient(create_app(tiny))


@pytest.fixture(scope="module")
def wide_client(wide):
    return TestClient(create_app(wide))


def _load_arrays(out):
    """Dense presence in tree row order and the built-in order permutations."""
    genomes = pl.read_parquet(out / "genomes.parquet").sort("genome_idx")
    og = pl.read_parquet(out / "orthogroups.parquet").sort("og_idx")
    n_og = og.height
    row_order = genomes["row_order"].to_numpy()
    inv = np.empty(len(row_order), np.int64)
    inv[row_order] = np.arange(len(row_order))
    packed = np.load(out / "presence.npy")
    dense = np.unpackbits(packed, axis=1, bitorder="big")[:, :n_og][inv]
    perms = {
        "tree": np.argsort(og["order_tree"].to_numpy(), kind="stable"),
        "frequency": np.argsort(og["order_frequency"].to_numpy(), kind="stable"),
    }
    return dense, perms


@pytest.fixture(scope="module")
def arrays(tiny):
    return _load_arrays(tiny)


@pytest.fixture(scope="module")
def wide_arrays(wide):
    return _load_arrays(wide)


def _read_scores(path):
    lines = path.read_text().strip().splitlines()[1:]
    return np.array(
        [np.nan if ln.split("\t")[1] == "" else float(ln.split("\t")[1]) for ln in lines]
    )


def _track_perm(values):
    missing = np.isnan(values)
    present = np.flatnonzero(~missing)
    ranked = np.argsort(values[present], kind="stable")
    return np.concatenate([present[ranked], np.flatnonzero(missing)]).astype(np.int64)


def _brute(dense_order, og0, og1, g0, g1, w, h):
    col_edges = og0 + (np.arange(w + 1) * (og1 - og0)) // w
    row_edges = g0 + (np.arange(h + 1) * (g1 - g0)) // h
    out = np.zeros((h, w), np.uint8)
    for i in range(h):
        for j in range(w):
            block = dense_order[row_edges[i]:row_edges[i + 1], col_edges[j]:col_edges[j + 1]]
            out[i, j] = int(np.rint(block.mean() * 255))
    return out


def test_meta(client):
    body = client.get("/api/meta").json()
    assert body["n_genomes"] == 4
    assert body["n_orthogroups"] == 7
    assert body["tree_source"] == "input"
    assert body["orders"] == ["tree", "frequency"]
    by_name = {t["name"]: t for t in body["tracks"]}
    assert by_name["score"] == {"name": "score", "min": 1.0, "max": 5.0, "missing": 2}
    assert by_name["copy_number_mean"]["missing"] == 1  # OG_empty


def test_genomes_row_order(client, tiny):
    rows = client.get("/api/genomes").json()
    assert [r["row"] for r in rows] == [0, 1, 2, 3]
    assert [r["name"] for r in rows] == ["A", "B", "C", "D"]
    genomes = pl.read_parquet(tiny / "genomes.parquet").sort("genome_idx")
    by_name = {g["name"]: g for g in genomes.iter_rows(named=True)}
    for r in rows:
        g = by_name[r["name"]]
        assert r["genome_idx"] == int(g["genome_idx"])
        assert r["n_orthogroups"] == int(g["n_orthogroups"])
        assert r["n_genes"] == int(g["n_genes"])


def test_tree_arrays_consistent(client, tiny):
    t = client.get("/api/tree").json()
    parent, length, label, row = t["parent"], t["length"], t["label"], t["row"]
    assert len(parent) == len(length) == len(label) == len(row)
    assert parent[0] == -1
    for i in range(1, len(parent)):
        assert 0 <= parent[i] < i  # parents precede children
    leaf_rows = sorted(r for r in row if r >= 0)
    assert leaf_rows == list(range(len(leaf_rows)))  # leaves cover all rows
    assert len(leaf_rows) == 4
    # each genome's leaf label sits at its row_order
    row_label = {r: label[i] for i, r in enumerate(row) if r >= 0}
    genomes = pl.read_parquet(tiny / "genomes.parquet")
    for g in genomes.iter_rows(named=True):
        assert row_label[int(g["row_order"])] == g["name"]


def test_matrix_full_detail_matches_presence(client, arrays):
    dense, perms = arrays
    n_og = dense.shape[1]
    for order, perm in list(perms.items()) + [("track:score", _track_perm(OG_SCORES))]:
        r = client.get("/api/matrix", params=dict(
            order=order, og0=0, og1=n_og, g0=0, g1=dense.shape[0], w=n_og, h=dense.shape[0]))
        assert r.status_code == 200
        assert r.headers["X-Width"] == str(n_og)
        assert r.headers["X-Height"] == str(dense.shape[0])
        got = np.frombuffer(r.content, "<u1").reshape(dense.shape[0], n_og)
        assert np.array_equal(got, dense[:, perm].astype(np.uint8) * 255)


def test_matrix_full_detail_wide_crosses_bytes(wide, wide_client, wide_arrays):
    dense, perms = wide_arrays
    n_og = dense.shape[1]
    assert n_og >= 17  # packed representation is three bytes wide
    scores = _read_scores(wide.parent / "track0.tsv")
    for order, perm in list(perms.items()) + [("track:score", _track_perm(scores))]:
        r = wide_client.get("/api/matrix", params=dict(
            order=order, og0=0, og1=n_og, g0=0, g1=4, w=n_og, h=4))
        got = np.frombuffer(r.content, "<u1").reshape(4, n_og)
        assert np.array_equal(got, dense[:, perm].astype(np.uint8) * 255)


def test_matrix_aggregation_matches_bruteforce(client, arrays):
    dense, perms = arrays
    for order, perm in perms.items():
        ordered = dense[:, perm]
        for og0, og1, w in [(0, 7, 3), (1, 7, 3), (0, 7, 7), (2, 5, 2)]:
            for g0, g1, h in [(0, 4, 2), (1, 4, 3), (0, 4, 4)]:
                r = client.get("/api/matrix", params=dict(
                    order=order, og0=og0, og1=og1, g0=g0, g1=g1, w=w, h=h))
                got = np.frombuffer(r.content, "<u1").reshape(h, w)
                assert np.array_equal(got, _brute(ordered, og0, og1, g0, g1, w, h)), \
                    (order, og0, og1, w, g0, g1, h)


def test_matrix_aggregation_wide_matches_bruteforce(wide_arrays, wide_client):
    dense, perms = wide_arrays
    n_og = dense.shape[1]
    for order, perm in perms.items():
        ordered = dense[:, perm]
        for og0, og1, w in [(0, n_og, 5), (1, n_og, 7), (3, n_og - 1, 3), (0, n_og, n_og)]:
            for g0, g1, h in [(0, 4, 2), (1, 4, 3), (0, 4, 4)]:
                r = wide_client.get("/api/matrix", params=dict(
                    order=order, og0=og0, og1=og1, g0=g0, g1=g1, w=w, h=h))
                got = np.frombuffer(r.content, "<u1").reshape(h, w)
                assert np.array_equal(got, _brute(ordered, og0, og1, g0, g1, w, h)), \
                    (order, og0, og1, w, g0, g1, h)


def test_tracks_mean_ignores_nan_and_all_nan_bin(client, arrays):
    # full detail over the track order: each bin is one OG, last two are missing.
    r = client.get("/api/tracks", params=dict(order="track:score", og0=0, og1=7, w=7,
                                              names="score"))
    assert r.status_code == 200
    assert r.headers["X-Tracks"] == "score"
    got = np.frombuffer(r.content, "<f4")
    assert np.allclose(got[:5], [1.0, 2.0, 3.0, 4.0, 5.0])
    assert np.isnan(got[5]) and np.isnan(got[6])

    # a bin mixing present and missing averages only the present values
    dense, perms = arrays
    r = client.get("/api/tracks", params=dict(order="tree", og0=0, og1=7, w=3,
                                              names="score"))
    got = np.frombuffer(r.content, "<f4")
    ordered = OG_SCORES[perms["tree"]]
    for j in range(3):
        lo = (j * 7) // 3
        hi = ((j + 1) * 7) // 3
        vals = ordered[lo:hi]
        expected = np.nan if np.all(np.isnan(vals)) else np.nanmean(vals)
        assert (np.isnan(got[j]) and np.isnan(expected)) or got[j] == pytest.approx(expected)


def test_tracks_multiple_names_track_major(client):
    r = client.get("/api/tracks", params=dict(order="tree", og0=0, og1=7, w=3,
                                              names="score,copy_number_mean"))
    assert r.headers["X-Tracks"] == "score,copy_number_mean"
    got = np.frombuffer(r.content, "<f4")
    assert got.shape == (6,)
    first = np.frombuffer(
        client.get("/api/tracks", params=dict(order="tree", og0=0, og1=7, w=3,
                                              names="score")).content, "<f4")
    assert np.allclose(got[:3], first)


def test_orthogroups_list(client):
    # default order is tree: DFS ranks [0,2,4,1,3,5,6] -> argsort
    assert client.get("/api/orthogroups", params=dict(og0=0, og1=7)).json() == [
        "OG_core", "OG_ac", "OG_ab", "OG_a", "OG_cd", "OG_c", "OG_empty"]
    sub = client.get("/api/orthogroups", params=dict(order="frequency", og0=0, og1=3)).json()
    assert sub == ["OG_core", "OG_ab", "OG_cd"]


def test_orthogroups_limit_exact(big_og):
    c = TestClient(create_app(big_og))
    assert len(c.get("/api/orthogroups", params=dict(og0=0, og1=5000)).json()) == 5000
    assert c.get("/api/orthogroups", params=dict(og0=0, og1=5001)).status_code == 400


def test_orthogroup_detail_members(client):
    body = client.get("/api/orthogroup/OG_a").json()
    assert body["orthogroup"] == "OG_a"
    assert body["og_idx"] == 4
    assert body["n_genomes"] == 1
    assert body["n_genes"] == 1
    assert body["tracks"]["score"] == 5.0
    assert body["members"] == [{"genome": "A", "genes": [{"id": "x", "raw": "x"}]}]

    empty = client.get("/api/orthogroup/OG_empty").json()
    assert empty["n_genomes"] == 0
    assert empty["members"] == []


def test_quoted_punctuation_genome_labels(tmp_path):
    labels = ["a b", "c,d's"]
    matrix = "Orthogroup\t" + "\t".join(labels) + "\nOG1\tx\tx\nOG2\tx\t\n"
    out = _import(tmp_path, matrix, "('a b':1,'c,d''s':1)N0;")
    c = TestClient(create_app(out))
    assert {g["name"] for g in c.get("/api/genomes").json()} == set(labels)
    t = c.get("/api/tree").json()
    row_label = {r: t["label"][i] for i, r in enumerate(t["row"]) if r >= 0}
    for g in pl.read_parquet(out / "genomes.parquet").iter_rows(named=True):
        assert row_label[int(g["row_order"])] == g["name"]
    r = c.get("/api/matrix", params=dict(order="tree", og0=0, og1=2, g0=0, g1=2, w=2, h=2))
    assert r.status_code == 200


def test_contrast_order(client):
    r = client.post("/api/orders", json={"by": "contrast", "genomes": ["A", "B"]})
    assert r.status_code == 200
    oid = r.json()["order"]
    ids = client.get("/api/orthogroups", params=dict(order=oid, og0=0, og1=7)).json()
    # ab and a only in selection (score 1); cd only outside (-1); c has one hit
    # outside so scores -0.5.
    assert ids == ["OG_ab", "OG_a", "OG_core", "OG_ac", "OG_empty", "OG_c", "OG_cd"]
    # the same id resolves for the matrix
    m = client.get("/api/matrix", params=dict(order=oid, og0=0, og1=7, g0=0, g1=4,
                                              w=7, h=4))
    assert m.status_code == 200


def test_order_id_no_collision(tmp_path):
    # comma-joining ["A","B,C"] and ["A,B","C"] would collide; json.dumps does not.
    labels = ["A", "B,C", "A,B", "C"]
    matrix = (
        "Orthogroup\t" + "\t".join(labels) + "\n"
        "OG1\tx\tx\t\t\n"
        "OG2\t\t\tx\tx\n"
        "OG3\tx\t\tx\t\n"
    )
    out = _import(tmp_path, matrix, "((A:1,'B,C':1)N1:0.5,('A,B':1,C:1)N2:0.5)N0;")
    c = TestClient(create_app(out))
    r1 = c.post("/api/orders", json={"by": "contrast", "genomes": ["A", "B,C"]})
    r2 = c.post("/api/orders", json={"by": "contrast", "genomes": ["A,B", "C"]})
    assert r1.status_code == r2.status_code == 200
    id1, id2 = r1.json()["order"], r2.json()["order"]
    assert id1 != id2
    before = c.get("/api/orthogroups", params=dict(order=id1, og0=0, og1=3)).json()
    # creating id2 must not overwrite id1's cached permutation
    after = c.get("/api/orthogroups", params=dict(order=id1, og0=0, og1=3)).json()
    assert before == after
    assert c.get("/api/matrix", params=dict(order=id1, og0=0, og1=3, g0=0, g1=4,
                                            w=3, h=4)).status_code == 200


def test_dynamic_lru_eviction_consistent(tiny, monkeypatch):
    monkeypatch.setattr(server, "_LRU_SIZE", 1)
    c = TestClient(create_app(tiny))
    id1 = c.post("/api/orders", json={"by": "contrast", "genomes": ["A", "B"]}).json()["order"]
    id2 = c.post("/api/orders", json={"by": "contrast", "genomes": ["C", "D"]}).json()["order"]
    assert c.get("/api/orthogroups", params=dict(order=id2, og0=0, og1=7)).status_code == 200
    # evicted order is gone for labels and matrix alike (no 200/404 skew)
    assert c.get("/api/orthogroups", params=dict(order=id1, og0=0, og1=7)).status_code == 404
    assert c.get("/api/matrix", params=dict(order=id1, og0=0, og1=7, g0=0, g1=4,
                                            w=7, h=4)).status_code == 404


def _distinct_selections():
    return [["A", "B"], ["A", "C"], ["A", "D"], ["B", "C"], ["B", "D"], ["C", "D"],
            ["A", "B", "C"], ["A", "B", "D"], ["A", "C", "D"], ["B", "C", "D"],
            ["A"], ["B"]]


def test_order_locks_bounded_under_churn(tiny, monkeypatch):
    monkeypatch.setattr(server, "_LRU_SIZE", 3)
    c = TestClient(create_app(tiny))
    pg = c.app.state.pg
    for sel in _distinct_selections():
        oid = c.post("/api/orders", json={"by": "contrast", "genomes": sel}).json()["order"]
        # force a cold packed build so per-order build locks would be allocated
        assert c.get("/api/matrix", params=dict(order=oid, og0=0, og1=7, g0=0, g1=4,
                                                w=7, h=4)).status_code == 200
    assert len(pg._dynamic) <= 3
    # dynamic churn must not grow the fixed lock map (built-ins only)
    assert set(pg._order_locks) <= {"tree", "frequency"}
    assert len({id(e["lock"]) for e in pg._dynamic.values()}) <= 3


def test_dynamic_order_eviction_recreation(tiny, monkeypatch):
    monkeypatch.setattr(server, "_LRU_SIZE", 1)
    c = TestClient(create_app(tiny))
    id1 = c.post("/api/orders", json={"by": "contrast", "genomes": ["A", "B"]}).json()["order"]
    c.post("/api/orders", json={"by": "contrast", "genomes": ["C", "D"]})  # evicts id1
    assert c.get("/api/orthogroups", params=dict(order=id1, og0=0, og1=7)).status_code == 404
    again = c.post("/api/orders", json={"by": "contrast", "genomes": ["A", "B"]}).json()["order"]
    assert again == id1
    assert c.get("/api/orthogroups", params=dict(order=again, og0=0, og1=7)).status_code == 200
    assert c.get("/api/matrix", params=dict(order=again, og0=0, og1=7, g0=0, g1=4,
                                            w=7, h=4)).status_code == 200


def test_repost_reuses_warm_cache(tiny, monkeypatch):
    monkeypatch.setattr(server, "_LRU_SIZE", 8)
    c = TestClient(create_app(tiny))
    pg = c.app.state.pg
    id1 = c.post("/api/orders", json={"by": "contrast", "genomes": ["A", "B"]}).json()["order"]
    assert c.get("/api/matrix", params=dict(order=id1, og0=0, og1=7, g0=0, g1=4,
                                            w=7, h=4)).status_code == 200
    warm = pg._dynamic[id1]["packed"]
    assert warm is not None
    # duplicate-equivalent repost reuses the entry instead of resetting packed=None
    id2 = c.post("/api/orders",
                 json={"by": "contrast", "genomes": ["B", "A", "A"]}).json()["order"]
    assert id2 == id1
    assert pg._dynamic[id1]["packed"] is warm


def test_orders_body_limit_chunked(client):
    # generator body -> httpx2 sends chunked, no Content-Length
    def chunks():
        yield b'{"by":"contrast","genomes":[' + b" " * server.MAX_BODY
        yield b"]}"

    r = client.post("/api/orders", content=chunks())
    assert r.status_code == 413


def test_orders_body_limit_crossed_on_later_chunk(client):
    # Raw ASGI: TestClient coalesces bodies, so feed chunks explicitly, no Content-Length.
    import asyncio

    half = server.MAX_BODY // 2
    chunks = [b" " * half, b" " * half, b" ", b"]}"]
    received, sent = [], []

    async def receive():
        body = chunks[len(received)]
        received.append(body)
        return {"type": "http.request", "body": body, "more_body": len(received) < len(chunks)}

    async def send(message):
        sent.append(message)

    scope = {"type": "http", "asgi": {"version": "3.0"}, "http_version": "1.1",
             "method": "POST", "scheme": "http", "path": "/api/orders", "raw_path": b"/api/orders",
             "query_string": b"", "root_path": "", "headers": [(b"content-type", b"application/json")],
             "client": ("127.0.0.1", 1), "server": ("127.0.0.1", 80)}
    asyncio.run(client.app(scope, receive, send))
    assert sent[0]["status"] == 413
    assert len(received) == 3  # rejected on the chunk crossing the limit; rest unread


def test_track_order_missing_last(client):
    ids = client.get("/api/orthogroups",
                     params=dict(order="track:score", og0=0, og1=7)).json()
    assert ids == ["OG_core", "OG_ac", "OG_cd", "OG_ab", "OG_a", "OG_c", "OG_empty"]
    ids_desc = client.get("/api/orthogroups",
                          params=dict(order="track:score:desc", og0=0, og1=7)).json()
    assert ids_desc == ["OG_a", "OG_ab", "OG_cd", "OG_ac", "OG_core", "OG_c", "OG_empty"]


def test_orders_body_validation(client):
    # non-list, or lists with non-string / unknown elements -> 400 (not 500)
    for genomes in ([[]], [{}], [1], [None], "A", 5, None, []):
        r = client.post("/api/orders", json={"by": "contrast", "genomes": genomes})
        assert r.status_code == 400, genomes
    assert client.post("/api/orders", json=[]).status_code == 400
    assert client.post("/api/orders", content=b"{not json",
                       headers={"content-type": "application/json"}).status_code == 400
    # more names than genomes -> 400
    assert client.post("/api/orders", json={"by": "contrast",
                                            "genomes": ["A", "B", "C", "D", "A"]}).status_code == 400
    # duplicates are canonicalized, so this is the same order as ["A", "B"]
    plain = client.post("/api/orders", json={"by": "contrast", "genomes": ["A", "B"]}).json()
    dupe = client.post("/api/orders", json={"by": "contrast",
                                            "genomes": ["A", "A", "B"]}).json()
    assert plain["order"] == dupe["order"]


def test_orders_body_limit(client):
    payload = b'{"by":"contrast","genomes":["A"]}' + b" " * (server.MAX_BODY + 1)
    r = client.post("/api/orders", content=payload,
                    headers={"content-type": "application/json"})
    assert r.status_code == 413


def test_bad_requests_400(client):
    p = dict(order="tree", og0=0, og1=7, g0=0, g1=4)
    assert client.get("/api/matrix", params={**p, "w": 8, "h": 4}).status_code == 400
    assert client.get("/api/matrix", params={**p, "w": 0, "h": 4}).status_code == 400
    assert client.get("/api/matrix", params={**p, "w": MAX_CELLS + 1, "h": 4}).status_code == 400
    assert client.get("/api/matrix", params=dict(order="tree", og0=5, og1=2, g0=0, g1=4,
                                                 w=2, h=4)).status_code == 400
    assert client.get("/api/matrix", params={**p, "w": 3, "h": 9}).status_code == 400
    assert client.get("/api/tracks", params=dict(order="tree", og0=0, og1=7, w=8,
                                                 names="score")).status_code == 400
    assert client.get("/api/orthogroups", params=dict(og0=0, og1=8)).status_code == 400
    assert client.post("/api/orders", json={"by": "nope"}).status_code == 400
    assert client.post(
        "/api/orders",
        json={"by": "contrast", "genomes": ["A", "B", "C", "D"]},
    ).status_code == 400
    assert client.post("/api/orders", json={"by": "contrast",
                                            "genomes": ["A", "ZZZ"]}).status_code == 400


def test_matrix_cell_cap(tiny, monkeypatch):
    monkeypatch.setattr(server, "MAX_CELLS_TOTAL", 6)
    c = TestClient(create_app(tiny))
    p = dict(order="tree", og0=0, og1=7, g0=0, g1=4)
    assert c.get("/api/matrix", params={**p, "w": 3, "h": 2}).status_code == 200
    assert c.get("/api/matrix", params={**p, "w": 4, "h": 2}).status_code == 400


def test_unknown_404(client):
    assert client.get("/api/orthogroup/NOPE").status_code == 404
    assert client.get("/api/matrix", params=dict(order="nope", og0=0, og1=7, g0=0, g1=4,
                                                 w=3, h=4)).status_code == 404
    assert client.get("/api/tracks", params=dict(order="tree", og0=0, og1=7, w=3,
                                                 names="nope")).status_code == 404
