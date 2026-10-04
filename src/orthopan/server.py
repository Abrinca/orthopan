"""`orthopan serve`: FastAPI backend for an imported pangenome directory.

Serves the JSON metadata endpoints, the bit-packed matrix and track binaries, and
(optionally) the built frontend under ``/``. The contract is ``docs/api.md``.

The presence matrix is loaded once at startup and reordered to tree row order.
For each order a column-permuted *packed* copy is cached lazily (small: one bit
per cell). Per request the requested window is unpacked on the fly and reduced
to output bins with a popcount prefix sum, never a Python loop over cells.
"""

from __future__ import annotations

import hashlib
import json
import threading
from collections import OrderedDict
from pathlib import Path

import numpy as np
import polars as pl
import uvicorn
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import Response
from fastapi.staticfiles import StaticFiles

from orthopan.importer import _parse_newick

DEFAULT_TRACKS = ("copy_number_mean", "copy_number_max")
MAX_CELLS = 4096
MAX_CELLS_TOTAL = 4_000_000  # w*h cap per matrix request
MAX_LIST = 5000
MAX_BODY = 64 * 1024  # request-body cap for the POST endpoint
_LRU_SIZE = 24
_ROW_CHUNK = 512
_DYNAMIC_PREFIX = "contrast:"  # order ids minted by POST /api/orders

# Top-r-bit masks for masking a partial byte; MASKS[0] == 0.
MASKS = np.array([((0xFF << (8 - r)) & 0xFF) for r in range(8)], np.uint8)


def _edges(n: int, k: int) -> np.ndarray:
    """floor boundaries of k bins over [0, n): floor(i*n/k)."""
    return (np.arange(k + 1, dtype=np.int64) * n) // k


class Pangenome:
    """Everything the endpoints need, loaded once at app creation."""

    def __init__(self, data_dir: str | Path):
        self.dir = Path(data_dir)
        manifest_path = self.dir / "manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"not an imported pangenome directory: {self.dir}")
        self.manifest = json.loads(manifest_path.read_text())

        self.genomes = pl.read_parquet(self.dir / "genomes.parquet").sort("genome_idx")
        self.orthogroups = pl.read_parquet(self.dir / "orthogroups.parquet").sort("og_idx")
        self.n_genomes = self.genomes.height
        self.n_og = self.orthogroups.height

        self.genome_names = self.genomes["name"].to_list()
        self._name_to_idx = {name: i for i, name in enumerate(self.genome_names)}
        row_order = self.genomes["row_order"].to_numpy()
        inv = np.empty(self.n_genomes, dtype=np.int64)
        inv[row_order] = np.arange(self.n_genomes, dtype=np.int64)
        self.row_order = row_order  # genome_idx -> tree row

        presence = np.load(self.dir / "presence.npy")
        self.presence_tree = presence[inv]  # tree row order
        self.nbytes = presence.shape[1]
        self.padded_nbytes = self.nbytes + 1

        # Order permutations: position -> og_idx.
        self.perm_tree = np.argsort(self.orthogroups["order_tree"].to_numpy(), kind="stable")
        self.perm_frequency = np.argsort(
            self.orthogroups["order_frequency"].to_numpy(), kind="stable"
        )

        self.track_names = list(DEFAULT_TRACKS)
        for name in self.manifest.get("track_columns", []):
            if name not in self.track_names:
                self.track_names.append(name)
        self._packed_cache: OrderedDict[str, np.ndarray] = OrderedDict()
        # Dynamic order entries hold {perm, packed}; one LRU entry = one order, so
        # eviction drops the permutation and its packed matrix together (no 200/404 skew).
        self._dynamic: OrderedDict[str, dict] = OrderedDict()
        self._order_locks: dict[str, threading.Lock] = {}
        # ponytail: one lock guards the cache bookkeeping; per-order locks guard builds.
        self._lock = threading.Lock()

        self.tree = self._build_tree()
        self._packed("tree")  # prewarm the default view

    # --- tree ---------------------------------------------------------------
    def _build_tree(self):
        text = (self.dir / "tree.nwk").read_text()
        root = _parse_newick(text)
        parent, length, label, row = [], [], [], []
        stack = [(root, -1, 0.0)]
        leaf = 0
        while stack:
            node, par, branch = stack.pop()
            i = len(parent)
            parent.append(par)
            length.append(branch)
            label.append(node.label or "")
            if node.children:
                row.append(-1)
                for child in reversed(node.children):
                    stack.append((child, i, float(child.branch) if child.branch is not None else 0.0))
            else:
                row.append(leaf)
                leaf += 1
        return {"parent": parent, "length": length, "label": label, "row": row}

    # --- orders -------------------------------------------------------------
    def order_perm(self, order: str) -> np.ndarray:
        if order == "tree":
            return self.perm_tree
        if order == "frequency":
            return self.perm_frequency
        if order.startswith("track:"):
            return self._track_perm(order[6:])
        entry = self._dynamic_get(order)
        if entry is None:
            raise HTTPException(status_code=404, detail=f"unknown order: {order!r}")
        return entry["perm"]

    def _dynamic_get(self, order: str):
        """Return the live dynamic entry (touch LRU) or None if unknown/evicted."""
        if not order.startswith(_DYNAMIC_PREFIX):
            return None
        with self._lock:
            entry = self._dynamic.get(order)
            if entry is not None:
                self._dynamic.move_to_end(order)
            return entry

    def _order_lock(self, order: str) -> threading.Lock:
        """Build lock for non-dynamic orders (bounded: built-ins + known tracks)."""
        with self._lock:
            lock = self._order_locks.get(order)
            if lock is None:
                lock = self._order_locks[order] = threading.Lock()
            return lock

    def _track_perm(self, spec: str) -> np.ndarray:
        desc = spec.endswith(":desc")
        name = spec[:-5] if desc else spec
        if name not in self.track_names:
            raise HTTPException(status_code=404, detail=f"unknown track: {name!r}")
        values = self.orthogroups[name].to_numpy()
        missing = np.isnan(values)
        present = np.flatnonzero(~missing)
        ranked = np.argsort(values[present], kind="stable")
        if desc:
            ranked = ranked[::-1]
        return np.concatenate([present[ranked], np.flatnonzero(missing)]).astype(np.int64)

    def _build_packed(self, perm: np.ndarray) -> np.ndarray:
        out = np.zeros((self.n_genomes, self.padded_nbytes), np.uint8)
        for r0 in range(0, self.n_genomes, _ROW_CHUNK):
            r1 = min(r0 + _ROW_CHUNK, self.n_genomes)
            dense = np.unpackbits(self.presence_tree[r0:r1], axis=1, bitorder="big")[:, : self.n_og]
            out[r0:r1, : self.nbytes] = np.packbits(dense[:, perm], axis=1, bitorder="big")
        return out

    def _packed(self, order: str) -> np.ndarray:
        if order.startswith(_DYNAMIC_PREFIX):
            return self._packed_dynamic(order)
        with self._lock:
            cached = self._packed_cache.get(order)
            if cached is not None:
                self._packed_cache.move_to_end(order)
                return cached
        perm = self.order_perm(order)
        with self._order_lock(order):
            with self._lock:
                cached = self._packed_cache.get(order)
                if cached is not None:
                    self._packed_cache.move_to_end(order)
                    return cached
            out = self._build_packed(perm)
            with self._lock:
                self._packed_cache[order] = out
                while len(self._packed_cache) > _LRU_SIZE:
                    self._packed_cache.popitem(last=False)
        return out

    def _packed_dynamic(self, order: str) -> np.ndarray:
        entry = self._dynamic_get(order)
        if entry is None:
            raise HTTPException(status_code=404, detail=f"unknown order: {order!r}")
        if entry["packed"] is not None:
            return entry["packed"]
        # Hold the entry's own lock: if it is evicted and recreated mid-build, the
        # waiter keeps using this reference instead of racing a fresh entry.
        with entry["lock"]:
            if entry["packed"] is not None:
                return entry["packed"]
            out = self._build_packed(entry["perm"])
            with self._lock:
                entry["packed"] = out
                live = self._dynamic.get(order)
                if live is not None and live is not entry:
                    live["packed"] = out
            return out

    # --- matrix -------------------------------------------------------------
    def matrix(self, order, og0, og1, g0, g1, w, h) -> np.ndarray:
        packed = self._packed(order)
        col_abs = og0 + _edges(og1 - og0, w)
        row_local = _edges(g1 - g0, h)

        b = col_abs >> 3
        r = col_abs & 7
        sub = packed[g0:g1]
        bc = np.bitwise_count(sub)
        cs = np.zeros((sub.shape[0], self.padded_nbytes + 1), np.int32)
        np.cumsum(bc, axis=1, dtype=np.int32, out=cs[:, 1:])
        counts = np.bitwise_count(sub[:, b] & MASKS[r])  # (H, w+1)
        prefix = cs[:, b] + counts  # (H, w+1)
        colsum = prefix[:, 1:] - prefix[:, :-1]  # (H, w)

        rowsum = np.add.reduceat(colsum, row_local[:-1], axis=0, dtype=np.int32)  # (h, w)
        row_sizes = np.diff(row_local)[:, None]
        col_sizes = np.diff(_edges(og1 - og0, w))[None, :]
        frac = rowsum / (row_sizes * col_sizes)
        return np.rint(frac * 255.0).clip(0, 255).astype("<u1")

    # --- tracks -------------------------------------------------------------
    def track_means(self, order, og0, og1, w, names) -> np.ndarray:
        perm = self.order_perm(order)
        starts = _edges(og1 - og0, w)[:-1]
        out = np.empty((len(names), w), np.float32)
        for i, name in enumerate(names):
            ordered = self.orthogroups[name].to_numpy()[perm][og0:og1]
            valid = ~np.isnan(ordered)
            sums = np.add.reduceat(np.where(valid, ordered, 0.0), starts, dtype=np.float64)
            counts = np.add.reduceat(valid.astype(np.float64), starts, dtype=np.float64)
            out[i] = np.divide(sums, counts, out=np.full(w, np.nan), where=counts > 0)
        return out

    # --- contrast -----------------------------------------------------------
    def contrast_perm(self, names) -> np.ndarray:
        rows = []
        for name in names:
            idx = self._name_to_idx.get(name)
            if idx is None:
                raise HTTPException(status_code=400, detail=f"unknown genome: {name!r}")
            rows.append(self.row_order[idx])
        if not rows or len(set(rows)) == self.n_genomes:
            raise HTTPException(status_code=400, detail="contrast needs 1..n-1 genomes")
        sel = np.zeros(self.n_genomes, dtype=bool)
        sel[rows] = True
        dense = np.unpackbits(self.presence_tree[sel], axis=1, bitorder="big")[:, : self.n_og]
        inside = dense.sum(axis=0, dtype=np.float64)
        total = self.orthogroups["n_genomes"].to_numpy().astype(np.float64)
        outside = total - inside
        nsel = sel.sum()
        score = inside / nsel - outside / (self.n_genomes - nsel)
        return np.argsort(-score, kind="stable").astype(np.int64)

    # --- metadata -----------------------------------------------------------
    def meta(self) -> dict:
        tracks = []
        for name in self.track_names:
            values = self.orthogroups[name].to_numpy()
            missing = int(np.isnan(values).sum())
            valid = values[~np.isnan(values)]
            tracks.append(
                {
                    "name": name,
                    "min": float(valid.min()) if valid.size else None,
                    "max": float(valid.max()) if valid.size else None,
                    "missing": missing,
                }
            )
        return {
            "n_genomes": self.n_genomes,
            "n_orthogroups": self.n_og,
            "tree_source": self.manifest.get("tree_source"),
            "orders": ["tree", "frequency"],
            "tracks": tracks,
        }


def _window(n_og, n_genomes, og0, og1, g0, g1, w, h):
    og1 = n_og if og1 is None else og1
    g1 = n_genomes if g1 is None else g1
    if not (0 <= og0 < og1 <= n_og) or not (0 <= g0 < g1 <= n_genomes):
        raise HTTPException(status_code=400, detail="window out of range")
    if not (1 <= w <= MAX_CELLS and 1 <= h <= MAX_CELLS):
        raise HTTPException(status_code=400, detail=f"w and h must be in 1..{MAX_CELLS}")
    if w > og1 - og0 or h > g1 - g0:
        raise HTTPException(status_code=400, detail="w/h exceed window size")
    if w * h > MAX_CELLS_TOTAL:
        raise HTTPException(
            status_code=400, detail=f"w*h exceeds {MAX_CELLS_TOTAL} output cells"
        )
    return og0, og1, g0, g1


def create_app(data_dir, static_dir=None) -> FastAPI:
    pg = Pangenome(data_dir)
    app = FastAPI(title="OrthoPan", docs_url=None, redoc_url=None)
    app.state.pg = pg

    @app.get("/api/meta")
    def meta():
        return pg.meta()

    @app.get("/api/genomes")
    def genomes():
        rows = []
        for row in pg.genomes.iter_rows(named=True):
            rows.append(
                {
                    "row": int(row["row_order"]),
                    "genome_idx": int(row["genome_idx"]),
                    "name": row["name"],
                    "n_orthogroups": int(row["n_orthogroups"]),
                    "n_genes": int(row["n_genes"]),
                }
            )
        rows.sort(key=lambda r: r["row"])
        return rows

    @app.get("/api/tree")
    def tree():
        return pg.tree

    @app.get("/api/matrix")
    def matrix(
        order: str = "tree",
        og0: int = 0,
        og1: int | None = None,
        g0: int = 0,
        g1: int | None = None,
        w: int = Query(...),
        h: int = Query(...),
    ):
        og0, og1, g0, g1 = _window(pg.n_og, pg.n_genomes, og0, og1, g0, g1, w, h)
        out = pg.matrix(order, og0, og1, g0, g1, w, h)
        return Response(
            content=out.tobytes(),
            media_type="application/octet-stream",
            headers={"X-Width": str(w), "X-Height": str(h)},
        )

    @app.get("/api/tracks")
    def tracks(
        order: str = "tree",
        og0: int = 0,
        og1: int | None = None,
        w: int = Query(...),
        names: str | None = None,
    ):
        og1 = pg.n_og if og1 is None else og1
        if not (0 <= og0 < og1 <= pg.n_og):
            raise HTTPException(status_code=400, detail="window out of range")
        if not (1 <= w <= MAX_CELLS):
            raise HTTPException(status_code=400, detail=f"w must be in 1..{MAX_CELLS}")
        if w > og1 - og0:
            raise HTTPException(status_code=400, detail="w exceeds window size")
        selected = names.split(",") if names else list(pg.track_names)
        for name in selected:
            if name not in pg.track_names:
                raise HTTPException(status_code=404, detail=f"unknown track: {name!r}")
        out = pg.track_means(order, og0, og1, w, selected)
        return Response(
            content=out.astype("<f4").tobytes(),
            media_type="application/octet-stream",
            headers={"X-Tracks": ",".join(selected)},
        )

    @app.get("/api/orthogroups")
    def orthogroups(order: str = "tree", og0: int = 0, og1: int | None = None):
        og1 = pg.n_og if og1 is None else og1
        if not (0 <= og0 < og1 <= pg.n_og):
            raise HTTPException(status_code=400, detail="window out of range")
        if og1 - og0 > MAX_LIST:
            raise HTTPException(status_code=400, detail=f"at most {MAX_LIST} per request")
        perm = pg.order_perm(order)
        ids = pg.orthogroups["orthogroup"].to_numpy()[perm[og0:og1]]
        return ids.tolist()

    @app.get("/api/orthogroup/{og}")
    def orthogroup(og: str):
        frame = pg.orthogroups.filter(pl.col("orthogroup") == og)
        if not frame.height:
            raise HTTPException(status_code=404, detail=f"unknown orthogroup: {og!r}")
        row = frame.row(0, named=True)
        og_idx = int(row["og_idx"])
        tracks = {}
        for name in pg.track_names:
            value = row[name]
            tracks[name] = None if value is None or np.isnan(value) else float(value)
        membership = pg.dir / "membership.parquet"
        members = []
        if membership.is_file():
            mem = (
                pl.scan_parquet(membership)
                .select(["og_idx", "genome_idx", "genes"])
                .filter(pl.col("og_idx") == og_idx)
                .collect()
                .sort("genome_idx")
            )
            for m in mem.iter_rows(named=True):
                members.append(
                    {
                        "genome": pg.genome_names[int(m["genome_idx"])],
                        "genes": list(m["genes"]),
                    }
                )
        return {
            "orthogroup": row["orthogroup"],
            "og_idx": og_idx,
            "n_genomes": int(row["n_genomes"]),
            "n_genes": int(row["n_genes"]),
            "mrca_label": row["mrca_label"],
            "purity": None if row["purity"] is None else float(row["purity"]),
            "tracks": tracks,
            "members": members,
        }

    @app.post("/api/orders")
    async def create_order(request: Request):
        length = request.headers.get("content-length")
        if length is not None and length.isdigit() and int(length) > MAX_BODY:
            raise HTTPException(status_code=413, detail="request body too large")
        raw = bytearray()
        async for chunk in request.stream():
            raw += chunk
            if len(raw) > MAX_BODY:
                raise HTTPException(status_code=413, detail="request body too large")
        try:
            body = json.loads(raw)
        except (ValueError, UnicodeDecodeError):
            raise HTTPException(status_code=400, detail="invalid JSON body")
        if not isinstance(body, dict) or body.get("by") != "contrast":
            raise HTTPException(status_code=400, detail="only 'contrast' orders are supported")
        names = body.get("genomes")
        if not isinstance(names, list) or not names or not all(isinstance(n, str) for n in names):
            raise HTTPException(status_code=400, detail="'genomes' must be a non-empty list of strings")
        if len(names) > pg.n_genomes:
            raise HTTPException(status_code=400, detail="more genomes than exist")
        unique = list(dict.fromkeys(names))  # canonicalize duplicates, keep order
        digest = hashlib.sha1(json.dumps(sorted(set(names))).encode()).hexdigest()[:12]
        order_id = _DYNAMIC_PREFIX + digest
        if pg._dynamic_get(order_id) is not None:
            return {"order": order_id}  # reuse the existing (possibly warm) entry
        perm = pg.contrast_perm(unique)
        with pg._lock:
            if order_id in pg._dynamic:
                pg._dynamic.move_to_end(order_id)
                return {"order": order_id}
            pg._dynamic[order_id] = {"perm": perm, "packed": None, "lock": threading.Lock()}
            pg._dynamic.move_to_end(order_id)
            while len(pg._dynamic) > _LRU_SIZE:
                pg._dynamic.popitem(last=False)
        return {"order": order_id}

    if static_dir is None:
        candidate = Path(__file__).resolve().parents[2] / "web" / "dist"
        static_dir = candidate if candidate.is_dir() else None
    if static_dir is not None:
        app.mount("/", StaticFiles(directory=str(static_dir), html=True), name="static")

    return app


def run(args) -> None:
    app = create_app(args.data_dir, args.static)
    uvicorn.run(app, host=args.host, port=args.port)


def add_parser(subparsers):
    p = subparsers.add_parser("serve", help="serve an imported pangenome directory")
    p.add_argument("data_dir", metavar="DIR", help="directory written by `orthopan import`")
    p.add_argument("--host", default="127.0.0.1", help="bind host (default 127.0.0.1)")
    p.add_argument("--port", type=int, default=8000, help="bind port (default 8000)")
    p.add_argument("--static", metavar="PATH", help="built frontend directory to serve at /")
    p.set_defaults(func=run)
