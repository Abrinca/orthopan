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
import os
import shutil
import subprocess
import tempfile
import threading
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote

import numpy as np
import polars as pl
import uvicorn
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import Response
from fastapi.staticfiles import StaticFiles

from orthopan import sequences as seqmod
from orthopan.importer import _parse_newick
from orthopan.io import read_matrix
from orthopan.tracks.msa import _expected_headers, _match_headers, read_fasta
from orthopan.tree import _quote

DEFAULT_TRACKS = ("copy_number_mean", "copy_number_max")
MAX_CELLS = 4096
MAX_CELLS_TOTAL = 4_000_000  # w*h cap per matrix request
MAX_LIST = 5000
MAX_BODY = 64 * 1024  # request-body cap for the POST endpoint
_LRU_SIZE = 8  # packed orders kept per cache (each is n_genomes x nbytes)
_ROW_CHUNK = 48  # rows per _build_packed task (small chunks load-balance well)
# order ids minted by POST /api/orders
_DYNAMIC_PREFIX = ("contrast:", "search:", "pattern:")
_FAMSA_CONCURRENCY = 2
_FAMSA_MAX_SEQS = 2000
_FAMSA_TIMEOUT = 120
# Shared worker pool for the block-parallel matrix/build paths (numpy releases
# the GIL, so real threads scale on these byte ops).
_MATRIX_WORKERS = min(os.cpu_count() or 1, 16)
_MATRIX_POOL: ThreadPoolExecutor | None = None
_MATRIX_POOL_LOCK = threading.Lock()


def _matrix_pool() -> ThreadPoolExecutor:
    global _MATRIX_POOL
    if _MATRIX_POOL is None:
        with _MATRIX_POOL_LOCK:
            if _MATRIX_POOL is None:
                _MATRIX_POOL = ThreadPoolExecutor(max_workers=_MATRIX_WORKERS)
    return _MATRIX_POOL


def _famsa_ok(path: str) -> bool:
    """True when the famsa binary resolves and aligns two tiny sequences."""
    exe = shutil.which(path)
    if exe is None:
        return False
    try:
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "in.fa", Path(tmp) / "out.fa"
            src.write_text(">a\nAAAA\n>b\nAAAT\n")
            subprocess.run([exe, str(src), str(dst), "-t", "1"],
                           check=True, capture_output=True, timeout=30)
            pairs = _parse_fasta_text(dst.read_text())
            return (len(pairs) == 2 and all(seq for _, seq in pairs)
                    and len({len(seq) for _, seq in pairs}) == 1)
    except (OSError, subprocess.SubprocessError, ValueError):
        return False


def _parse_fasta_text(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into ``[(header, sequence), ...]``."""
    out, name, buf = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line[0] == ">":
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def _to_newick(root, name_of) -> str:
    """Serialize a parsed Newick node tree, naming leaves via ``name_of``."""
    out: list[str] = []
    stack = [("node", root)]
    while stack:
        kind, obj = stack.pop()
        if kind == "str":
            out.append(obj)
            continue
        node = obj
        if not node.children:
            out.append(_quote(name_of(node)))
            if node.branch is not None:
                out.append(f":{node.branch:g}")
            continue
        out.append("(")
        suffix = ")"
        if node.label:
            suffix += _quote(node.label)
        if node.branch is not None:
            suffix += f":{node.branch:g}"
        stack.append(("str", suffix))
        items = []
        for i, child in enumerate(node.children):
            if i:
                items.append(("str", ","))
            items.append(("node", child))
        for item in reversed(items):
            stack.append(item)
    return "".join(out) + ";"


def _content_disposition(filename: str) -> str:
    """RFC 6266/5987 attachment header with an ASCII fallback and UTF-8 name."""
    if any(ord(c) < 32 or ord(c) == 127 for c in filename):
        raise HTTPException(status_code=400, detail="invalid filename")
    fallback = filename.encode("ascii", "replace").decode("ascii")
    fallback = fallback.replace("\\", "_").replace('"', "_")
    return (f'attachment; filename="{fallback}"; '
            f"filename*=UTF-8''{quote(filename, safe='')}")


# Top-r-bit masks for masking a partial byte; MASKS[0] == 0.
MASKS = np.array([((0xFF << (8 - r)) & 0xFF) for r in range(8)], np.uint8)


def _edges(n: int, k: int) -> np.ndarray:
    """floor boundaries of k bins over [0, n): floor(i*n/k)."""
    return (np.arange(k + 1, dtype=np.int64) * n) // k


class Pangenome:
    """Everything the endpoints need, loaded once at app creation."""

    def __init__(self, data_dir: str | Path, famsa: str = "famsa"):
        self.dir = Path(data_dir)
        manifest_path = self.dir / "manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"not an imported pangenome directory: {self.dir}")
        self.manifest = json.loads(manifest_path.read_text())
        version = self.manifest.get("format_version", 1)
        if version < 2:
            raise ValueError(
                f"dataset format {version} is outdated; re-run orthopan import"
            )

        self.genomes = pl.read_parquet(self.dir / "genomes.parquet").sort("genome_idx")
        self.orthogroups = pl.read_parquet(self.dir / "orthogroups.parquet").sort("og_idx")
        self.n_genomes = self.genomes.height
        self.n_og = self.orthogroups.height

        self.genome_names = self.genomes["name"].to_list()
        self._name_to_idx = {name: i for i, name in enumerate(self.genome_names)}
        self._og_to_idx = {
            og: int(i) for og, i in zip(
                self.orthogroups["orthogroup"].to_list(), self.orthogroups["og_idx"].to_list())
        }
        row_order = self.genomes["row_order"].to_numpy()
        inv = np.empty(self.n_genomes, dtype=np.int64)
        inv[row_order] = np.arange(self.n_genomes, dtype=np.int64)
        self.row_order = row_order  # genome_idx -> tree row

        self.pattern_ids = self.orthogroups["pattern_id"].to_numpy()
        self.n_same_pattern = self.orthogroups["n_same_pattern"].to_numpy()

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

        # --- v1 additions: annotations, genome metadata, multicopy ----------
        self.annotation_specs = self.manifest.get("annotations", [])
        self.annotation_names = [s["name"] for s in self.annotation_specs]
        self._ann_kind = {s["name"]: s["kind"] for s in self.annotation_specs}
        self._ann_categories = {
            s["name"]: list(s.get("categories", [])) for s in self.annotation_specs
        }
        self.annotations = None
        ann_path = self.dir / "annotations.parquet"
        if self.annotation_names and ann_path.is_file():
            self.annotations = pl.read_parquet(ann_path).sort("og_idx")
        self._ann_codes: dict[str, np.ndarray] = {}
        if self.annotations is not None:
            for spec in self.annotation_specs:
                if spec["kind"] != "categorical":
                    continue
                index = {c: i + 1 for i, c in enumerate(self._ann_categories[spec["name"]])}
                column = self.annotations[spec["name"]].to_list()
                self._ann_codes[spec["name"]] = np.array(
                    [index.get(v, 0) if v is not None else 0 for v in column], np.uint16
                )
        # (column name, frame) pairs searched by /api/search.
        self._search_fields = [("orthogroup", self.orthogroups)]
        if self.annotations is not None:
            for spec in self.annotation_specs:
                self._search_fields.append((spec["name"], self.annotations))

        self.gm_specs = self.manifest.get("genome_metadata", [])
        self.gm_names = [s["name"] for s in self.gm_specs]
        gm_path = self.dir / "genome_metadata.parquet"
        self.genome_metadata = None
        self._genome_meta = [dict() for _ in range(self.n_genomes)]
        if self.gm_names and gm_path.is_file():
            self.genome_metadata = pl.read_parquet(gm_path).sort("genome_idx")
            for row in self.genome_metadata.iter_rows(named=True):
                self._genome_meta[int(row["genome_idx"])] = {
                    name: row[name] for name in self.gm_names
                }

        mc_path = self.dir / "multicopy.npz"
        if mc_path.is_file():
            with np.load(mc_path) as z:
                self.mc_row = z["row"].astype(np.int64)
                self.mc_og = z["og"].astype(np.int64)
                self.mc_extra = z["extra"].astype(np.int64)
        else:
            self.mc_row = np.empty(0, np.int64)
            self.mc_og = np.empty(0, np.int64)
            self.mc_extra = np.empty(0, np.int64)

        self.detail = self.manifest.get("detail", {}) or {}

        # --- v2: single-orthogroup detail sources ---------------------------
        self.fasta_files = [str(p) for p in self.manifest.get("fasta", [])]
        self._seq_lookup: dict[tuple[str, str], tuple[int, int, int]] = {}
        seq_path = self.dir / "sequences.parquet"
        if self.fasta_files and seq_path.is_file():
            for r in pl.read_parquet(seq_path).iter_rows(named=True):
                self._seq_lookup[(r["gene"], r["genome"])] = (
                    int(r["file"]), int(r["offset"]), int(r["length"]))
        self.msa_dir = Path(self.detail["msa"]) if self.detail.get("msa") else None
        self.seq_dir = Path(self.detail["sequences"]) if self.detail.get("sequences") else None
        self.gene_tree_file = (
            Path(self.detail["gene_trees"]) if self.detail.get("gene_trees") else None)
        if self.msa_dir is not None:
            self.results_dir = self.msa_dir.parent
        elif self.seq_dir is not None:
            self.results_dir = self.seq_dir.parent
        else:
            self.results_dir = None
        self._has_sequences = bool(self.seq_dir) or bool(self._seq_lookup)
        self.famsa = famsa
        self.famsa_ok = _famsa_ok(famsa)
        self.cache_dir = self.dir / "cache" / "msa"
        self._famsa_sem = threading.BoundedSemaphore(_FAMSA_CONCURRENCY)
        self._famsa_locks: dict[str, threading.Lock] = {}
        self._gt_index = None
        self._gt_lock = threading.Lock()
        self._dup_index = None
        self._dup_lock = threading.Lock()

        self.tree = self._build_tree()
        self._genomes_json = self._build_genomes_json()
        self._packed("tree")  # prewarm the default view (fast, ~row-chunked)
        # Prebuild the frequency order in the background so the first request
        # for it is warm without delaying startup.
        threading.Thread(
            target=self._prewarm, args=("frequency",), daemon=True
        ).start()

    def _prewarm(self, order: str) -> None:
        try:
            self._packed(order)
        except Exception:  # pragma: no cover - background best effort
            pass

    def _build_genomes_json(self) -> bytes:
        rows = []
        for row in self.genomes.iter_rows(named=True):
            genome_idx = int(row["genome_idx"])
            rows.append({
                "row": int(row["row_order"]),
                "genome_idx": genome_idx,
                "name": row["name"],
                "n_orthogroups": int(row["n_orthogroups"]),
                "n_genes": int(row["n_genes"]),
                "metadata": self._genome_meta[genome_idx],
            })
        rows.sort(key=lambda r: r["row"])
        # compact like Starlette's JSONResponse (no ", "/": " padding)
        return json.dumps(rows, ensure_ascii=False, separators=(",", ":")).encode("utf-8")

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

    # --- v2: single-orthogroup detail ---------------------------------------
    def og_members(self, og_idx: int):
        """Return ``[(genome_idx, genome, display_genes, raw_genes), ...]``."""
        mem = (
            pl.scan_parquet(self.dir / "membership.parquet")
            .select(["og_idx", "genome_idx", "genes", "genes_raw"])
            .filter(pl.col("og_idx") == og_idx)
            .collect()
            .sort("genome_idx")
        )
        names = self.genome_names
        genome_idx = mem["genome_idx"].to_numpy()
        genes = mem["genes"].to_list()
        raws = mem["genes_raw"].to_list()
        return [
            (int(g), names[int(g)], gg, rr)
            for g, gg, rr in zip(genome_idx, genes, raws)
        ]

    def _member_context(self, og_idx: int):
        members = self.og_members(og_idx)
        member_map: dict[str, list[str]] = {}
        raw_to_display: dict[str, str] = {}
        for _, gname, genes, raws in members:
            member_map[gname] = raws
            for d, r in zip(genes, raws):
                raw_to_display[r] = d
        return members, member_map, raw_to_display, _expected_headers(member_map)

    def _row_of(self, genome) -> int:
        idx = self._name_to_idx.get(genome)
        return int(self.row_order[idx]) if idx is not None else -1

    def _trimmed_value(self, og_idx: int):
        if "trimmed" not in self.orthogroups.columns:
            return None
        value = self.orthogroups["trimmed"].to_numpy()[og_idx]
        return None if value is None or np.isnan(value) else bool(value)

    def og_sequences(self, og_idx: int):
        """Unaligned sequences for an OG as ``[(raw, display, genome, seq), ...]``."""
        og = self.orthogroups["orthogroup"].to_numpy()[og_idx]
        members, _, raw_to_display, _ = self._member_context(og_idx)
        raw_genome = {
            r: gname for _, gname, _, raws in members for r in raws
        }
        rows = []
        path = (self.seq_dir / f"{og}.fa") if self.seq_dir else None
        if path is not None and path.is_file():
            for header, seq in read_fasta(path):
                rows.append((header, raw_to_display.get(header, header),
                             raw_genome.get(header), seq))
        else:
            for _, gname, genes, raws in members:
                for display, raw in zip(genes, raws):
                    record = self._seq_lookup.get((raw, gname))
                    if record is None:
                        continue
                    _, seq = seqmod.read_record(self.fasta_files, record)
                    rows.append((raw, display, gname, seq))
        rows.sort(key=lambda t: (self._row_of(t[2]) if t[2] else 1 << 30, t[0]))
        return rows

    def _msa_text(self, og: str, og_idx: int):
        """Return ``(fasta_text, source)``, generating with famsa if needed."""
        if self.msa_dir is not None:
            path = self.msa_dir / f"{og}.fa"
            if path.is_file():
                return path.read_text(), "orthofinder"
        cache = self.cache_dir / f"{og}.fa"
        if cache.is_file():
            return cache.read_text(), "famsa"
        seqs = self.og_sequences(og_idx)
        if not seqs:
            raise HTTPException(status_code=404, detail="no MSA and no sequences available")
        if len(seqs) > _FAMSA_MAX_SEQS:
            raise HTTPException(
                status_code=404,
                detail=f"no MSA; {len(seqs)} sequences exceed the {_FAMSA_MAX_SEQS} limit",
            )
        if not self.famsa_ok:
            raise HTTPException(status_code=404, detail="no MSA and famsa is not available")
        if not self._famsa_sem.acquire(blocking=False):
            raise HTTPException(status_code=503, detail="alignment queue is busy")
        try:
            with self._famsa_lock(og):
                if cache.is_file():  # another request finished while we waited
                    return cache.read_text(), "famsa"
                raw_text = self._run_famsa(og, seqs)
                text = self._canonicalize_famsa(raw_text, seqs)
                self._write_cache(og, text)
                return text, "famsa"
        finally:
            self._famsa_sem.release()

    def _famsa_lock(self, og: str) -> threading.Lock:
        with self._lock:
            lock = self._famsa_locks.get(og)
            if lock is None:
                lock = self._famsa_locks[og] = threading.Lock()
            return lock

    @staticmethod
    def _famsa_headers(n: int) -> list[str]:
        return [f"m{i}" for i in range(n)]

    def _run_famsa(self, og: str, seqs) -> str:
        # Monkeypatched in tests; returns raw famsa output (indexed m<i> headers).
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        fd_in, tmp_in = tempfile.mkstemp(suffix=".fa", dir=self.cache_dir)
        fd_out, tmp_out = tempfile.mkstemp(suffix=".fa", dir=self.cache_dir)
        os.close(fd_in)
        os.close(fd_out)
        try:
            with open(tmp_in, "w") as fh:
                for header, (_, _, _, seq) in zip(self._famsa_headers(len(seqs)), seqs):
                    fh.write(f">{header}\n{seq}\n")
            try:
                subprocess.run(
                    [self.famsa, tmp_in, tmp_out, "-t", "1"],
                    check=True, capture_output=True, timeout=_FAMSA_TIMEOUT,
                )
            except subprocess.TimeoutExpired:
                raise HTTPException(status_code=500, detail="famsa timed out") from None
            except subprocess.CalledProcessError as exc:
                raise HTTPException(
                    status_code=500,
                    detail=f"famsa failed: {exc.stderr.decode(errors='replace')[:200]}",
                ) from None
            return Path(tmp_out).read_text()
        finally:
            for path in (tmp_in, tmp_out):
                try:
                    os.unlink(path)
                except OSError:
                    pass

    def _canonicalize_famsa(self, text: str, seqs) -> str:
        """Map famsa's indexed headers back to membership IDs; validate output.

        Rejects incomplete/wrong-identity output, empty records and unequal
        alignment lengths, so bad output is never cached.
        """
        pairs = _parse_fasta_text(text)
        mapping = {f"m{i}": (raw, genome)
                   for i, (raw, _, genome, _) in enumerate(seqs)}
        if sorted(h for h, _ in pairs) != sorted(mapping):  # exact permutation
            raise HTTPException(
                status_code=500, detail="famsa produced an incomplete alignment")
        lengths, out = set(), []
        for header, seq in pairs:
            if header not in mapping or not seq:
                raise HTTPException(
                    status_code=500, detail="famsa produced an invalid alignment")
            lengths.add(len(seq))
            raw, genome = mapping[header]
            canonical = f"{genome.replace('.', '_')}_{raw}" if genome else raw
            out.append(f">{canonical}\n{seq}\n")
        if len(lengths) != 1:
            raise HTTPException(
                status_code=500, detail="famsa produced an invalid alignment")
        return "".join(out)

    def _write_cache(self, og: str, text: str) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(suffix=".fa.tmp", dir=self.cache_dir)
        os.close(fd)
        try:
            Path(tmp).write_text(text)
            os.replace(tmp, self.cache_dir / f"{og}.fa")  # unique temp -> atomic
        finally:
            try:
                os.unlink(tmp)
            except OSError:
                pass

    def og_msa(self, og: str) -> dict:
        og_idx = self._og_to_idx.get(og)
        if og_idx is None:
            raise HTTPException(status_code=404, detail=f"unknown orthogroup: {og!r}")
        text, source = self._msa_text(og, og_idx)
        _, _, raw_to_display, expected = self._member_context(og_idx)
        pairs = _parse_fasta_text(text)
        matched, _, _ = _match_headers([h for h, _ in pairs], expected)
        rows = []
        for header, seq in pairs:
            match = matched.get(header)
            if match:
                genome, raw = match
                display = raw_to_display.get(raw, raw)
            else:
                genome, raw, display = None, header, header
            rows.append({
                "genome": genome, "row": self._row_of(genome),
                "gene": display, "raw": raw, "seq": seq,
            })
        rows.sort(key=lambda r: (r["row"] if r["row"] >= 0 else 1 << 30, r["gene"]))
        return {
            "orthogroup": og,
            "source": source,
            "length": len(rows[0]["seq"]) if rows else 0,
            "trimmed": self._trimmed_value(og_idx),
            "rows": rows,
        }

    def _gene_tree_index(self) -> dict:
        if self._gt_index is not None:
            return self._gt_index
        with self._gt_lock:
            if self._gt_index is None:
                index = {}
                offset = 0
                with open(self.gene_tree_file, "rb") as fh:
                    for line in fh:
                        if b":" in line:
                            og = line.split(b":", 1)[0].strip().decode()
                            index[og] = (offset, len(line))
                        offset += len(line)
                self._gt_index = index
        return self._gt_index

    def _gene_tree_raw(self, og: str):
        if self.gene_tree_file is None or not self.gene_tree_file.is_file():
            return None
        record = self._gene_tree_index().get(og)
        if record is None:
            return None
        with open(self.gene_tree_file, "rb") as fh:
            fh.seek(record[0])
            line = fh.read(record[1]).decode()
        _, sep, newick = line.strip().partition(":")
        return newick.strip() if sep else None

    def _duplications(self, og: str) -> dict:
        if self._dup_index is None:
            with self._dup_lock:
                if self._dup_index is None:
                    index: dict[str, dict[str, float]] = {}
                    path = (self.results_dir / "Gene_Duplication_Events" / "Duplications.tsv"
                            if self.results_dir else None)
                    if path is not None and path.is_file():
                        labels, rows = read_matrix(path)
                        node_col = labels.index("Gene Tree Node") if "Gene Tree Node" in labels else None
                        support_col = labels.index("Support") if "Support" in labels else None
                        if node_col is not None:
                            for og_id, cells in rows:
                                node = cells[node_col]
                                if not node:
                                    continue
                                support = None
                                if support_col is not None and support_col < len(cells):
                                    try:
                                        support = float(cells[support_col])
                                    except ValueError:
                                        support = None
                                index.setdefault(og_id, {}).setdefault(node, support)
                    self._dup_index = index
        return self._dup_index.get(og, {})

    def og_gene_tree(self, og: str) -> dict:
        og_idx = self._og_to_idx.get(og)
        if og_idx is None:
            raise HTTPException(status_code=404, detail=f"unknown orthogroup: {og!r}")
        raw = self._gene_tree_raw(og)
        if raw is None:
            raise HTTPException(
                status_code=404,
                detail="no gene tree for this orthogroup (OrthoFinder builds trees "
                       "only for orthogroups with at least 4 genes)",
            )
        root = _parse_newick(raw)
        _, _, raw_to_display, expected = self._member_context(og_idx)
        dup = self._duplications(og)
        return self._gene_tree_arrays(og, root, expected, raw_to_display, dup)

    def _gene_tree_arrays(self, og, root, expected, raw_to_display, dup) -> dict:
        parent, length, label, gene, genome_row, duplication, support = (
            [], [], [], [], [], [], [])
        stack = [(root, -1, 0.0)]
        while stack:
            node, par, branch = stack.pop()
            i = len(parent)
            parent.append(par)
            length.append(branch)
            lab = node.label or ""
            label.append(lab)
            # A node is a duplication iff the Duplications table names it; a
            # missing Support value leaves the flag true but support null.
            duplication.append(lab in dup)
            value = dup.get(lab)
            support.append(float(value) if value is not None else None)
            if node.children:
                gene.append(None)
                genome_row.append(-1)
                for child in reversed(node.children):
                    stack.append((
                        child, i,
                        float(child.branch) if child.branch is not None else 0.0))
            else:
                cands = expected.get(lab)
                if cands and len(cands) == 1:
                    gname, raw = next(iter(cands))
                    gene.append(raw_to_display.get(raw, raw))
                    genome_row.append(self._row_of(gname))
                else:
                    gene.append(None)
                    genome_row.append(-1)
        return {
            "orthogroup": og, "parent": parent, "length": length, "label": label,
            "gene": gene, "genome_row": genome_row,
            "duplication": duplication, "support": support,
        }

    def _download_name(self, node, expected, raw_to_display, ids) -> str:
        lab = node.label or ""
        cands = expected.get(lab)
        if cands and len(cands) == 1:
            _, raw = next(iter(cands))
            return raw if ids == "raw" else raw_to_display.get(raw, raw)
        return lab

    def download_fasta(self, og: str, aligned: bool, ids: str) -> tuple[str, str]:
        og_idx = self._og_to_idx.get(og)
        if og_idx is None:
            raise HTTPException(status_code=404, detail=f"unknown orthogroup: {og!r}")
        if aligned:
            data = self.og_msa(og)
            body = "".join(
                f">{r['raw'] if ids == 'raw' else r['gene']}\n{r['seq']}\n"
                for r in data["rows"]
            )
            return body, f"{og}.aligned.fa"
        seqs = self.og_sequences(og_idx)
        if not seqs:
            raise HTTPException(status_code=404, detail="no sequences available")
        body = "".join(
            f">{raw if ids == 'raw' else display}\n{seq}\n"
            for raw, display, _, seq in seqs
        )
        return body, f"{og}.fa"

    def download_newick(self, og: str, ids: str) -> tuple[str, str]:
        og_idx = self._og_to_idx.get(og)
        if og_idx is None:
            raise HTTPException(status_code=404, detail=f"unknown orthogroup: {og!r}")
        raw = self._gene_tree_raw(og)
        if raw is None:
            raise HTTPException(status_code=404, detail="no gene tree for this orthogroup")
        root = _parse_newick(raw)
        _, _, raw_to_display, expected = self._member_context(og_idx)
        name = lambda node: self._download_name(node, expected, raw_to_display, ids)
        return _to_newick(root, name), f"{og}.nwk"

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

        def work(r0, r1):
            dense = np.unpackbits(
                self.presence_tree[r0:r1], axis=1, bitorder="big")[:, : self.n_og]
            out[r0:r1, : self.nbytes] = np.packbits(dense[:, perm], axis=1, bitorder="big")

        ranges = [
            (r0, min(r0 + _ROW_CHUNK, self.n_genomes))
            for r0 in range(0, self.n_genomes, _ROW_CHUNK)
        ]
        if len(ranges) <= 1 or _MATRIX_WORKERS <= 1:
            for r0, r1 in ranges:
                work(r0, r1)
        else:
            list(_matrix_pool().map(lambda ab: work(*ab), ranges))
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
    def matrix_counts(self, order, og0, og1, g0, g1, w, h) -> np.ndarray:
        """Present-cell count per output bin, ``(h, w)`` int32.

        The cumsum/popcount are restricted to the window's byte span and split
        across row bins so the work runs on several cores; the output is
        bit-identical to the whole-matrix version.
        """
        packed = self._packed(order)
        col_abs = og0 + _edges(og1 - og0, w)
        row_local = _edges(g1 - g0, h)
        b = col_abs >> 3
        r = col_abs & 7
        bmin, bmax = int(b[0]), int(b[-1])
        br = b - bmin
        span = bmax - bmin + 1
        H = g1 - g0

        def work(ha, hb):
            r0, r1 = int(row_local[ha]), int(row_local[hb])
            sub = packed[g0 + r0:g0 + r1, bmin:bmax + 1]
            rl = row_local[ha:hb + 1] - r0
            cs = np.zeros((sub.shape[0], span + 1), np.int32)
            np.cumsum(np.bitwise_count(sub), axis=1, dtype=np.int32, out=cs[:, 1:])
            counts = np.bitwise_count(sub[:, br] & MASKS[r])  # (rows, w+1)
            prefix = cs[:, br] + counts
            colsum = prefix[:, 1:] - prefix[:, :-1]
            return ha, hb, np.add.reduceat(colsum, rl[:-1], axis=0, dtype=np.int32)

        tasks = min(_MATRIX_WORKERS, h)
        if tasks <= 1 or H < 512:
            return work(0, h)[2]
        bounds = np.arange(tasks + 1) * h // tasks
        out = np.empty((h, w), np.int32)
        futures = [
            _matrix_pool().submit(work, int(bounds[i]), int(bounds[i + 1]))
            for i in range(tasks) if bounds[i + 1] > bounds[i]
        ]
        for ha, hb, part in (f.result() for f in futures):
            out[ha:hb] = part
        return out

    def matrix(self, order, og0, og1, g0, g1, w, h) -> np.ndarray:
        rowsum = self.matrix_counts(order, og0, og1, g0, g1, w, h)
        row_sizes = np.diff(_edges(g1 - g0, h))[:, None]
        col_sizes = np.diff(_edges(og1 - og0, w))[None, :]
        frac = rowsum / (row_sizes * col_sizes)
        return np.rint(frac * 255.0).clip(0, 255).astype("<u1")

    def _inverse_perm(self, order) -> np.ndarray:
        perm = self.order_perm(order)
        inv = np.empty(self.n_og, np.int64)
        inv[perm] = np.arange(self.n_og)
        return inv

    def matrix_copies(self, order, og0, og1, g0, g1, w, h) -> np.ndarray:
        """Copy-number colour: 10*mean over present cells; 0 = no present cell."""
        counts = self.matrix_counts(order, og0, og1, g0, g1, w, h)
        if not counts.any():
            return np.zeros((h, w), np.uint8)
        # float64 on purpose: float32 shifts values across .5 rounding boundaries
        counts = counts.astype(np.float64)
        extra = np.zeros((h, w), np.float64)
        if self.mc_og.size:
            inv = self._inverse_perm(order)
            pos = inv[self.mc_og]
            row = self.mc_row
            sel = (pos >= og0) & (pos < og1) & (row >= g0) & (row < g1)
            if sel.any():
                cb = np.searchsorted(_edges(og1 - og0, w), pos[sel] - og0, side="right") - 1
                rb = np.searchsorted(_edges(g1 - g0, h), row[sel] - g0, side="right") - 1
                extra = np.bincount(
                    rb * w + cb, weights=self.mc_extra[sel], minlength=h * w
                ).reshape(h, w)
        return _encode_copies(counts, extra)

    # --- annotation track ---------------------------------------------------
    def annotation_track(self, order, og0, og1, w, name) -> np.ndarray:
        """Most common category per bin, as ``1 + category index`` (0 = missing)."""
        if name not in self._ann_codes:
            raise HTTPException(
                status_code=400, detail=f"not a categorical annotation: {name!r}"
            )
        perm = self.order_perm(order)
        codes = self._ann_codes[name][perm[og0:og1]]
        edges = _edges(og1 - og0, w)
        n_cat = len(self._ann_categories[name])
        out = np.zeros(w, np.uint16)
        for i in range(w):
            seg = codes[edges[i]:edges[i + 1]]
            if seg.size:
                counts = np.bincount(seg, minlength=n_cat + 1)[1:]
                if counts.any():
                    out[i] = int(np.argmax(counts)) + 1
        return out

    # --- search -------------------------------------------------------------
    @staticmethod
    def _field_mask(frame, name, q) -> np.ndarray:
        """Case-insensitive literal substring mask for one column.

        Nulls are filled to ``False`` before the numpy conversion, so an
        all-null annotation column masks to all-``False``.
        """
        return frame.select(
            pl.col(name).str.to_lowercase().str.contains(q, literal=True).fill_null(False)
        )[name].to_numpy()

    def _hit_idxs(self, q) -> np.ndarray:
        if not q:
            return np.empty(0, np.int64)
        q = q.lower()
        hit = np.zeros(self.n_og, bool)
        for name, frame in self._search_fields:
            hit |= self._field_mask(frame, name, q)
        return np.flatnonzero(hit)

    def search(self, q, order, limit) -> dict:
        inv = self._inverse_perm(order)
        og_ids = self.orthogroups["orthogroup"].to_numpy()
        hits = []
        if q:
            q = q.lower()
            for name, frame in self._search_fields:
                idxs = np.flatnonzero(self._field_mask(frame, name, q))
                if idxs.size == 0:
                    continue
                values = (og_ids if name == "orthogroup"
                          else frame[name].to_numpy())
                for i in idxs:
                    hits.append({
                        "orthogroup": og_ids[i],
                        "og_idx": int(i),
                        "position": int(inv[i]),
                        "field": name,
                        "text": values[i],
                    })
        hits.sort(key=lambda h: h["position"])
        return {"hits": hits[:limit], "total": len(hits)}

    def search_perm(self, q) -> np.ndarray:
        hit = np.zeros(self.n_og, bool)
        hit[self._hit_idxs(q)] = True
        mask = hit[self.perm_tree]
        return np.concatenate([self.perm_tree[mask], self.perm_tree[~mask]]).astype(np.int64)

    def pattern_perm(self, og_idx) -> np.ndarray:
        mask = self.pattern_ids == self.pattern_ids[og_idx]
        ordered = mask[self.perm_tree]
        return np.concatenate([self.perm_tree[ordered], self.perm_tree[~ordered]]).astype(np.int64)

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
            "annotations": self.annotation_specs,
            "genome_metadata": self.gm_specs,
            "detail": {
                "msa": bool(self.detail.get("msa")),
                "gene_trees": bool(self.detail.get("gene_trees")),
                "sequences": self._has_sequences,
                "famsa": self.famsa_ok,
            },
        }


def _encode_copies(counts: np.ndarray, extra: np.ndarray) -> np.ndarray:
    """uint8 = round(10 * (present + extra) / present); 0 where nothing is present."""
    counts = np.asarray(counts, np.float64)
    mean = np.divide(counts + extra, counts, out=np.zeros(counts.shape), where=counts > 0)
    return np.rint(mean * 10.0).clip(0, 255).astype(np.uint8)


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


def _annotation_values(pg, og_idx: int) -> dict:
    """``{annotation name: value}`` for one OG (row index == og_idx)."""
    if pg.annotations is None:
        return {}
    row = pg.annotations.row(og_idx, named=True)
    return {name: (None if row[name] is None else str(row[name]))
            for name in pg.annotation_names}


def create_app(data_dir, static_dir=None, famsa="famsa") -> FastAPI:
    pg = Pangenome(data_dir, famsa=famsa)
    app = FastAPI(title="OrthoPan", docs_url=None, redoc_url=None)
    app.state.pg = pg

    @app.get("/api/meta")
    def meta():
        return pg.meta()

    @app.get("/api/genomes")
    def genomes():
        # Built and serialized once at startup; just stream the cached bytes.
        return Response(content=pg._genomes_json, media_type="application/json")

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
        value: str = "presence",
    ):
        if value not in ("presence", "copies"):
            raise HTTPException(status_code=400, detail=f"unknown value mode: {value!r}")
        og0, og1, g0, g1 = _window(pg.n_og, pg.n_genomes, og0, og1, g0, g1, w, h)
        if value == "copies":
            out = pg.matrix_copies(order, og0, og1, g0, g1, w, h)
        else:
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

    @app.get("/api/annotation-track")
    def annotation_track(
        order: str = "tree",
        og0: int = 0,
        og1: int | None = None,
        w: int = Query(...),
        name: str = Query(...),
    ):
        og1 = pg.n_og if og1 is None else og1
        if not (0 <= og0 < og1 <= pg.n_og):
            raise HTTPException(status_code=400, detail="window out of range")
        if not (1 <= w <= MAX_CELLS):
            raise HTTPException(status_code=400, detail=f"w must be in 1..{MAX_CELLS}")
        if w > og1 - og0:
            raise HTTPException(status_code=400, detail="w exceeds window size")
        out = pg.annotation_track(order, og0, og1, w, name)
        return Response(
            content=out.astype("<u2").tobytes(),
            media_type="application/octet-stream",
            headers={"X-Width": str(w), "X-Categories": name},
        )

    @app.get("/api/search")
    def search(
        q: str = Query(...),
        order: str = "tree",
        limit: int = 50,
    ):
        limit = max(1, min(limit, 500))
        return pg.search(q, order, limit)

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
        annotations = _annotation_values(pg, og_idx)
        members = [
            {
                "genome": genome,
                "genes": [{"id": g, "raw": r} for g, r in zip(genes, raw)],
            }
            for _, genome, genes, raw in pg.og_members(og_idx)
        ]
        body = {
            "orthogroup": row["orthogroup"],
            "og_idx": og_idx,
            "n_genomes": int(row["n_genomes"]),
            "n_genes": int(row["n_genes"]),
            "mrca_label": row["mrca_label"],
            "purity": None if row["purity"] is None else float(row["purity"]),
            "tracks": tracks,
            "annotations": annotations,
            "pattern_id": int(row["pattern_id"]),
            "n_same_pattern": int(row["n_same_pattern"]),
            "members": members,
        }
        # Serialize here: FastAPI's jsonable_encoder is very slow on 10k members.
        return Response(
            content=json.dumps(body, ensure_ascii=False, separators=(",", ":")),
            media_type="application/json",
        )

    @app.get("/api/orthogroup/{og}/pattern")
    def orthogroup_pattern(og: str, limit: int = 200):
        og_idx = pg._og_to_idx.get(og)
        if og_idx is None:
            raise HTTPException(status_code=404, detail=f"unknown orthogroup: {og!r}")
        limit = max(1, min(limit, MAX_LIST))
        pid = int(pg.pattern_ids[og_idx])
        og_ids = pg.orthogroups["orthogroup"].to_numpy()
        pids = pg.pattern_ids
        idxs = np.flatnonzero(pids == pid)[:limit]
        orthogroups_out = [
            {
                "orthogroup": og_ids[i],
                "og_idx": int(i),
                "annotations": _annotation_values(pg, int(i)),
            }
            for i in idxs
        ]
        return {
            "pattern_id": pid,
            "total": int(pg.n_same_pattern[og_idx]),
            "orthogroups": orthogroups_out,
        }

    @app.get("/api/orthogroup/{og}/msa")
    def orthogroup_msa(og: str):
        return pg.og_msa(og)

    @app.get("/api/orthogroup/{og}/gene-tree")
    def orthogroup_gene_tree(og: str):
        return pg.og_gene_tree(og)

    @app.get("/api/orthogroup/{og}/fasta")
    def orthogroup_fasta(og: str, aligned: int = 0, ids: str = "display"):
        if aligned not in (0, 1):
            raise HTTPException(status_code=400, detail="aligned must be 0 or 1")
        if ids not in ("display", "raw"):
            raise HTTPException(status_code=400, detail="ids must be display or raw")
        body, filename = pg.download_fasta(og, bool(aligned), ids)
        return Response(
            content=body,
            media_type="text/plain",
            headers={"Content-Disposition": _content_disposition(filename)},
        )

    @app.get("/api/orthogroup/{og}/gene-tree.nwk")
    def orthogroup_newick(og: str, ids: str = "display"):
        if ids not in ("display", "raw"):
            raise HTTPException(status_code=400, detail="ids must be display or raw")
        body, filename = pg.download_newick(og, ids)
        return Response(
            content=body,
            media_type="text/plain",
            headers={"Content-Disposition": _content_disposition(filename)},
        )

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
        if not isinstance(body, dict):
            raise HTTPException(status_code=400, detail="invalid JSON body")
        by = body.get("by")
        if by == "contrast":
            names = body.get("genomes")
            if not isinstance(names, list) or not names or not all(isinstance(n, str) for n in names):
                raise HTTPException(
                    status_code=400, detail="'genomes' must be a non-empty list of strings")
            if len(names) > pg.n_genomes:
                raise HTTPException(status_code=400, detail="more genomes than exist")
            unique = list(dict.fromkeys(names))  # canonicalize duplicates, keep order
            digest = hashlib.sha1(json.dumps(sorted(set(names))).encode()).hexdigest()[:12]
            order_id = "contrast:" + digest
        elif by == "search":
            q = body.get("q")
            if not isinstance(q, str) or not q:
                raise HTTPException(status_code=400, detail="'q' must be a non-empty string")
            digest = hashlib.sha1(json.dumps(q).encode()).hexdigest()[:12]
            order_id = "search:" + digest
        elif by == "pattern":
            og = body.get("og")
            if not isinstance(og, str):
                raise HTTPException(status_code=400, detail="'og' must be a string")
            if og not in pg._og_to_idx:
                raise HTTPException(status_code=404, detail=f"unknown orthogroup: {og!r}")
            digest = hashlib.sha1(og.encode()).hexdigest()[:12]
            order_id = "pattern:" + digest
        else:
            raise HTTPException(status_code=400, detail="unknown order type")
        if pg._dynamic_get(order_id) is not None:
            return {"order": order_id}  # reuse the existing (possibly warm) entry
        if by == "contrast":
            perm = pg.contrast_perm(unique)
        elif by == "search":
            perm = pg.search_perm(q)
        else:
            perm = pg.pattern_perm(pg._og_to_idx[og])
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
        packaged = Path(__file__).resolve().parent / "static"
        if packaged.is_dir():
            static_dir = packaged
        else:
            candidate = Path(__file__).resolve().parents[2] / "web" / "dist"
            static_dir = candidate if candidate.is_dir() else None
    if static_dir is not None:
        app.mount("/", StaticFiles(directory=str(static_dir), html=True), name="static")

    return app


def run(args) -> None:
    app = create_app(args.data_dir, args.static, famsa=args.famsa)
    uvicorn.run(app, host=args.host, port=args.port)


def add_parser(subparsers):
    p = subparsers.add_parser("serve", help="serve an imported pangenome directory")
    p.add_argument("data_dir", metavar="DIR", help="directory written by `orthopan import`")
    p.add_argument("--host", default="127.0.0.1", help="bind host (default 127.0.0.1)")
    p.add_argument("--port", type=int, default=8000, help="bind port (default 8000)")
    p.add_argument("--static", metavar="PATH", help="built frontend directory to serve at /")
    p.add_argument("--famsa", default="famsa", metavar="PATH",
                   help="famsa binary for on-demand alignment (default: famsa on PATH)")
    p.set_defaults(func=run)
