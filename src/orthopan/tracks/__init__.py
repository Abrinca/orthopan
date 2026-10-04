"""Per-orthogroup track preprocessors.

Every track file is a TSV: header `orthogroup<TAB><track>[<TAB><track>...]`,
one row per orthogroup, float values, empty cell = missing.
"""

import math
import sys
from pathlib import Path


def write_tracks(rows: dict[str, dict[str, float]], columns: list[str], out: str | Path | None = None) -> None:
    """Write {orthogroup: {track: value}} as a track TSV (stdout if out is None)."""
    fh = open(out, "w") if out else sys.stdout
    try:
        fh.write("\t".join(["orthogroup", *columns]) + "\n")
        for og, values in rows.items():
            cells = [values.get(c) for c in columns]
            fh.write("\t".join([og, *("" if v is None or not math.isfinite(v) else f"{v:.6g}" for v in cells)]) + "\n")
    finally:
        if out:
            fh.close()
