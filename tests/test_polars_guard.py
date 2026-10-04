"""Guard against the broken polars Series dispatch on Python 3.15.

In this environment ``pl.Series(...).str`` / ``.list`` / ``.fill_null`` /
``.is_null`` silently return ``None`` (a bytecode check in polars'
``expr_dispatch``), so the source must use DataFrame *expressions* instead.
This test fails with an explanation if any such Series call appears in ``src/``.
The scan is best-effort: it only catches literal ``pl.Series(...).<method>`` chains,
not methods called on a Series held in a variable.
"""

import re
from pathlib import Path

import polars as pl

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"

# pl.Series(...) followed by an affected method/namespace access.
OFFENDING = re.compile(r"pl\.Series\([^)]*\)\s*\.(?:str|list|fill_null|is_null)\b")


def test_report_polars_series_dispatch():
    # Informational: skips once polars dispatch works again (guard can then be dropped).
    if pl.Series(["A"]).str.to_lowercase() is not None:
        import pytest

        pytest.skip("polars Series dispatch works; this guard can be removed")


def test_src_avoids_broken_polars_series_methods():
    offenders = []
    for path in SRC.rglob("*.py"):
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            if OFFENDING.search(line):
                offenders.append(f"{path.relative_to(REPO)}:{lineno}: {line.strip()}")
    assert not offenders, (
        "polars Series methods silently return None on Python 3.15; "
        "use DataFrame expressions instead. Offending lines:\n" + "\n".join(offenders)
    )
