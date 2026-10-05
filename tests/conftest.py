from pathlib import Path
import os
import urllib.request
import pytest

DATA_DIR = Path(__file__).parent / "data"

REMOTE_FILES = {
    "orthologs.tsv": "https://demo.abrinca.com/files_html/orthologs/orthologs.tsv?download",
    "orthologs-descriptions.tsv": "https://demo.abrinca.com/files_html/orthologs/HOG_best_name.tsv?download",
}


@pytest.fixture(scope="session")
def test_data():
    if os.environ.get("ORTHOPAN_TEST_NO_RESOURCES") == "1":
        pytest.skip("ORTHOPAN_TEST_NO_RESOURCES set: tests/data unavailable")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    downloaded_paths = {}

    for filename, url in REMOTE_FILES.items():
        dest = DATA_DIR / filename
        if not dest.is_file():
            req = urllib.request.Request(url, headers={"User-Agent": "OrthoPan-TestRunner"})
            tmp = dest.with_name(dest.name + ".part")
            try:
                with urllib.request.urlopen(req) as resp, open(tmp, "wb") as f:
                    f.write(resp.read())
            except OSError as exc:
                tmp.unlink(missing_ok=True)
                pytest.skip(
                    f"tests/data is unavailable (download of {filename} failed: {exc}); "
                    "run with network access once to cache it in tests/data/"
                )
            tmp.replace(dest)
        downloaded_paths[filename] = dest

    return downloaded_paths
