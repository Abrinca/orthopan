from pathlib import Path
import urllib.request
import pytest

DATA_DIR = Path(__file__).parent / "data"

REMOTE_FILES = {
    "orthologs.tsv": "https://demo.abrinca.com/files_html/orthologs/orthologs.tsv?download",
    "orthologs-descriptions.tsv": "https://demo.abrinca.com/files_html/orthologs/HOG_best_name.tsv?download",
}


@pytest.fixture(scope="session")
def test_data():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    downloaded_paths = {}

    for filename, url in REMOTE_FILES.items():
        dest = DATA_DIR / filename
        if not dest.is_file():
            req = urllib.request.Request(url, headers={"User-Agent": "OrthoPan-TestRunner"})
            with urllib.request.urlopen(req) as resp, open(dest, "wb") as f:
                f.write(resp.read())
        downloaded_paths[filename] = dest

    return downloaded_paths