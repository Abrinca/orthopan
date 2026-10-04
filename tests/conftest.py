from pathlib import Path
import urllib.request
import pytest

DATA_DIR = Path(__file__).parent / "data"

REMOTE_FILES = {
    "N0.tsv": "https://opengenomebrowser.bioinformatics.unibe.ch/files_html/orthologs/N0.tsv",
    "N0.descriptions.tsv": "https://opengenomebrowser.bioinformatics.unibe.ch/files_html/annotation-descriptions/OL.tsv",
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