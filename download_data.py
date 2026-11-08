import requests
import zipfile
from pathlib import Path

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

FILES = {
    "expression.zip": (
        "https://www.proteinatlas.org/download/"
        "rna_tissue_consensus_tissues.tsv.zip"
    ),
    "Recon3D.xml": (
        "https://bigg.ucsd.edu/static/models/Recon3D.xml"
    ),
}


def download_file(filename, url):
    """Download a dataset if it is not already available."""

    output = DATA_DIR / filename

    if output.exists():
        print(f"Already downloaded: {filename}")
        return output

    print(f"Downloading {filename}...")

    response = requests.get(url, timeout=180)
    response.raise_for_status()

    output.write_bytes(response.content)

    print(f"Saved: {output}")
    return output


def main():
    for filename, url in FILES.items():
        download_file(filename, url)

    zip_path = DATA_DIR / "expression.zip"

    with zipfile.ZipFile(zip_path, "r") as archive:
        archive.extractall(DATA_DIR)

    print("Dataset download completed.")


if __name__ == "__main__":
    main()
