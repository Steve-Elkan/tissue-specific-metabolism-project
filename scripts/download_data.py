
import requests
import zipfile
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

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
    """Download a file if it does not already exist."""

    output = DATA_DIR / filename

    if output.exists() and output.stat().st_size > 0:
        print(f"Already downloaded: {filename}")
        return output

    print(f"Downloading {filename}...")

    temporary_file = output.with_suffix(output.suffix + ".part")

    try:
        with requests.get(
            url,
            stream=True,
            timeout=(30, 180)
        ) as response:

            response.raise_for_status()

            with open(temporary_file, "wb") as file:
                for chunk in response.iter_content(
                    chunk_size=1024 * 1024
                ):
                    if chunk:
                        file.write(chunk)

        temporary_file.replace(output)

    except Exception:
        temporary_file.unlink(missing_ok=True)
        raise

    print(f"Saved: {output}")
    return output


def extract_expression(zip_path):
    """Safely extract the expression dataset."""

    print("Extracting expression dataset...")

    with zipfile.ZipFile(zip_path, "r") as archive:

        for member in archive.infolist():

            destination = (DATA_DIR / member.filename).resolve()

            if not destination.is_relative_to(DATA_DIR.resolve()):
                raise ValueError(
                    f"Unsafe archive path: {member.filename}"
                )

        archive.extractall(DATA_DIR)

    print("Expression dataset extracted successfully.")


def main():

    for filename, url in FILES.items():
        download_file(filename, url)

    zip_path = DATA_DIR / "expression.zip"

    if not zipfile.is_zipfile(zip_path):
        raise ValueError("Downloaded expression file is not a valid ZIP.")

    extract_expression(zip_path)

    print("\nAll datasets are ready.")
    print(f"Data directory: {DATA_DIR}")


if __name__ == "__main__":
    main()
