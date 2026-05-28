import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

URL = "https://ndownloader.figshare.com/files/55209713"
DB_TARBALL = Path("db.tar.gz")
USER_AGENT = "Mozilla/5.0"
MAX_ATTEMPTS = 30
RETRY_WAIT_SECONDS = 20


def download_database():
    """Download the database tarball from FigShare, validating and retrying."""
    print("Downloading database from FigShare...")
    request = urllib.request.Request(URL, headers={"User-Agent": USER_AGENT})
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(request) as response, DB_TARBALL.open("wb") as handle:
                shutil.copyfileobj(response, handle)
            if tarfile.is_tarfile(DB_TARBALL):
                print(f"Database downloaded (attempt {attempt}).")
                return
            reason = "response was not a valid tarball (FigShare may be staging the file)"
        except urllib.error.URLError as err:
            reason = str(err)
        if attempt < MAX_ATTEMPTS:
            print(
                f"Attempt {attempt}/{MAX_ATTEMPTS} failed: {reason}; "
                f"retrying in {RETRY_WAIT_SECONDS}s..."
            )
            time.sleep(RETRY_WAIT_SECONDS)
    sys.exit(
        f"ERROR: could not download a valid database from {URL} after "
        f"{MAX_ATTEMPTS} attempts."
    )


def main():
    download_database()

    print("Extracting database")
    with tarfile.open(DB_TARBALL, "r:gz") as tar:
        # filter="data" (Python >= 3.12) blocks unsafe members; older
        # interpreters don't accept the argument.
        try:
            tar.extractall(filter="data")
        except TypeError:
            tar.extractall()
    DB_TARBALL.unlink()
    print("Database complete!")


if __name__ == "__main__":
    main()
