import gzip
import os
import shutil
import subprocess
import tarfile
import tempfile

from appdirs import user_cache_dir, user_data_dir

from ..external import download, _red, _green, _TICK, _CROSS


def check_pfam_db(pfam_dir: str):
    result = _check_folder(pfam_dir, ["fwd.ecurve", "idmap", "prot_thresh_e2", "prot_thresh_e3", "rev.ecurve"])

    if result == "found":
        return 0, "pfam", "  " + _TICK + " Pfam database\t" + _green(result)
    else:
        return 1 if result == "mising" else 2, "pfam", _CROSS + " Pfam database\t" + _red(result)


def check_model(model_dir: str):
    result = _check_folder(model_dir,
                           ["aa_probs", "alphabet", "codon_scores", "orf_thresh_e1", "orf_thresh_e2", "substmat"])

    if result == "found":
        return 0, "model", "  " + _TICK + " Models\t\t" + _green(result)
    else:
        return 1 if result == "missing" else 2, "model", _CROSS + " Models\t\t" + _red(result)


def _check_folder(folder: str, files):
    try:
        content = os.listdir(folder)
    except FileNotFoundError:
        return "not found"

    return "found" if all([file in content for file in files]) else "error"


def download_pfam_db(url: str, import_bin: str):
    # not using /tmp, because of the large file size
    with tempfile.TemporaryDirectory(prefix="libcc_", dir=user_cache_dir(None)) as tmpdir:
        download(
            url,
            tmpdir,
            "pfam.uprocdb.gz",
            "- Downloading UProC Pfam database"
        )
        print("- Downloading UProc Pfam database ✓")

        print("- Extracting database. This may take a while.", end="", flush=True)
        with gzip.open(os.path.join(tmpdir, "pfam.uprocdb.gz"), "rb") as f_in:
            with open(os.path.join(tmpdir, "pfam.uprocdb"), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        print("\r- Extracting database ✓                                      ")

        print("- Importing database. This may take a while.", end="", flush=True)
        uproc_import = subprocess.Popen([import_bin,
                                         os.path.join(tmpdir, "pfam.uprocdb"),
                                         os.path.join(user_data_dir("libcc"), "pfam_db")],
                                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        uproc_import.wait()
        print("\r- Importing database ✓                             \n")


def download_model(url: str):
    with tempfile.TemporaryDirectory() as tmpdir:
        download(
            url,
            tmpdir,
            "model.tar.gz",
            "- Downloading UProC Model"
        )
        print("- Downloading UProc Model ✓")

        print("- Extracting UProC model")
        tar = tarfile.open(os.path.join(tmpdir, "model.tar.gz"))
        tar.extractall(user_data_dir("libcc"))
        print("\r- Extracting UProC model ✓\n")
