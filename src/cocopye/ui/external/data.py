# Copyright 2023 Niklas Birth
#
# This file is part of CoCoPyE.
#
# CoCoPyE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CoCoPyE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CoCoPyE. If not, see <https://www.gnu.org/licenses/>.

import gzip
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
import requests
from typing import Tuple, List, Optional

from appdirs import user_cache_dir, user_data_dir
from packaging.version import Version

from ..external import download, _red, _yellow, _green, _TICK, _CROSS
from ... import constants


def check_pfam_db(pfam_dir: str, version: str = "28") -> Tuple[int, str, str]:
    result = _check_folder(os.path.join(pfam_dir, version),
                           ["fwd.ecurve", "idmap", "prot_thresh_e2", "prot_thresh_e3", "rev.ecurve"])

    if result == "found":
        return 0, "pfam", "  " + _TICK + " Pfam database\t" + _green("v" + version)
    else:
        return 1 if result == "not found" else 2, "pfam", "  " + _CROSS + " Pfam database\t" + _red(result)


def check_model(model_dir: str) -> Tuple[int, str, str]:
    result = _check_folder(model_dir,
                           ["aa_probs", "alphabet", "codon_scores", "orf_thresh_e1", "orf_thresh_e2", "substmat"])

    if result == "found":
        return 0, "model", "  " + _TICK + " Models\t\t" + _green(result)
    else:
        return 1 if result == "not found" else 2, "model", "  " + _CROSS + " Models\t\t" + _red(result)


def check_cocopye_db(db_dir: str, offline: bool = False) -> Tuple[int, str, str]:
    db_files = [
        "universal_Bacteria.npy",
        "universal_Archaea.npy",
        "metadata.csv",
        "count_matrix.npz",
        "model_comp.pickle",
        "model_cont.pickle"
    ]

    result_28 = _check_folder(os.path.join(db_dir, "28"), db_files)
    result_24 = _check_folder(os.path.join(db_dir, "24"), db_files)

    if result_24 == "found" and result_28 == "found":
        version = open(os.path.join(db_dir, "version.txt"), "r").read().strip()
        if not offline and new_db_version_available(version):
            return 0, "cocopye_db", _TICK + " CoCoPyE database\t" + _yellow(version + " (outdated)")
        return 0, "cocopye_db", _TICK + " CoCoPyE database\t" + _green(version)
    else:
        return 1 if result_28 == "not found" else 2, "cocopye_db", _CROSS + " CoCoPyE database\t" + _red(result_28)


def new_db_version_available(current_version):
    try:
        version = requests.get(constants.COCOPYE_DB_LATEST_RELEASE).json()["tag_name"]
    except requests.exceptions.ConnectionError:
        return False
    except KeyError:
        return False

    versions = [version[1:], current_version[1:]]
    versions.sort(key=Version)

    return "v" + versions[-1] != current_version


def _check_folder(folder: str, files: List[str]) -> str:
    try:
        content = os.listdir(folder)
    except FileNotFoundError:
        return "not found"

    return "found" if all([file in content for file in files]) else "error"


def download_pfam_db(url: str, import_bin: str, version: int = 28, verbose: bool = False) -> None:
    output = subprocess.DEVNULL if not verbose else None

    # not using /tmp, because of the large file size
    with tempfile.TemporaryDirectory(prefix="cocopye_", dir=user_cache_dir(None)) as tmpdir:
        download(
            url[version],
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
        os.makedirs(os.path.join(user_data_dir("cocopye"), "pfam_db", str(version)), exist_ok=True)
        uproc_import = subprocess.Popen([import_bin,
                                         os.path.join(tmpdir, "pfam.uprocdb"),
                                         os.path.join(user_data_dir("cocopye"), "pfam_db", str(version))],
                                        stdout=output, stderr=output)
        if uproc_import.wait() != 0:
            print("\n\nError while running uproc-import. Rerun the command with '--verbose' for subprocess output.")
            sys.exit(1)

        print("\r- Importing database ✓                             \n")


def download_model(url: str) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        download(
            url,
            tmpdir,
            "model.tar.gz",
            "- Downloading UProC Model"
        )
        print("- Downloading UProc Model ✓")

        print("- Extracting UProC model", end="", flush=True)
        with tarfile.open(os.path.join(tmpdir, "model.tar.gz")) as tar:
            tar.extractall(user_data_dir("cocopye"))
        print("\r- Extracting UProC model ✓\n")


def download_cocopye_db(url: str, db_dir: Optional[str] = None) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        download(
            url,
            tmpdir,
            "cocopye_db.zip",
            "- Downloading CoCoPyE database"
        )
        print("- Downloading CoCoPyE database ✓")

        print("- Extracting database", end="", flush=True)
        with zipfile.ZipFile(os.path.join(tmpdir, "cocopye_db.zip"), 'r') as zip_ref:
            if db_dir is None:
                db_dir = os.path.join(user_data_dir("cocopye"), "cocopye_db")
            zip_ref.extractall(db_dir)
        print("\r- Extracting database ✓\n")


def update_cocopye_db(url: str, db_dir: str) -> None:
    print("- Removing old database", end="")
    shutil.rmtree(db_dir)
    print("\r- Removing old database ✓")

    download_cocopye_db(url, db_dir)

    print("Database update successful.\n")
