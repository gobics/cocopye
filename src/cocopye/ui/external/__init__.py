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

"""
This module contains functions for checking the applications dependencies (binaries and data) and downloading them if
needed.

Please note that these functions are part of the cli and not intended for use as a library. In particular, the functions
may output to stdout and sterr, expect input from stdin or even terminate the program.

### Submodules:

- **bin**: Functions for checking and downloading Prodigal and UProC
- **data**: Functions for checking and downloading the Pfam database and UProC models
"""
import os
import platform
import stat
import sys
from typing import Tuple, List

import requests
from tqdm import tqdm
from appdirs import user_data_dir

from .. import config
from ..config import change_config
from ... import constants


def check_and_download_dependencies(check_only: bool = False) -> None:
    """
    This function checks if all neccessary dependencies are available. If there's something missing, it will offer the
    user to download the missing files automatically. The function writes some status messages to stdout and waits for
    user input on stdin when asking if it should download the missing parts. Furthermore, it may terminate the program
    if there are corrupted dependencies or errors during download. If this function returns successfully, it is
    (hopefully)  guaranteed that all required files are present and working.
    """
    ok, missing, error = check_dependencies()

    if len(error) > 0:
        # It might be useful to provide automatic solutions for common errors (at least try to redownload the file or
        # provide a more useful error message).
        print("\nError: Some dependencies were found, but they don't seem to work. Please check them manually.\n")
        sys.exit(1)

    if len(missing) > 0:
        if check_only:
            print("\nSome dependencies are missing.")
            print("If you have already downlaoded them, you can specify their path in the configuration file.")
            print("You can also use 'cocopye toolbox download-dependencies' to automatically download missing files.")
            sys.exit(1)
        if config.ARGS.offline:
            print("You are in offline mode. Disable it if you want to use the automatic download.")
            sys.exit(1)

        download_dependencies(missing)

        print("Downloads finished. Rechecking dependencies...\n")
        ok, missing, error = check_dependencies()
        if len(missing) > 0 or len(error) > 0:
            print("There are still errors or missing files. This should not happen and you have likely encountered a")
            print("bug in the application. Exiting.")
            sys.exit(1)

        print("Check successful.\n")


def check_dependencies() -> Tuple[List[str], List[str], List[str]]:
    """
    This function checks the status of dependencies. It writes an overview of the result to stdout. This function is
    mainly for use by `check_and_download_dependencies`.

    :return: 3-tuple (ok, missing, error) where each element is a list of dependencies that have the respective status
    Example: `(["prodigal", "uproc"], ["pfam"], ["model"])` means that Prodigal and UProC are present and working, the
    pfam database is missing and the model is there but corrupted.
    """
    from .bin import check_uproc
    from .data import check_pfam_db, check_model, check_cocopye_db

    checks = [
        check_uproc(config.CONFIG["external"]["uproc_prot_bin"]),
        check_pfam_db(config.CONFIG["external"]["uproc_pfam_db"], "24" if config.ARGS.pfam24 else "28"),
        check_model(config.CONFIG["external"]["uproc_models"]),
        check_cocopye_db(config.CONFIG["external"]["cocopye_db"], config.ARGS.offline)
    ]

    status: Tuple[List[str], List[str], List[str]] = ([], [], [])  # ok, missing, error

    print("External dependencies:")
    for lst, identifier, output in checks:
        print(output)
        status[lst].append(identifier)
    print()

    if "outdated" in checks[3][2]:
        print("New CoCoPyE database release available.")
        print("You can update the database by running \"cocopye toolbox --update-database\".\n")

    return status[0], status[1], status[2]


def download_dependencies(missing: List[str]) -> None:
    """
    This tries to download and install missing dependencies. During the process it prints a bunch of status and progress
    messages to stdout. If the detected operating system is not Windows or Linux it will terminate the application.
    This function is mainly for use by `check_and_download_dependencies`.

    :param missing: A list of missing dependencies as generated by `check_dependencies`
    """
    from .bin import build_uproc_prot, download_uproc_win
    from .data import download_pfam_db, download_model, download_cocopye_db

    if "uproc" in missing:
        opsys = platform.system()
        if opsys == "Linux":
            build_uproc_prot(
                constants.UPROC["SRC"],
                os.path.join(user_data_dir("cocopye"), "uproc"),
                config.ARGS.verbose
            )

            change_config(
                "external", "uproc_prot_bin",
                os.path.join(user_data_dir("cocopye"), "uproc", "bin", "uproc-prot")
            )

            change_config(
                "external", "uproc_import_bin",
                os.path.join(user_data_dir("cocopye"), "uproc", "bin", "uproc-import")
            )

            change_config(
                "external", "uproc_orf_bin",
                os.path.join(user_data_dir("cocopye"), "uproc", "bin", "uproc-orf")
            )
        elif opsys == "Windows":
            download_uproc_win(constants.UPROC["WIN"], os.path.join(user_data_dir("cocopye"), "uproc"))

            change_config(
                "external", "uproc_prot_bin",
                os.path.join(user_data_dir("cocopye"), "uproc", "uproc-prot.exe")
            )

            change_config(
                "external", "uproc_import_bin",
                os.path.join(user_data_dir("cocopye"), "uproc", "uproc-import.exe")
            )

            change_config(
                "external", "uproc_orf_bin",
                os.path.join(user_data_dir("cocopye"), "uproc", "uproc-orf.exe")
            )
        else:
            print("Automatic installation of UProC is currently only supported on Windows and Linux.")
            print("See http://uproc.gobics.de for more information on how to install UProC on your system.")
            print("Exiting.\n")
            sys.exit(1)

    if "pfam" in missing:
        download_pfam_db(
            constants.PFAM_DB,
            config.CONFIG["external"]["uproc_import_bin"],
            24 if config.ARGS.pfam24 else 28,
            config.ARGS.verbose
        )

        change_config(
            "external", "uproc_pfam_db",
            os.path.join(user_data_dir("cocopye"), "pfam_db")
        )

    if "model" in missing:
        download_model(constants.UPROC_MODEL)

        change_config(
            "external", "uproc_models",
            os.path.join(user_data_dir("cocopye"), "model")
        )

    if "cocopye_db" in missing:
        download_cocopye_db(constants.COCOPYE_DB)

        change_config(
            "external", "cocopye_db",
            os.path.join(user_data_dir("cocopye"), "cocopye_db")
        )


def download(url: str, dirname: str, fname: str, label: str, chunk_size: int = 1024) -> None:
    """
    An auxiliary function to download a file. It shows a progress bar which gets removed once the download is complete.

    :param url: Download URL
    :param dirname: Destination directory (without filename)
    :param fname: Destination file (without directory path)
    :param label: Progress bar label
    :param chunk_size: Download chunk size
    """
    os.makedirs(dirname, exist_ok=True)

    # Adapted from https://stackoverflow.com/questions/37573483/progress-bar-while-download-file-over-http-with- \
    # requests/62113263#62113263
    # Author: Yan QiDong
    # License: CC BY-SA 4.0
    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    with open(os.path.join(dirname, fname), 'wb') as file, tqdm(
        desc=label,
        total=total,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
        ncols=100,
        leave=False
    ) as bar:
        for content in resp.iter_content(chunk_size=chunk_size):
            size = file.write(content)
            bar.update(size)

    # This is probably only necessary for Prodigal, but I assume it won't hurt in other cases
    os.chmod(os.path.join(dirname, fname), stat.S_IRUSR ^ stat.S_IWUSR ^ stat.S_IXUSR)


def _green(s: str) -> str:
    return "\033[92m" + s + "\033[0m"


def _yellow(s: str) -> str:
    return "\033[93m" + s + "\033[0m"


def _red(s: str) -> str:
    return "\033[91m" + s + "\033[0m"


_TICK = _green("✓")
_CROSS = _red("✗")
