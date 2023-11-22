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

import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from typing import Tuple

from ..external import download, _green, _red, _TICK, _CROSS


def check_uproc(uproc_bin: str) -> Tuple[int, str, str]:
    try:
        process = subprocess.Popen(
            [uproc_bin, '-v'],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
        )
    except FileNotFoundError:
        return 1, "uproc", _CROSS + " UProC\t\t\t" + _red("not found")

    assert process.stdout is not None

    version = "v" + process.stdout.read().strip().partition("\n")[0].rpartition("version ")[2]

    return 0, "uproc", _TICK + " UProC\t\t\t" + _green(version)


def build_uproc_prot(url: str, install_dir: str, verbose: bool = False) -> None:
    output = subprocess.DEVNULL if not verbose else None

    with tempfile.TemporaryDirectory() as tmpdir:
        download(url, tmpdir, "uproc.tar.gz", "- Downloading UProC repository")
        print("- Downloading UProC repository ✓")

        print("- Extracting repository", end="", flush=True)
        tar = tarfile.open(os.path.join(tmpdir, "uproc.tar.gz"))
        tar.extractall(tmpdir)
        print("\r- Extracting repository ✓")

        print("- Running configure", end="", flush=True)
        configure = subprocess.Popen(["./configure", "--prefix", install_dir],
                                     cwd=os.path.join(tmpdir, "uproc-1.2.0"),
                                     stdout=output, stderr=output)
        if configure.wait() != 0:
            print("\n\nError while running configure. Rerun the command with '--verbose' for subprocess output.")
            sys.exit(1)

        print("\r- Running configure ✓")

        print("- Running make", end="", flush=True)
        make = subprocess.Popen("make", cwd=os.path.join(tmpdir, "uproc-1.2.0"),
                                stdout=output, stderr=output)
        if make.wait() != 0:
            print("\n\nError while running make. Rerun the command with '--verbose' for subprocess output.")
            sys.exit(1)

        print("\r- Running make ✓")

        print("- Running make install", end="", flush=True)
        install = subprocess.Popen(["make", "install"], cwd=os.path.join(tmpdir, "uproc-1.2.0"),
                                   stdout=output, stderr=output)
        if install.wait() != 0:
            print("\n\nError while running make install. Rerun the command with '--verbose' for subprocess output.")
            sys.exit(1)

        print("\r- Running make install ✓\n")


def download_uproc_win(url: str, install_dir: str) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        download(url, tmpdir, "uproc.zip", "- Downloading UProC")
        print("- Downloading UProC ✓")

        print("- Extracting archive", end="", flush=True)
        with zipfile.ZipFile(os.path.join(tmpdir, "uproc.zip"), 'r') as zip_ref:
            zip_ref.extractall(tmpdir)
        shutil.move(os.path.join(tmpdir, "uproc-1.2.0-win-x86_64"), install_dir)
        print("\r- Extracting archive ✓\n")
