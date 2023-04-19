import os
import subprocess
import tarfile
import tempfile

from ..external import download, _green, _red, _TICK, _CROSS


def check_prodigal(prodigal_bin: str):
    try:
        process = subprocess.Popen(
            [prodigal_bin, '-v'],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )
    except FileNotFoundError:
        return 1, "prodigal", _CROSS + " Prodigal\t\t" + _red("not found")

    version = "v" + process.stderr.read().strip().rpartition("V")[2].rpartition(":")[0]

    return 0, "prodigal", _TICK + " Prodigal\t\t" + _green(version)


def check_uproc(uproc_bin: str):
    try:
        process = subprocess.Popen(
            [uproc_bin, '-v'],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
        )
    except FileNotFoundError:
        return 1, "uproc", _CROSS + " UProC\t\t\t" + _red("not found")

    version = "v" + process.stdout.read().strip().partition("\n")[0].rpartition("version ")[2]

    return 0, "uproc", _TICK + " UProC\t\t\t" + _green(version)


def build_uproc_prot(url: str, install_dir: str):
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
                                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # TODO: hardcoded
        configure.wait()
        print("\r- Running configure ✓")

        print("- Running make", end="", flush=True)
        make = subprocess.Popen("make", cwd=os.path.join(tmpdir, "uproc-1.2.0"),
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # TODO: hardcoded
        make.wait()
        print("\r- Running make ✓")

        print("- Running make install", end="", flush=True)
        install = subprocess.Popen(["make", "install"], cwd=os.path.join(tmpdir, "uproc-1.2.0"),
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # TODO: hardcoded
        install.wait()
        print("\r- Running make install ✓\n")
