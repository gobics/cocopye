import os
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from urllib.request import urlopen
from tomlkit import dumps

from appdirs import user_data_dir


def check_dependencies(config, config_file):
    tick = _green("✓")
    cross = _red("✗")

    prodigal_status = check_prodigal(config["external"]["prodigal_bin"])
    prodigal_missing = prodigal_status == "not found"
    prodigal_error = prodigal_status == "error"
    prodigal_output = (cross, _red(prodigal_status)) if prodigal_missing or prodigal_error \
        else (tick, _green(prodigal_status))

    uproc_status = check_uproc(config["external"]["uproc_bin"])
    uproc_missing = uproc_status == "not found"
    uproc_error = uproc_status == "error"
    uproc_output = (cross, _red(uproc_status)) if uproc_missing or uproc_error \
        else (tick, _green(uproc_status))

    print("External dependencies:")
    print(prodigal_output[0] + " Prodigal\t" + prodigal_output[1])
    print(uproc_output[0] + " UProC Prot\t" + uproc_output[1])
    print()

    if prodigal_error or uproc_error:
        # This is not a good solution, but probably sufficient for now, as I assume that this shouldn't really happen
        # during normal use anyway
        print("Some dependencies were found, but they don't seem to work. Please check them manually.")
        sys.exit(1)

    if prodigal_missing or uproc_missing:
        print("Some dependencies are missing.")
        print("If you have already downlaoded them, you can specify their path in the configuration file.")
        print("Apart from that, libcc can try to download them automatically.")
        print()
        download_now = input("Download missing dependencies now? [Y/n] ")

        if download_now.strip() == "N" or download_now.strip() == "n":
            print("Missing dependencies were not downloaded. Exiting.")
            sys.exit(1)

        opsys = platform.system()
        if opsys != "Windows" and opsys != "Linux":
            print("No supported operating system detected. Exiting.")
            sys.exit(1)

        if prodigal_missing:
            url = config["download"]["prodigal_url_linux"] if opsys == "Linux"\
                else config["download"]["prodigal_url_windows"]

            print("Downloading Prodigal", end="")
            sys.stdout.flush()

            download_file(
                url,
                os.path.join(user_data_dir("libcc"), "prodigal"),
                "prodigal"
            )

            config["external"]["prodigal_bin"] = os.path.join(user_data_dir("libcc"), "prodigal", "prodigal")
            f = open(config_file, "w")
            f.write(dumps(config))
            f.close()

            print("\rProdigal for " + opsys + " downloaded")

        if uproc_missing:
            if opsys == "Linux":
                url = config["download"]["uproc_src"]

                print("Downloading and compiling UProC", end="")
                sys.stdout.flush()

                build_uproc_prot(url, os.path.join(user_data_dir("libcc"), "uproc"))

                config["external"]["uproc_bin"] = os.path.join(user_data_dir("libcc"), "uproc", "bin", "uproc-prot")
                f = open(config_file, "w")
                f.write(dumps(config))
                f.close()

                print("\rUProC for " + opsys + " downloaded and compiled")
            elif opsys == "Windows":
                print("TODO")  # TODO

        print("Please restart the application.")
        sys.exit(0)


def _green(s: str) -> str:
    return "\033[92m" + s + "\033[0m"


def _red(s: str) -> str:
    return "\033[91m" + s + "\033[0m"


def check_prodigal(prodigal_bin: str):
    try:
        process = subprocess.Popen(
            [prodigal_bin, '-v'],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )
    except FileNotFoundError:
        return "not found"
    #except:
    #    return "error"

    version = "v" + process.stderr.read().strip().rpartition("V")[2].rpartition(":")[0]

    return version


def check_uproc(uproc_bin: str):
    try:
        process = subprocess.Popen(
            [uproc_bin, '-v'],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
        )
    except FileNotFoundError:
        return "not found"
    #except:
    #    return "error"

    version = "v" + process.stdout.read().strip().partition("\n")[0].rpartition("version ")[2]

    return version


def download_file(url: str, path: str, filename: str):
    os.makedirs(path, exist_ok=True)
    with urlopen(url) as f:
        content = f.read()
    with open(os.path.join(path, filename), "wb") as target:
        target.write(content)
    os.chmod(os.path.join(path, filename), stat.S_IRUSR ^ stat.S_IWUSR ^ stat.S_IXUSR)


def build_uproc_prot(url: str, install_dir: str):
    with tempfile.TemporaryDirectory() as tmpdir:
        download_file(url, tmpdir, "uproc.tar.gz")
        tar = tarfile.open(os.path.join(tmpdir, "uproc.tar.gz"))
        tar.extractall(tmpdir)
        configure = subprocess.Popen(["./configure", "--prefix", install_dir],
                                     cwd=os.path.join(tmpdir, "uproc-1.2.0"),
                                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # TODO: hardcoded
        configure.wait()
        make = subprocess.Popen("make", cwd=os.path.join(tmpdir, "uproc-1.2.0"),
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # TODO: hardcoded
        make.wait()
        install = subprocess.Popen(["make", "install"], cwd=os.path.join(tmpdir, "uproc-1.2.0"),
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # TODO: hardcoded
        install.wait()
