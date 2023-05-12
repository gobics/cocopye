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

import requests
from tqdm import tqdm
from appdirs import user_data_dir

from ..config import change_config


def check_and_download_dependencies(config, config_file):
    """
    This function checks if all neccessary dependencies are available. If there's something missing, it will offer the
    user to download the missing files automatically. The function writes some status messages to stdout and waits for
    user input on stdin when asking if it should download the missing parts. Furthermore, it may terminate the program
    if there are corrupted dependencies or errors during download. If this function returns successfully, it is
    (hopefully)  guaranteed that all required files are present and working.

    :param config: Dictionary generated by `config.parse_config`.
    :param config_file: Location of the configuration file
    """
    ok, missing, error = check_dependencies(config)

    if len(error) > 0:
        # It might be useful to provide automatic solutions for common errors (at least try to redownload the file or
        # provide a more useful error message).
        print("\nError: Some dependencies were found, but they don't seem to work. Please check them manually.\n")
        sys.exit(1)

    if len(missing) > 0:
        print("\nSome dependencies are missing.")
        print("If you have already downlaoded them, you can specify their path in the configuration file.")
        print("Apart from that, libcc can try to download them automatically.")
        print()
        download_now = input("Download missing dependencies now? [Y/n] ")
        print()
        if download_now.strip() == "N" or download_now.strip() == "n":
            print("Missing dependencies were not downloaded. Exiting.\n")
            sys.exit(1)

        download_dependencies(config, config_file, missing)

        print("Downloads finished. Rechecking dependencies...\n")
        ok, missing, error = check_dependencies(config)
        if len(missing) > 0 or len(error) > 0:
            print("There are still errors or missing files. This should not happen and you have likely encountered a")
            print("bug in the application. Exiting.")
            sys.exit(1)

        print("Check successful. Continung execution of the application.\n")


def check_dependencies(config):
    """
    This function checks the status of dependencies. It writes an overview of the result to stdout. This function is
    mainly for use by `check_and_download_dependencies`.

    :param config: Dictionary generated by `config.parse_config`.
    :return: 3-tuple (ok, missing, error) where each element is a list of dependencies that have the respective status
    Example: `(["prodigal", "uproc"], ["pfam"], ["model"])` means that Prodigal and UProC are present and working, the
    pfam database is missing and the model is there but corrupted.
    """
    from .bin import check_prodigal, check_uproc
    from .data import check_pfam_db, check_model

    checks = [
        check_prodigal(config["external"]["prodigal_bin"]),
        check_uproc(config["external"]["uproc_bin"]),
        check_pfam_db(config["external"]["uproc_db"]),
        check_model(config["external"]["uproc_models"])
    ]

    status = [[], [], []]  # ok, missing, error

    print("External dependencies:")
    for lst, identifier, output in checks:
        print(output)
        status[lst].append(identifier)
    print()

    return status[0], status[1], status[2]


def download_dependencies(config, config_file, missing):
    """
    This tries to download and install missing dependencies. During the process it prints a bunch of status and progress
    messages to stdout. If the detected operating system is not Windows or Linux it will terminate the application.
    This function is mainly for use by `check_and_download_dependencies`.

    :param config: Dictionary generated by `config.parse_config`.
    :param config_file: Configuration file location
    :param missing: A list of missing dependencies as generated by `check_dependencies`
    """
    from .bin import build_uproc_prot, download_uproc_win
    from .data import download_pfam_db, download_model

    opsys = platform.system()
    if opsys != "Windows" and opsys != "Linux":
        print("No supported operating system detected. Exiting.\n")
        sys.exit(1)

    if "prodigal" in missing:
        url = config["download"]["prodigal_url_linux"] if opsys == "Linux" \
            else config["download"]["prodigal_url_windows"]

        download(url, os.path.join(user_data_dir("libcc"), "prodigal"), "prodigal", "- Downloading Prodigal")

        config = change_config(
            config, config_file, "external", "prodigal_bin",
            os.path.join(user_data_dir("libcc"), "prodigal", "prodigal")
        )

        print("- Downloading Prodigal ✓\n")

    if "uproc" in missing:
        if opsys == "Linux":
            build_uproc_prot(config["download"]["uproc_src"], os.path.join(user_data_dir("libcc"), "uproc"))

            config = change_config(
                config, config_file, "external", "uproc_bin",
                os.path.join(user_data_dir("libcc"), "uproc", "bin", "uproc-prot")
            )

            config = change_config(
                config, config_file, "external", "uproc_import_bin",
                os.path.join(user_data_dir("libcc"), "uproc", "bin", "uproc-import")
            )

            config = change_config(
                config, config_file, "external", "uproc_orf_bin",
                os.path.join(user_data_dir("libcc"), "uproc", "bin", "uproc-orf")
            )
        elif opsys == "Windows":
            download_uproc_win(config["download"]["uproc_win"], os.path.join(user_data_dir("libcc"), "uproc"))

            config = change_config(
                config, config_file, "external", "uproc_bin",
                os.path.join(user_data_dir("libcc"), "uproc", "uproc-prot.exe")
            )

            config = change_config(
                config, config_file, "external", "uproc_import_bin",
                os.path.join(user_data_dir("libcc"), "uproc", "uproc-import.exe")
            )

            config = change_config(
                config, config_file, "external", "uproc_orf_bin",
                os.path.join(user_data_dir("libcc"), "uproc", "uproc-orf.exe")
            )

    if "pfam" in missing:
        download_pfam_db(config["download"]["pfam_db"], config["external"]["uproc_import_bin"])

        config = change_config(
            config, config_file, "external", "uproc_db",
            os.path.join(user_data_dir("libcc"), "pfam_db")
        )

    if "model" in missing:
        download_model(config["download"]["model"])

        change_config(
            config, config_file, "external", "uproc_models",
            os.path.join(user_data_dir("libcc"), "model")
        )


def download(url: str, dirname: str, fname: str, label: str, chunk_size=1024):
    """
    An auxiliary function to download a file. It shows a progress bar which gets removed once the download is complete.

    :param url: Download URL
    :param dirname: Destination directory (without filename)
    :param fname: Destination file (without directory path)
    :param label: Progress bar label
    :param chunk_size: Download chunk size
    """
    os.makedirs(dirname, exist_ok=True)
    # https://gist.github.com/yanqd0/c13ed29e29432e3cf3e7c38467f42f51
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


def _red(s: str) -> str:
    return "\033[91m" + s + "\033[0m"


_TICK = _green("✓")
_CROSS = _red("✗")
