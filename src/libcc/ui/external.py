import gzip
import os
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from tomlkit import dumps
import requests
from tqdm import tqdm

from appdirs import user_data_dir, user_cache_dir


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

    pfam_status = check_pfam_db(config["external"]["uproc_db"])
    pfam_missing = pfam_status == "not found"
    pfam_error = pfam_status == "error"
    pfam_output = (cross, _red(pfam_status)) if pfam_missing or pfam_error \
        else (tick, _green(pfam_status))

    model_status = check_model(config["external"]["uproc_models"])
    model_missing = model_status == "not found"
    model_error = model_status == "error"
    model_output = (cross, _red(model_status)) if model_missing or model_error \
        else (tick, _green(model_status))

    print("External dependencies:")
    print(prodigal_output[0] + " Prodigal\t\t" + prodigal_output[1])
    print(uproc_output[0] + " UProC\t\t\t" + uproc_output[1])
    print("  " + pfam_output[0] + " Pfam database\t" + pfam_output[1])
    print("  " + model_output[0] + " Models\t\t" + model_output[1])
    print(cross + " libcc database\t" + _red("todo"))
    print()

    if prodigal_error or uproc_error or pfam_error or model_error:
        # This is not a good solution, but probably sufficient for now, as I assume that this shouldn't really happen
        # during normal use anyway
        print("Error: Some dependencies were found, but they don't seem to work. Please check them manually.\n")
        sys.exit(1)

    if prodigal_missing or uproc_missing or pfam_missing or model_missing:
        print("Some dependencies are missing.")
        print("If you have already downlaoded them, you can specify their path in the configuration file.")
        print("Apart from that, libcc can try to download them automatically.")
        print()
        download_now = input("Download missing dependencies now? [Y/n] ")
        print()
        if download_now.strip() == "N" or download_now.strip() == "n":
            print("Missing dependencies were not downloaded. Exiting.\n")
            sys.exit(1)

        opsys = platform.system()
        if opsys != "Windows" and opsys != "Linux":
            print("No supported operating system detected. Exiting.\n")
            sys.exit(1)

        if prodigal_missing:
            url = config["download"]["prodigal_url_linux"] if opsys == "Linux"\
                else config["download"]["prodigal_url_windows"]

            download(url, os.path.join(user_data_dir("libcc"), "prodigal"), "prodigal", "- Downloading Prodigal")

            config["external"]["prodigal_bin"] = os.path.join(user_data_dir("libcc"), "prodigal", "prodigal")
            f = open(config_file, "w")
            f.write(dumps(config))
            f.close()

            print("- Downloading Prodigal ✓\n")

        if uproc_missing:
            if opsys == "Linux":
                build_uproc_prot(config["download"]["uproc_src"], os.path.join(user_data_dir("libcc"), "uproc"))

                config["external"]["uproc_bin"] = os.path.join(user_data_dir("libcc"), "uproc", "bin", "uproc-prot")
                config["external"]["uproc_import_bin"] = os.path.join(user_data_dir("libcc"), "uproc", "bin", "uproc-import")
                f = open(config_file, "w")
                f.write(dumps(config))
                f.close()
            elif opsys == "Windows":
                print("TODO")  # TODO

        if pfam_missing:
            # not using /tmp, because of the large file size
            with tempfile.TemporaryDirectory(prefix="libcc_", dir=user_cache_dir(None)) as tmpdir:
                download(
                    config["download"]["pfam_db"],
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
                uproc_import = subprocess.Popen([config["external"]["uproc_import_bin"],
                                                os.path.join(tmpdir, "pfam.uprocdb"),
                                                os.path.join(user_data_dir("libcc"), "pfam_db")],
                                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                uproc_import.wait()
                print("\r- Importing database ✓                             \n")

                config["external"]["uproc_db"] = os.path.join(user_data_dir("libcc"), "pfam_db")
                f = open(config_file, "w")
                f.write(dumps(config))
                f.close()

        if model_missing:
            with tempfile.TemporaryDirectory() as tmpdir:
                download(
                    config["download"]["model"],
                    tmpdir,
                    "model.tar.gz",
                    "- Downloading UProC Model"
                )
                print("- Downloading UProc Model ✓")

                print("- Extracting UProC model")
                tar = tarfile.open(os.path.join(tmpdir, "model.tar.gz"))
                tar.extractall(user_data_dir("libcc"))
                print("\r- Extracting UProC model ✓\n")

                config["external"]["uproc_models"] = os.path.join(user_data_dir("libcc"), "model")
                f = open(config_file, "w")
                f.write(dumps(config))
                f.close()

        print("Downloads finished. Please restart the application.\n")
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


def check_pfam_db(pfam_dir: str):
    try:
        content = os.listdir(pfam_dir)
    except FileNotFoundError:
        return "not found"

    files = ["fwd.ecurve", "idmap", "prot_thresh_e2", "prot_thresh_e3", "rev.ecurve"]
    complete = all([file in content for file in files])

    return "found" if complete else "error"


def check_model(model_dir: str):
    try:
        content = os.listdir(model_dir)
    except FileNotFoundError:
        return "not found"

    files = ["aa_probs", "alphabet", "codon_scores", "orf_thresh_e1", "orf_thresh_e2", "substmat"]
    complete = all([file in content for file in files])

    return "found" if complete else "error"


def download(url: str, dirname: str, fname: str, label: str, chunk_size=1024):
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
        for data in resp.iter_content(chunk_size=chunk_size):
            size = file.write(data)
            bar.update(size)

    # This is probably only necessary for Prodigal, but I assume it won't hurt in other cases
    os.chmod(os.path.join(dirname, fname), stat.S_IRUSR ^ stat.S_IWUSR ^ stat.S_IXUSR)


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
