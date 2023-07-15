import os
import shutil
import sys
import importlib.util
from typing import List

import pkg_resources

from appdirs import user_data_dir, user_config_dir
# from numba import set_num_threads

from ..config import ARGS, CONFIG
from ..external import check_and_download_dependencies
from ...matrices import DatabaseMatrix, load_u8mat_from_file, QueryMatrix
from ...preprocessing.kmer import count_kmers
from ...preprocessing.pfam import count_pfams
from ...preprocessing import kmer, pfam


def main() -> None:
    """
    Entry point of the terminal user interface. This is called by `src/cli.py`.
    """
    try:
        print("Welcome to CoCoPyE v" + pkg_resources.get_distribution('CoCoPyE').version + ".\n")
    except pkg_resources.DistributionNotFound:
        print("Welcome to CoCoPyE.\n")

    if ARGS.subcommand == "cleanup":
        cleanup()
        sys.exit(0)

    check_and_download_dependencies()

    if ARGS.subcommand == "database":
        create_database()

    if ARGS.subcommand == "run":
        if not os.path.isdir(ARGS.infolder):
            print("Error: The specified input folder does not exist. Exiting.")
            sys.exit(1)

        if ARGS.kmer:
            run_kmer()
        else:
            run_pfam()

    if ARGS.subcommand == "web":
        web()

    if ARGS.subcommand is None:
        print("No subcommand selected. If you don't know which subcommand to use, you will most likely")
        print("want to use \"run\". Execute \"cocopye -h\" or \"cocopye run -h\" for usage instructions.")
        sys.exit(1)


def web():
    # This is not optimal, but currently there seems to be no way to check if an extra was selected during
    # package installation
    for dep in ["fastapi", "celery", "uvicorn"]:
        if importlib.util.find_spec(dep) is None:
            print("Please run 'pip install cocopye[web] for webserver support.")
            sys.exit(1)

    from ..web.server import run_server
    run_server()


def create_database() -> None:
    filter_list = None
    if ARGS.filter is not None:
        filter_file = open(ARGS.filter, "r")
        filter_list = [seq.strip() for seq in filter_file.read().split("\n")]
        filter_file.close()

    if ARGS.kmer:
        db_mat = kmer.create_database_matrix(CONFIG["external"]["prodigal_bin"], ARGS.infile, sequences=filter_list)
    else:
        db_mat = pfam.create_database_matrix(
            CONFIG["external"]["uproc_orf_bin"],
            CONFIG["external"]["uproc_bin"],
            CONFIG["external"]["uproc_db"],
            CONFIG["external"]["uproc_models"],
            ARGS.infile,
            filter_list
        )

    db_mat.save_to_file(ARGS.outfile)


def run_pfam() -> None:
    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(CONFIG["external"]["cocopye_db"], "mat1234.npy")))
    query_mat, bin_ids = count_pfams(
        CONFIG["external"]["uproc_orf_bin"],
        CONFIG["external"]["uproc_bin"],
        CONFIG["external"]["uproc_db"],
        CONFIG["external"]["uproc_models"],
        ARGS.infolder,
        ARGS.file_extension
    )

    run(db_mat, query_mat, bin_ids, ARGS.outfile, k=ARGS.k)


def run_kmer():
    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(CONFIG["external"]["cocopye_db"], "kmerdb1234.npy")))
    query_mat, bin_ids = count_kmers(CONFIG["external"]["prodigal_bin"], ARGS.infolder, ARGS.file_extension)

    run(db_mat, query_mat, bin_ids, ARGS.outfile, 0.15, ARGS.k)  # TODO: knn mit range


def run(
    db_mat: DatabaseMatrix,
    query_mat: QueryMatrix,
    bin_ids: List[str],
    outfile_path: str,
    var_thresh: float = None,
    k: int = 30
):
    estimates = query_mat.estimates(db_mat, k, var_thresh=var_thresh)  # TODO: Remove hardcoded k

    assert len(bin_ids) == query_mat.mat().shape[0]

    outfile = open(outfile_path, "w")
    for idx in range(len(bin_ids)):
        outfile.write(
            bin_ids[idx] + "," +
            str(estimates[idx, 0]) + "," +
            str(estimates[idx, 1]) + "," +
            str(estimates[idx, 2]) + "\n")
    outfile.close()


def cleanup() -> None:
    print("This will remove the following directory with all its contents:")
    print("\t" + user_data_dir("cocopye") + "\n")
    proceed = input("Proceed? [y/N] ")
    if proceed == "y" or proceed == "Y":
        shutil.rmtree(user_data_dir("cocopye"))
        print(user_data_dir("cocopye") + " removed.")
    else:
        print("Aborted.")
        sys.exit(0)
    proceed = input("\nDo you also want to remove the user configuration at " + user_config_dir("cocopye") + "? [y/N] ")
    if proceed == "y" or proceed == "Y":
        shutil.rmtree(user_config_dir("cocopye"))
        print(user_config_dir("cocopye") + " removed.")
