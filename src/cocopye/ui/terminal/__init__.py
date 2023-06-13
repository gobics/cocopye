import argparse
import os
import shutil
import sys
from typing import List

import pkg_resources

from appdirs import user_data_dir, user_config_dir
from numba import set_num_threads
from tomlkit import TOMLDocument

from ..config import parse_args, parse_config
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

    args = parse_args()
    config_file, configuration = parse_config(args.config)

    if args.subcommand == "cleanup":
        cleanup()
        sys.exit(0)

    check_and_download_dependencies(configuration, config_file)

    if args.subcommand == "database":
        create_database(args, configuration, kmer_flag=args.kmer)

    if args.subcommand == "run":
        if not os.path.isdir(args.infolder):
            print("Error: The specified input folder does not exist. Exiting.")
            sys.exit(1)

        if args.kmer:
            run_kmer(args, configuration)
        else:
            run_pfam(args, configuration)

    if args.subcommand is None:
        print("No subcommand selected. If you don't know which subcommand to use, you will most likely")
        print("want to use \"run\". Execute \"cocopye -h\" or \"cocopye run -h\" for usage instructions.")
        sys.exit(1)


def create_database(args: argparse.Namespace, config: TOMLDocument, kmer_flag: bool = False) -> None:
    filter_list = None
    if args.filter is not None:
        filter_file = open(args.filter, "r")
        filter_list = [seq.strip() for seq in filter_file.read().split("\n")]
        filter_file.close()

    if kmer_flag:
        db_mat = kmer.create_database_matrix(args.infile, sequences=filter_list)
    else:
        db_mat = pfam.create_database_matrix(
            config["external"]["uproc_orf_bin"],
            config["external"]["uproc_bin"],
            config["external"]["uproc_db"],
            config["external"]["uproc_models"],
            args.infile,
            filter_list
        )

    db_mat.save_to_file(args.outfile)


def run_pfam(args: argparse.Namespace, config: TOMLDocument) -> None:
    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(config["external"]["cocopye_db"], "mat1234.npy")))
    query_mat, bin_ids = count_pfams(
        config["external"]["uproc_orf_bin"],
        config["external"]["uproc_bin"],
        config["external"]["uproc_db"],
        config["external"]["uproc_models"],
        args.infolder,
        args.file_extension
    )

    run(db_mat, query_mat, bin_ids, args.outfile)


def run_kmer(args: argparse.Namespace, config: TOMLDocument):
    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(config["external"]["cocopye_db"], "kmer_mat1234.npy")))
    query_mat, bin_ids = count_kmers(args.infolder, args.file_extension)

    run(db_mat, query_mat, bin_ids, args.outfile, 0.3)  # TODO: knn mit range


def run(db_mat: DatabaseMatrix, query_mat: QueryMatrix, bin_ids: List[str], outfile_path: str, var_thresh: float = None):
    estimates = query_mat.estimates(db_mat, 30, var_thresh=var_thresh)  # TODO: Remove hardcoded k

    assert len(bin_ids) == query_mat.mat().shape[0]

    outfile = open(outfile_path, "w")
    for idx in range(len(bin_ids)):
        outfile.write(bin_ids[idx] + "," + str(estimates[idx, 0]) + "," + str(estimates[idx, 1]) + "\n")
    outfile.close()


def cleanup() -> None:
    print("This will remove the following directory with all its contents:")
    print("\t" + user_data_dir("cocopye") + "\n")
    proceed = input("Proceed? [y/N]")
    if proceed == "y" or proceed == "Y":
        shutil.rmtree(user_data_dir("cocopye"))
        print(user_data_dir("cocopye") + "removed.")
    else:
        print("Aborted.")
        sys.exit(0)
    proceed = input("\nDo you also want to remove the user configuration at " + user_config_dir("cocopye") + "? [y/N]")
    if proceed == "y" or proceed == "Y":
        shutil.rmtree(user_config_dir("cocopye"))
        print(user_config_dir("cocopye") + "removed.")
