import argparse
import os
import shutil
import sys
import pkg_resources

from appdirs import user_data_dir, user_config_dir
from tomlkit import TOMLDocument

from ..config import parse_args, parse_config
from ..external import check_and_download_dependencies
from ...matrices import DatabaseMatrix, load_u8mat_from_file
from ...preprocessing.pfam import create_database_matrix, count_pfams


def main() -> None:
    """
    Entry point of the terminal user interface. This is called by `src/cli.py`.
    """
    print("Welcome to CoCoPyE v" + pkg_resources.get_distribution('CoCoPyE').version + ".\n")

    args = parse_args()
    config_file, configuration = parse_config(args.config)

    if args.subcommand == "cleanup":
        cleanup()
        sys.exit(0)

    check_and_download_dependencies(configuration, config_file)

    if args.subcommand == "database":
        create_database(args, configuration)

    if args.subcommand == "run":
        run_pfam(args, configuration)

    if args.subcommand is None:
        print("No subcommand selected. If you don't know which subcommand to use, you will most likely")
        print("want to use \"run\". Execute \"cocopye -h\" or \"cocopye run -h\" for usage instructions.")
        sys.exit(1)


def create_database(args: argparse.Namespace, config: TOMLDocument) -> None:
    filter_list = None
    if args.filter is not None:
        filter_file = open(args.filter, "r")
        filter_list = [seq.strip() for seq in filter_file.read().split("\n")]
        filter_file.close()

    db_mat = create_database_matrix(
        config["external"]["uproc_orf_bin"],
        config["external"]["uproc_bin"],
        config["external"]["uproc_db"],
        config["external"]["uproc_models"],
        args.infile,
        filter_list
    )

    db_mat.save_to_file(args.outfile)


def run_pfam(args: argparse.Namespace, config: TOMLDocument) -> None:
    if not os.path.isdir(args.infolder):
        print("Error: The specified input folder does not exist. Exiting.")
        sys.exit(1)

    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(config["external"]["cocopye_db"], "mat1234.npy")))
    query_mat, bin_ids = count_pfams(
        config["external"]["uproc_orf_bin"],
        config["external"]["uproc_bin"],
        config["external"]["uproc_db"],
        config["external"]["uproc_models"],
        args.infolder,
        args.file_extension
    )

    estimates = query_mat.estimates(db_mat, 30)  # TODO: Remove hardcoded k

    assert len(bin_ids) == query_mat.mat().shape[0]

    outfile = open(args.outfile, "w")
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
