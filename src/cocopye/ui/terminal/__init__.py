import argparse

from tomlkit import TOMLDocument

from ..config import parse_args, parse_config
from ..external import check_and_download_dependencies
from ...matrices import DatabaseMatrix, load_u8mat_from_file
from ...preprocessing.pfam import create_database_matrix, count_pfams


def main() -> None:
    """
    Entry point of the terminal user interface. This is called by `src/cli.py`.
    """
    # print("Hello from the terminal user interface.\n")

    args = parse_args()
    config_file, configuration = parse_config(args.config)

    check_and_download_dependencies(configuration, config_file)

    if args.subcommand == "database":
        create_database(args, configuration)

    if args.subcommand == "run":
        run_pfam(args, configuration)

    if args.subcommand is None:
        print("No subcommand selected. Exiting.")


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
    db_mat = DatabaseMatrix(load_u8mat_from_file(config["external"]["pfam_db_mat"]))  # TODO: Add config key (external module)
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
