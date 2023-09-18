import os
import shutil
import sys
import importlib.util

import numpy as np
import pandas as pd
import pkg_resources

from appdirs import user_data_dir, user_config_dir
# from numba import set_num_threads

from .. import config
from ..external import check_and_download_dependencies
from ...matrices import DatabaseMatrix, load_u8mat_from_file
from ...pfam import count_pfams
from ... import pfam, constants


def main() -> None:
    """
    Entry point of the terminal user interface. This is called by `src/cli.py`.
    """
    try:
        print("Welcome to CoCoPyE v" + pkg_resources.get_distribution('CoCoPyE').version + ".\n")
    except pkg_resources.DistributionNotFound:
        print("Welcome to CoCoPyE.\n")

    config.init()

    if config.ARGS.subcommand == "cleanup":
        cleanup()
        sys.exit(0)

    check_and_download_dependencies()

    if config.ARGS.subcommand == "database":
        create_database()

    if config.ARGS.subcommand == "run":
        if not os.path.isdir(config.ARGS.infolder):
            print("Error: The specified input folder does not exist. Exiting.")
            sys.exit(1)

        run()

    if config.ARGS.subcommand == "web":
        web()

    if config.ARGS.subcommand is None:
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
    if config.ARGS.filter is not None:
        filter_file = open(config.ARGS.filter, "r")
        filter_list = [seq.strip() for seq in filter_file.read().split("\n")]
        filter_file.close()

    db_mat = pfam.create_database_matrix(
        config.CONFIG["external"]["uproc_orf_bin"],
        config.CONFIG["external"]["uproc_prot_bin"],
        os.path.join(config.CONFIG["external"]["uproc_pfam_db"], "24" if config.ARGS.pfam24 else "28"),
        config.CONFIG["external"]["uproc_models"],
        config.ARGS.infile,
        filter_list
    )

    db_mat.save_to_file(config.ARGS.outfile)


def run():
    pfam_version = "24" if config.ARGS.pfam24 else "28"

    db_mat = DatabaseMatrix(load_u8mat_from_file(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "mat_pfam.npy")))
    #metadata = pd.read_csv(os.path.join(config.CONFIG["external"]["cocopye_db"], "metadata.csv"))
    query_mat, bin_ids = count_pfams(
        config.CONFIG["external"]["uproc_orf_bin"],
        config.CONFIG["external"]["uproc_prot_bin"],
        os.path.join(config.CONFIG["external"]["uproc_pfam_db"], pfam_version),
        config.CONFIG["external"]["uproc_models"],
        config.ARGS.infolder,
        config.ARGS.file_extension
    )
    query_mat = query_mat.with_database(db_mat, constants.K)

    assert len(bin_ids) == query_mat.mat().shape[0]

    universal_arc = np.load(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "universal_arc.npy"))
    universal_bac = np.load(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "universal_bac.npy"))

    preestimates_arc = query_mat.preestimates(universal_arc)
    preestimates_bac = query_mat.preestimates(universal_bac)

    estimates = query_mat.estimates()

    feature_mat_comp = query_mat.into_feature_mat(estimates, constants.RESOLUTION)
    feature_mat_cont = query_mat.into_feature_mat(estimates, constants.RESOLUTION)
    # TODO difference 24 and 28; maybe use two completely different database folders
    ml_estimates_comp = feature_mat_comp.ml_estimates(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "model_comp.pickle"))
    ml_estimates_cont = feature_mat_cont.ml_estimates(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "model_cont.pickle"))

    #taxonomy = query_mat.taxonomy(metadata)

    outfile = open(config.ARGS.outfile, "w")
    outfile.write(
        "bin," +
        "0_completeness_arc," +
        "0_contamination_arc," +
        "0_completeness_bac," +
        "0_contamination_bac," +
        "1_completeness," +
        "1_contamination," +
        "1_num_markers," +
        "2_completeness," +
        "2_contamination\n"
        #"2_contamination," +
        #"taxonomy\n"
    )
    for idx in range(len(bin_ids)):
        outfile.write(
            bin_ids[idx] + "," +
            str(preestimates_arc[idx, 0]) + "," +
            str(preestimates_arc[idx, 1]) + "," +
            str(preestimates_bac[idx, 0]) + "," +
            str(preestimates_bac[idx, 1]) + "," +
            str(estimates[idx, 0]) + "," +
            str(estimates[idx, 1]) + "," +
            str(estimates[idx, 2]) + "," +
            str(ml_estimates_comp[idx]) + "," +
            str(ml_estimates_cont[idx]) + "\n"
            #str(ml_estimates_cont[idx]) + "," +
            #taxonomy[idx] + "\n"
        )
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
