import os
import shutil
import sys
import importlib.util

import numpy as np
import pandas as pd
import pkg_resources

from appdirs import user_data_dir, user_config_dir
from numba import set_num_threads

from .. import config
from ..external import check_and_download_dependencies
from ..external.data import update_cocopye_db
from ...matrices import DatabaseMatrix, load_u8mat_from_file, QueryMatrix
from ...pfam import count_pfams
from ... import constants


def main() -> None:
    """
    Entry point of the terminal user interface. This is called by `src/cli.py`.
    """
    try:
        print("Welcome to CoCoPyE v" + pkg_resources.get_distribution('CoCoPyE').version + ".\n")
    except pkg_resources.DistributionNotFound:
        print("Welcome to CoCoPyE.\n")

    config.init()
    if config.ARGS.subcommand in ["run", "database"]:
        set_num_threads(int(config.ARGS.threads))

    if config.ARGS.subcommand == "cleanup":
        cleanup()
        sys.exit(0)

    if config.ARGS.subcommand == "toolbox":
        if config.ARGS.update_database:
            update_cocopye_db(constants.COCOPYE_DB, config.CONFIG["external"]["cocopye_db"])
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
    count_mat, seq_list, _ = count_pfams(
        config.CONFIG["external"]["uproc_orf_bin"],
        config.CONFIG["external"]["uproc_prot_bin"],
        os.path.join(config.CONFIG["external"]["uproc_pfam_db"], "24" if config.ARGS.pfam24 else "28"),
        config.CONFIG["external"]["uproc_models"],
        config.ARGS.infolder,
        num_threads=config.ARGS.threads
    )

    # Read and sort metadata
    metadata = pd.read_csv(config.ARGS.metadata, sep=",")
    metadata = metadata.set_index("sequence").loc[seq_list].reset_index()

    # Create universal markers, one set for each superkingdom
    superkingdoms = np.unique(metadata["superkingdom"].to_numpy())
    all_markers = {}
    for superkingdom in superkingdoms:
        inds = metadata.index[metadata["superkingdom"] == superkingdom].to_list()
        submatrix = count_mat[inds]
        universal_markers = DatabaseMatrix(submatrix).universal_markers(threshold=0.95)
        all_markers[superkingdom] = universal_markers

    # Save everything to files
    os.makedirs(config.ARGS.outfolder)
    for superkingdom in all_markers:
        np.save(os.path.join(config.ARGS.outfolder, "universal_" + superkingdom + ".npy"), all_markers[superkingdom])
    np.savez_compressed(os.path.join(config.ARGS.outfolder, "count_matrix.npz"), count_mat)
    metadata.to_csv(os.path.join(config.ARGS.outfolder, "metadata.csv"), sep=",", index=False)


def run():
    pfam_version = "24" if config.ARGS.pfam24 else "28"

    db_mat = DatabaseMatrix(
        load_u8mat_from_file(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "count_matrix.npz")),
        pd.read_csv(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "metadata.csv"), sep=",")
    )
    query_mat, bin_ids, count_ratio = count_pfams(
        config.CONFIG["external"]["uproc_orf_bin"],
        config.CONFIG["external"]["uproc_prot_bin"],
        os.path.join(config.CONFIG["external"]["uproc_pfam_db"], pfam_version),
        config.CONFIG["external"]["uproc_models"],
        config.ARGS.infolder,
        config.ARGS.file_extension,
        config.ARGS.threads
    )
    query_mat = QueryMatrix(query_mat).with_database(db_mat, constants.K)

    assert len(bin_ids) == query_mat.mat().shape[0]

    universal_arc = np.load(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "universal_Archaea.npy"))
    universal_bac = np.load(os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "universal_Bacteria.npy"))

    preestimates_arc = query_mat.preestimates(universal_arc)
    preestimates_bac = query_mat.preestimates(universal_bac)

    estimates = query_mat.estimates()

    feature_mat_comp = query_mat.into_feature_mat(estimates, constants.RESOLUTION_COMP)
    feature_mat_cont = query_mat.into_feature_mat(estimates, constants.RESOLUTION_CONT)

    ml_estimates_comp = feature_mat_comp.ml_estimates(
        os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "model_comp.pickle")).clip(0, 1)
    ml_estimates_cont = feature_mat_cont.ml_estimates(
        os.path.join(config.CONFIG["external"]["cocopye_db"], pfam_version, "model_cont.pickle")).clip(0, 1000000)

    taxonomy = query_mat.taxonomy()
    knn_scores = query_mat.knn_scores()

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
        "2_contamination," +
        "count_length_ratio," +
        "knn_score," +
        "taxonomy\n"
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
            str(ml_estimates_cont[idx]) + "," +
            str(count_ratio[idx]) + "," +
            str(knn_scores[idx]) + "," +
            taxonomy[idx] + "\n"
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
