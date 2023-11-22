# Copyright 2023 Niklas Birth
#
# This file is part of CoCoPyE.
#
# CoCoPyE is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CoCoPyE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CoCoPyE. If not, see <https://www.gnu.org/licenses/>.

import os
import shutil
import subprocess
import sys
import importlib.util
import tempfile

import numpy as np
import pandas as pd
import pkg_resources

from appdirs import user_data_dir, user_config_dir
from numba import set_num_threads

from .. import config
from ..external import check_and_download_dependencies
from ..external.data import update_cocopye_db
from ...core import log
from ...matrices import DatabaseMatrix
from ...pfam import count_pfams
from ... import constants, core


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

    if config.ARGS.subcommand == "toolbox":
        if config.ARGS.subcommand_toolbox == "update-database":
            update_cocopye_db(constants.COCOPYE_DB, config.CONFIG["external"]["cocopye_db"])
            sys.exit(0)
        if config.ARGS.subcommand_toolbox == "cleanup":
            cleanup()
            sys.exit(0)
        if config.ARGS.subcommand_toolbox == "testrun":
            with tempfile.TemporaryDirectory() as tmpdir:
                command = (["cocopye"]
                           + (["--offline"] if config.ARGS.offline else [])
                           + (["--pfam24"] if config.ARGS.pfam24 else [])
                           + ["run", "-i", os.path.join(config.CONFIG["external"]["cocopye_db"], "testdata"),
                              "-o", os.path.join(tmpdir, "tempfile.csv")])

                print("Starting testrun.\n\033[37m(" + " ".join(command) + ")\033[0m\n")
                print("=======================================================================\n")
                cocopye_process = subprocess.Popen(command)
                cocopye_process.wait()
                print("\n=======================================================================")

                if cocopye_process.returncode == 0:
                    print("\n\033[92mTestrun sucessful.\033[0m")
                else:
                    print("\n\033[91mTestrun failed.\033[0m")

            sys.exit(0)
        if config.ARGS.subcommand_toolbox == "download-dependencies":
            check_and_download_dependencies()
            sys.exit(0)
        else:
            print("Unrecognized subcommand. Exiting.")
            sys.exit(1)

    check_and_download_dependencies(check_only=True)

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
    results = core.core(config.CONFIG["external"]["cocopye_db"],
                        config.CONFIG["external"]["uproc_orf_bin"],
                        config.CONFIG["external"]["uproc_prot_bin"],
                        config.CONFIG["external"]["uproc_pfam_db"],
                        config.CONFIG["external"]["uproc_models"],
                        config.ARGS.infolder,
                        24 if config.ARGS.pfam24 else 28,
                        config.ARGS.file_extension,
                        config.ARGS.threads
                        )

    log("Saving results to file")
    outfile = open(config.ARGS.outfile, "w")

    if config.ARGS.verbosity == "full":
        outfile.write("bin,stage,method,1_completeness_arc,1_contamination_arc,1_completeness_bac,1_contamination_bac,"
                      "2_completeness,2_contamination,2_num_markers,3_completeness,3_contamination,"
                      "coding_density,knn_score,taxonomy,taxonomy_level,notes\n")
    elif config.ARGS.verbosity == "extended":
        outfile.write("bin,completeness,contamination,stage,method,num_markers,coding_density,knn_score,"
                      "taxonomy,taxonomy_level,notes\n")
    else:
        outfile.write("bin,completeness,contamination,method,taxonomy,taxonomy_level,notes\n")

    for result in results:
        outfile.write(result.to_csv(config.ARGS.verbosity) + "\n")

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
