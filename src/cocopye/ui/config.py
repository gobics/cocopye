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
import sys

import argparse
from importlib import resources
from typing import Tuple

import numba
from appdirs import user_config_dir, user_cache_dir, user_log_dir
from tomlkit import parse, dumps, TOMLDocument


ARGS: argparse.Namespace
"""Description"""
CONFIG_FILE: str
"""Description"""
CONFIG: TOMLDocument
"""Description"""


def init() -> None:
    global ARGS, CONFIG_FILE, CONFIG
    CONFIG_FILE, CONFIG = parse_config()
    ARGS = parse_args()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="cocopye",
        description="Feature-based prediction of genome quality indices",
        epilog="See https://github.com/gobics/cocopye/wiki for more detailed usage instructions."
    )

    parser.add_argument("--pfam24", action='store_true',
                        help="Use Pfam database version 24 (instead of 28)")

    parser.add_argument("--offline", action='store_true',
                        help="Run CoCoPyE in offline mode")

    subparsers = parser.add_subparsers(title="subcommands", dest="subcommand")

    # Subparser run

    run_parser = subparsers.add_parser(
        "run",
        help="Calculate contamination and completeness",
        description="Calculate contamination and completeness"
    )

    run_parser.add_argument("-i", "--infolder", required=True,
                            help="Input folder containing bins in FASTA format")
    run_parser.add_argument("-o", "--outfile", default="cocopye_output.csv",
                            help="Output file (default: cocopye_output.csv)")
    run_parser.add_argument("--file-extensions", default="fasta,fna,fa",
                            help="Allowed file extensions for the FASTA files (default: fasta,fna,fa)")
    run_parser.add_argument("-t", "--threads", default=str(min(8, numba.config.NUMBA_NUM_THREADS)),
                            help="Number of threads")
    run_parser.add_argument("-v", "--verbosity", default="standard",
                            help="Output verbosity (standard, extended, full; default: standard)")
    run_parser.add_argument("--pfam24", action='store_true', dest="pfam24_run",
                            help="Use Pfam database version 24 (instead of 28)")

    # Subparser database

    if CONFIG["advanced"]["enable_db_creator"]:
        db_parser = subparsers.add_parser(
            "database",
            help="Create a new database matrix",
            description="Create a new database matrix"
        )

        db_parser.add_argument("-i", "--infolder", required=True)
        db_parser.add_argument("-m", "--metadata", required=True)
        db_parser.add_argument("-o", "--outfolder", required=True)
        db_parser.add_argument("-t", "--threads", default=str(min(8, numba.config.NUMBA_NUM_THREADS)),
                               help="Number of threads")

    # Subparser toolbox

    setup_parser = subparsers.add_parser(
        "setup",
        help="Tools to setup or repair CoCoPyE",
        description="Tools to setuo or repair CoCoPyE"
    )

    subparsers_toolbox = setup_parser.add_subparsers(title="subcommands", dest="subcommand_setup")

    subparsers_toolbox.add_parser("update-database", help="Update CoCoPyE database")
    subparsers_toolbox.add_parser("cleanup", help="Remove automatically downloaded/generated files")

    testrun_subparser = subparsers_toolbox.add_parser("testrun",
                                                      help="Start a testrun to ensure everything works as expected")
    testrun_subparser.add_argument("--pfam24", action='store_true', dest="pfam24_testrun",
                                   help="Use Pfam database version 24 (instead of 28)")

    dl_subparser = subparsers_toolbox.add_parser("download-dependencies", help="Download missing dependencies")
    dl_subparser.add_argument("-v", "--verbose", action='store_true', help="Show output of subprocesses")
    dl_subparser.add_argument("--pfam24", action='store_true', dest="pfam24_dl",
                              help="Use Pfam database version 24 (instead of 28)")

    # Subparser web

    if CONFIG["advanced"]["enable_webserver"]:
        subparsers.add_parser(
            "web",
            help="Start webserver",
            description="Start webserver"
        )

    return _merge_flag(parser.parse_args(), "pfam24")


def _merge_flag(parsed, flag):
    parsed_dict = vars(parsed)
    keys = [key for key in parsed_dict.keys() if flag in key]

    setattr(parsed, flag, any([parsed_dict[key] for key in keys]))

    for key in keys:
        if key == flag:
            continue
        delattr(parsed, key)

    return parsed


def parse_config() -> Tuple[str, TOMLDocument]:
    # Check several locations for an existing config file
    locations = [
        os.environ.get("COCOPYE_CONFIG", ""),
        "cocopye.toml",
        os.path.join(user_config_dir("cocopye"), "cocopye.toml")
    ]

    for config_file in locations:
        try:
            with open(config_file) as config:
                return config_file, parse(config.read())
                pass
        except IOError:
            pass

    # If we cannot find any, we create a new default configuration file
    print("Unable to find config file.")

    default_config = parse(resources.read_text("cocopye.ui", "config.toml"))
    default_config["server"]["tmpdir"] = os.path.join(user_cache_dir("cocopye"), "server")
    default_config["server"]["logdir"] = os.path.join(user_log_dir("cocopye"), "server")

    os.makedirs(user_config_dir("cocopye"), exist_ok=True)
    f = open(os.path.join(user_config_dir("cocopye"), "cocopye.toml"), "w")
    f.write(dumps(default_config))
    f.close()

    print("Created a new default configuration file at ", os.path.join(user_config_dir("cocopye"), "cocopye.toml.\n"))

    # And then we read it
    try:
        with open(os.path.join(user_config_dir("cocopye"), "cocopye.toml")) as config:
            return os.path.join(user_config_dir("cocopye"), "cocopye.toml"), parse(config.read())
    except IOError:
        print("I wasn't able to read the config file I just created. This shouldn't happen.")
        sys.exit(1)


def change_config(table: str, elem: str, new_value: str) -> None:
    global CONFIG
    CONFIG[table][elem] = new_value

    f = open(CONFIG_FILE, "w")
    f.write(dumps(CONFIG))
    f.close()
