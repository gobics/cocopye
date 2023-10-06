import os
import sys

import argparse
from importlib import resources
from typing import Tuple

from appdirs import user_config_dir, user_cache_dir
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
        description="A description",
        epilog="Text at the bottom of help"
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
    run_parser.add_argument("-o", "--outfile", required=True, help="Output file")
    run_parser.add_argument("--file-extension", default="fna",
                            help="File extension of the bin FASTA files (default: fna)")
    run_parser.add_argument("-t", "--threads", default="8", help="Number of threads")
    run_parser.add_argument("-v", "--verbosity", default="standard",
                            help="Output verbosity (standard, extended, everything; default: standard)")

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
        db_parser.add_argument("-t", "--threads", default="8", help="Number of threads")

    # Subparser toolbox

    toolbox_parser = subparsers.add_parser(
        "toolbox",
        help="Tools to update or repair CoCoPyE (TODO: better description)",
        description="Tools to update or repair CoCoPyE (TODO: better description)"
    )

    toolbox_parser.add_argument("--update-database", action='store_true', help="Update CoCoPyE database")

    # Other subparsers

    subparsers.add_parser(
        "cleanup",
        help="Remove automatically downloaded/generated files",
        description="Remove automatically downloaded/generated files"
    )

    if CONFIG["advanced"]["enable_webserver"]:
        subparsers.add_parser(
            "web",
            help="Start webserver",
            description="Start webserver"
        )

    return parser.parse_args()


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
