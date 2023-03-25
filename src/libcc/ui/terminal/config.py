import os
import sys

import argparse
from appdirs import user_config_dir
from tomlkit import parse


def parse_args():
    parser = argparse.ArgumentParser(
        prog="libcc",
        description="A description",
        epilog="Text at the bottom of help"
    )

    parser.add_argument('-c', '--config')

    return parser.parse_args()


def parse_config(config_file: str = None):
    # If there is an explicit path given, we only check this one and exit if the file
    # cannot be found in the given location
    if config_file is not None:
        try:
            with open(config_file) as config:
                return config_file, parse(config.read())
        except IOError:
            print("Configuration file specified in parameter not found. Exiting.")
            sys.exit(1)

    # Otherwise we check several locations for an existing config file
    for config_file in "libcc.toml", os.path.join(user_config_dir("libcc"), "libcc.toml"):
        try:
            with open(config_file) as config:
                return config_file, parse(config.read())
                pass
        except IOError:
            pass

    # If we cannot find any, we create a new default configuration file
    print("Unable to find config file.")
    os.makedirs(user_config_dir("libcc"), exist_ok=True)
    f = open(os.path.join(user_config_dir("libcc"), "libcc.toml"), "w")
    f.write(DEFAULT_CONFIG)
    f.close()
    print("Created a new default configuration file at ", os.path.join(user_config_dir("libcc"), "libcc.toml.\n"))

    # And then we read it
    try:
        with open(os.path.join(user_config_dir("libcc"), "libcc.toml")) as config:
            return os.path.join(user_config_dir("libcc"), "libcc.toml"), parse(config.read())
    except IOError:
        print("I wasn't able to read the config file I just created. This shouldn't happen.")
        sys.exit(1)


DEFAULT_CONFIG = """[external]
prodigal_bin = "prodigal"
uproc_bin = "uproc-prot"
uproc_db = "none"
uproc_models = "none"
libcc_db = "none"

#####################################################################################################
# During normal use of the tool, it shouldn't be necessary to change the variables below this line. #
#####################################################################################################

[download]
prodigal_url_windows = "https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.windows.exe"
prodigal_url_linux = "https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux"
uproc_src = "http://uproc.gobics.de/downloads/uproc/uproc-1.2.0.tar.gz"
uproc_win = "http://uproc.gobics.de/downloads/uproc/uproc-1.2.0-win-x86_64.zip"
"""


