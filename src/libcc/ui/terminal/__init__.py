from ..config import parse_args, parse_config
from ..external import check_and_download_dependencies


def main():
    """
    Entry point of the terminal user interface. This is called by `src/cctui.py`.
    """
    print("Hello from the terminal user interface.\n")

    args = parse_args()
    config_file, configuration = parse_config(args.config)

    check_and_download_dependencies(configuration, config_file)

    print("At this point the execution of the tool's main part (the calculation of completeness and contamination)")
    print("would start. As this is not yet finished, there is currently nothing to see here.")
