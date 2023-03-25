from .config import parse_args, parse_config
from ..external import check_dependencies


def main():
    print("Hello from the terminal user interface.\n")

    args = parse_args()
    config_file, configuration = parse_config(args.config)

    check_dependencies(configuration, config_file)
