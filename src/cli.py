import sys

from cocopye.ui import terminal


def main() -> None:
    terminal.main()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted.")
        sys.exit(130)
