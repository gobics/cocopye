import os
import sys

from cocopye.ui import terminal


def main() -> None:
    try:
        terminal.main()
    except KeyboardInterrupt:
        print("\n\nInterrupted.")
        try:
            sys.exit(130)
        except SystemExit:
            os._exit(130)


if __name__ == "__main__":
    main()
