#!/opt/venv/bin/python
"""Entry point of the `dmrdetector` command in the image.

It is installed as /opt/venv/bin/dmrdetector. The package is
/opt/dmrdetector/src/dmrdetector; the venv finds it through dmrdetector.pth.
"""
import sys

from dmrdetector.cli import main


if __name__ == "__main__":
    sys.exit(main())
