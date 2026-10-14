"""`python -m dmrdetector` is the same as the `dmrdetector` command."""

import sys

from dmrdetector import cli


sys.exit(cli.main())
