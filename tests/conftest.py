"""pytest setup: puts src/ and the test helpers on sys.path.

The project has no pytest.ini, setup.cfg or tox.ini (pyproject.toml holds no
pytest configuration), so the tests import the checkout directly and need no
installation.
"""

from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
for path in (HERE.parent / "src", HERE / "unit", HERE / "unit" / "oracles"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
