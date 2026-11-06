"""
Root pytest configuration.

WHY THIS FILE EXISTS
--------------------
`pytest.ini` sets `pythonpath = src .`, which lets *pytest itself* import
`telotrade` from `src/`. It does not help **child processes**.

Several tests shell out to the CLI like this:

    subprocess.run([sys.executable, "-m", "telotrade.validate"], ...)

That spawns a fresh interpreter with a fresh `sys.path`. The `pythonpath` ini
option has no effect on it, and the package is not installed into the
environment, so the child dies with:

    ModuleNotFoundError: No module named 'telotrade'

which made `test_cli_execution_module` and `test_cli_fails_on_corrupt_data`
fail with exit code 1 for *both* the passing and the corrupt-data case. The
corrupt-data test was worse than broken -- it asserted `returncode == 1` and
was passing for the wrong reason, since the CLI never ran at all.

Setting PYTHONPATH here fixes every subprocess call in one place, rather than
patching each test to add its own `env=`. It also keeps the tests honest: the
CLI genuinely executes, so their exit-code assertions mean what they claim.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"

# Make the package importable in child processes spawned by the CLI tests.
_existing = os.environ.get("PYTHONPATH", "")
_parts = [str(SRC), str(ROOT)]
for p in _parts:
    if p not in _existing.split(os.pathsep):
        _existing = f"{_existing}{os.pathsep}{p}" if _existing else p
os.environ["PYTHONPATH"] = _existing

# Belt and braces for this interpreter, in case pytest was started without
# honouring pytest.ini's pythonpath (e.g. `pytest -c /dev/null`).
for p in (str(SRC), str(ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)
