#!/usr/bin/env bash
# One-command entry point for a Colab (or any Linux) agent.
#
#   bash colab/run_colab.sh            # self-test, then the full download + analysis
#   bash colab/run_colab.sh --selftest # offline verification only
#
# Exits non-zero if the self-test fails, so an agent can gate on it.
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1

# --- pick an interpreter that can see numpy/scipy/gemmi -------------------------
pick_python() {
  for cand in "${PYTHON:-}" .venv/bin/python python3 python; do
    [ -z "$cand" ] && continue
    if command -v "$cand" >/dev/null 2>&1 || [ -x "$cand" ]; then
      if "$cand" -c "import numpy, scipy, gemmi" >/dev/null 2>&1; then
        echo "$cand"; return 0
      fi
    fi
  done
  for cand in "${PYTHON:-}" python3 python; do
    [ -z "$cand" ] && continue
    if command -v "$cand" >/dev/null 2>&1; then
      echo "  attempting to install deps into $cand" >&2
      "$cand" -m pip -q install --upgrade numpy scipy gemmi >/dev/null 2>&1 \
        || "$cand" -m pip -q install --upgrade --break-system-packages numpy scipy gemmi >/dev/null 2>&1
      if "$cand" -c "import numpy, scipy, gemmi" >/dev/null 2>&1; then
        echo "$cand"; return 0
      fi
    fi
  done
  return 1
}

PY=$(pick_python)
if [ -z "${PY:-}" ]; then
  echo "ERROR: no python interpreter with numpy/scipy/gemmi could be prepared." >&2
  exit 2
fi

echo "=== environment ==="
echo "interpreter: $PY"
"$PY" -V
"$PY" -c "import numpy, scipy, gemmi; print('numpy', numpy.__version__, '| scipy', scipy.__version__, '| gemmi', gemmi.__version__)"

echo
echo "=== self-test (offline) ==="
"$PY" colab/tpp1_colab.py --selftest
SELF=$?
if [ "$SELF" -ne 0 ]; then
  echo
  echo "SELF-TEST FAILED (exit $SELF). Not continuing to the download stage." >&2
  exit "$SELF"
fi

if [ "${1:-}" = "--selftest" ]; then
  echo
  echo "self-test only requested; stopping."
  exit 0
fi

echo
echo "=== full run (needs internet) ==="
"$PY" colab/tpp1_colab.py
RUN=$?

echo
echo "=== report ==="
if [ -f colab_out/report.md ]; then
  cat colab_out/report.md
else
  echo "no report produced"
fi

exit "$RUN"
