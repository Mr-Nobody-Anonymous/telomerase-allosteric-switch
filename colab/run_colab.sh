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
echo "=== self-test 1/3: geometry against the bundled CIFs (offline) ==="
"$PY" colab/tpp1_colab.py --selftest
SELF=$?
if [ "$SELF" -ne 0 ]; then
  echo
  echo "SELF-TEST 1 FAILED (exit $SELF). Not continuing to the download stage." >&2
  exit "$SELF"
fi

echo
echo "=== self-test 2/3: verdict rule against real coordinate subsets (offline) ==="
"$PY" colab/verdict.py --selftest
SELF2=$?
if [ "$SELF2" -ne 0 ]; then
  echo
  echo "SELF-TEST 2 FAILED (exit $SELF2). The verdict rule does not reproduce the" >&2
  echo "documented numbers. Not continuing." >&2
  exit "$SELF2"
fi

echo
echo "=== self-test 3/3: does AGENT_START_HERE.md match the code? ==="
"$PY" colab/test_colab.py
SELF3=$?
if [ "$SELF3" -ne 0 ]; then
  echo
  echo "SELF-TEST 3 FAILED (exit $SELF3). The document and the code disagree." >&2
  exit "$SELF3"
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
