#!/usr/bin/env bash
# Rebuild the analysis environment and restore the branch tip.
#
# The venv lives OUTSIDE the repository root (/home/user/.venv-t) and is
# therefore not part of any persisted snapshot. The branch ref can also be
# reset to the base commit between sessions. Both have happened, so this
# script exists to make recovery one command instead of exploration.
#
#   bash scripts/setup_environment.sh
#
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV="/home/user/.venv-t"
BRANCH="arena/7b6ab982-telomerase-allosteric-switch"

echo "== repo: $REPO"
cd "$REPO"

# ---- 1. restore the branch tip if the local ref was reset underneath us
if git rev-parse --verify --quiet "origin/$BRANCH" >/dev/null 2>&1; then
    if ! git merge-base --is-ancestor "origin/$BRANCH" HEAD 2>/dev/null; then
        echo "== local HEAD is behind origin/$BRANCH - restoring"
        git fetch -q origin "$BRANCH"
        git reset --hard "origin/$BRANCH"
    fi
fi
echo "== HEAD: $(git log --oneline -1)"

# ---- 2. build the venv if it is missing or broken
if [ ! -x "$VENV/bin/python" ] || ! "$VENV/bin/python" -c "import rdkit" 2>/dev/null; then
    echo "== building $VENV"
    python3 -m venv "$VENV"
    "$VENV/bin/pip" install -q --upgrade pip
    "$VENV/bin/pip" install -q \
        rdkit biopython numpy scipy gemmi pandas matplotlib pytest \
        requests pyyaml sympy
fi
echo "== python: $("$VENV/bin/python" -V), rdkit $("$VENV/bin/python" -c 'import rdkit;print(rdkit.__version__)')"

# ---- 3. confirm the tree is coherent
echo "== running the three ledgers' self-tests"
"$VENV/bin/python" scripts/solve/wp38_archive.py --self-test | tail -1
"$VENV/bin/python" scripts/solve/wp46_closure.py --self-test | tail -1
"$VENV/bin/python" scripts/solve/wp48_new_designs.py --self-test | tail -1

echo "== ready"
