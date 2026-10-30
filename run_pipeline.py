#!/usr/bin/env python3
"""
Single-Entry Automated Pipeline Runner for Human Telomerase TERT–TPP1 Study.

Executes all pipeline stages:
1. Validation and SHA-256 integrity verification (7TRE, 7TRD, 7TRF, 7TRC, 2I46).
2. Sequence-aware residue correspondence and construct offset reconciliation.
3. Quantitative interface geometry profiling (contacts, BSA, polar, TEL-patch).
4. Cross-state comparative structural analysis and robustness testing.
5. Standardized results CSV table generation in results/tables/.
6. Publication-grade figure generation in results/figures/.

Usage:
    python run_pipeline.py
    python run_pipeline.py --raw-dir data/raw --tables-dir results/tables --figures-dir results/figures
"""

import sys
from pathlib import Path

# Ensure src/ is on PYTHONPATH for direct execution
_REPO_ROOT = Path(__file__).resolve().parent
_SRC_DIR = _REPO_ROOT / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from telotrade.pipeline import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
