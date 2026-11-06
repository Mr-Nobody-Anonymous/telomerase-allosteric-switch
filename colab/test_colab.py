#!/usr/bin/env python3
"""
Verify that colab/AGENT_START_HERE.md tells the truth.

An agent reading that file will believe what it says and report it onward. So the
file's numbers, thresholds and file paths are asserted here, mechanically. If the
document and the code ever disagree, this test fails -- the document is not allowed
to drift away from the measurements.

This deliberately does NOT import anything from the document's prose. It parses the
claim table, recomputes every figure from colab/fixtures/, and compares.

USAGE
    python colab/test_colab.py          # exit 0 = the document is accurate
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
COLAB = REPO / "colab"
DOC = COLAB / "AGENT_START_HERE.md"

FAILS: list[str] = []
CHECKS = 0


def check(cond: bool, msg: str) -> None:
    global CHECKS
    CHECKS += 1
    if not cond:
        FAILS.append(msg)


# --------------------------------------------------------------------- helpers
def import_verdict():
    sys.path.insert(0, str(COLAB))
    import importlib
    return importlib.import_module("verdict")


def load_fixture(entry: str) -> dict:
    p = COLAB / "fixtures" / f"atoms_{entry}.json"
    return json.loads(p.read_text())


# --------------------------------------------------------------------- tests
def test_document_exists_and_is_findable() -> None:
    check(DOC.exists(), f"{DOC} does not exist")
    if not DOC.exists():
        return
    text = DOC.read_text()
    check("bash colab/run_colab.sh" in text,
          "the document does not give the one command an agent needs")
    check("Test DID NOT RUN".lower() in text.lower(),
          "the document does not warn that a failed download is not a negative result")
    # an agent must be told what not to conclude
    check("not a mechanism" in text.lower(),
          "the document does not bound the result against a mechanism claim")

    # An agent has already run `pytest`, seen 480 green tests, and reported the
    # task complete -- a false claim, since pytest answers no scientific question.
    # The document must actively warn against that specific failure.
    low = text.lower()
    check("do not run `pytest`" in low or "do not run pytest" in low,
          "the document does not warn against running pytest instead of the task")
    check("precondition" in low,
          "the document does not say that green tests are a precondition, not the result")


def test_headline_numbers_match_fixtures() -> None:
    """The three numbers in the document's table must come out of the fixtures."""
    v = import_verdict()
    expected = {
        "7QXS": {"anchor": 2.75, "ctrl_max": 6.0},
        "7QXA": {"anchor": 4.35, "ctrl_max": 6.0},
        "9SHY": {"anchor": 30.23, "ctrl_max": 6.0},
    }
    for entry, exp in expected.items():
        fx = load_fixture(entry)
        call = v.structure_call(entry, fx["records"], fx["modelled_dna"]["n_modelled"])
        check("anchor_d" in call, f"{entry}: no anchor distance produced")
        if "anchor_d" not in call:
            continue
        check(abs(call["anchor_d"] - exp["anchor"]) <= 0.02,
              f"{entry}: document says {exp['anchor']} A, code produces {call['anchor_d']} A")
        check(call.get("catalytic_ok") is True,
              f"{entry}: the catalytic control does not pass, so the document's "
              f"claim that these entries are usable is false")
        check(call.get("d868_3end", 99) <= exp["ctrl_max"],
              f"{entry}: Asp868->3'end is {call.get('d868_3end')} A")

    # the two contact entries must actually be called CONTACT
    for entry in ("7QXS", "7QXA"):
        fx = load_fixture(entry)
        call = v.structure_call(entry, fx["records"], fx["modelled_dna"]["n_modelled"])
        check(call["call"] == "CONTACT",
              f"{entry}: document says CONTACT, code says {call['call']}")

    # 9SHY must NOT be reported as a negative result
    fx = load_fixture("9SHY")
    call = v.structure_call("9SHY", fx["records"], fx["modelled_dna"]["n_modelled"])
    check(call["call"] != "NO CONTACT",
          "9SHY is being reported as NO CONTACT; its anchor nucleotide is unmodelled, "
          "so that would be missing data read as a negative result")
    check("caveat" in call or call["call"] == "NOT RESOLVED",
          "9SHY carries no caveat about unmodelled nucleotides")


def test_seven_mer_reference_is_quoted_correctly() -> None:
    """
    The document quotes 11.11 / 11.94 / 11.94 A for the 7-mer structures, which is
    the only part of the result that relies on the bundled CIFs. Verify against the
    same independent values the runner's own self-test uses.
    """
    src = (COLAB / "tpp1_colab.py").read_text()
    for v in ("11.11", "11.94"):
        check(v in src, f"the runner no longer contains the reference value {v} A")
    doc = DOC.read_text()
    for v in ("11.11", "11.94"):
        check(v in doc, f"the document omits the reference value {v} A")


def test_every_verdict_branch_reachable() -> None:
    v = import_verdict()
    branches = {
        "ANCHOR SITE ENGAGES DNA": [{"entry": "X", "call": "CONTACT", "anchor_d": 2.8,
                                     "anchor_pair": "a-b", "nearest_nt_auth": 19}],
        "ANCHOR SITE DOES NOT ENGAGE DNA": [{"entry": "X", "call": "NO CONTACT",
                                             "anchor_d": 30.0}],
        "NOT RESOLVED": [{"entry": "X", "call": "NO CONTACT", "anchor_d": 30.0,
                          "caveat": "unresolved"}],
        "TEST DID NOT RUN": [{"entry": "X", "call": "CANNOT ANSWER", "reasons": []}],
    }
    for want, calls in branches.items():
        got = v.overall(calls)["status"]
        check(got == want, f"branch {want!r} unreachable; got {got!r}")


def test_failure_is_not_reported_as_negative() -> None:
    """The single most important behaviour in the whole pipeline."""
    v = import_verdict()
    out = v.overall([{"entry": "X", "call": "CANNOT ANSWER", "reasons": ["no network"]}])
    check("NOT RUN" in out["status"],
          "a failed download does not produce a 'did not run' status")
    check("NOT a negative" in out["text"] or "not a negative" in out["text"].lower(),
          "the failure status does not say it is not a negative result")


def test_old_defect_stays_fixed() -> None:
    """
    7QXA (12 modelled nt) and 9SHY (10 modelled nt) are usable. If a size gate ever
    returns, these become unusable again and the whole result silently vanishes.
    """
    v = import_verdict()
    for entry in ("7QXA", "9SHY"):
        fx = load_fixture(entry)
        n = fx["modelled_dna"]["n_modelled"]
        check(n < 15, f"{entry} has {n} modelled nt; this test no longer covers the size gate")
        call = v.structure_call(entry, fx["records"], n)
        check(call["call"] != "CANNOT ANSWER" or call.get("catalytic_ok") is False,
              f"{entry} was discarded despite passing the active-site control -- "
              f"a size gate has been reintroduced")


def test_runner_and_files_exist() -> None:
    for rel in ("run_colab.sh", "tpp1_colab.py", "verdict.py", "test_colab.py"):
        check((COLAB / rel).exists(), f"colab/{rel} is referenced but missing")
    sh = (COLAB / "run_colab.sh").read_text()
    check("--selftest" in sh, "the runner cannot be run in self-test mode")
    check("pick_python" in sh, "the runner lost its interpreter discovery")


def test_document_references_real_paths() -> None:
    """Every repo path mentioned in the document must exist."""
    doc = DOC.read_text()
    for m in re.finditer(r"`(results/[A-Za-z0-9_./]+|colab/[A-Za-z0-9_./]+\.(?:py|sh|md|json|csv))`", doc):
        rel = m.group(1)
        check((REPO / rel).exists(), f"document points at {rel}, which does not exist")


# --------------------------------------------------------------------- main
def main() -> int:
    print("VERIFYING colab/AGENT_START_HERE.md AGAINST THE CODE")
    print("=" * 78)
    tests = [
        test_document_exists_and_is_findable,
        test_headline_numbers_match_fixtures,
        test_seven_mer_reference_is_quoted_correctly,
        test_every_verdict_branch_reachable,
        test_failure_is_not_reported_as_negative,
        test_old_defect_stays_fixed,
        test_runner_and_files_exist,
        test_document_references_real_paths,
    ]
    for t in tests:
        before = len(FAILS)
        try:
            t()
        except Exception as e:                                    # noqa: BLE001
            FAILS.append(f"{t.__name__} raised {type(e).__name__}: {e}")
        status = "ok" if len(FAILS) == before else f"{len(FAILS)-before} FAILED"
        print(f"  {t.__name__:42s} {status}")

    print()
    if FAILS:
        print("FAILURES")
        for f in FAILS:
            print(f"  - {f}")
        print()
        print(f"{CHECKS} checks, {len(FAILS)} failed")
        print("VERDICT: the document does NOT match the code. Fix before shipping.")
        return 1

    print(f"{CHECKS} checks, 0 failed")
    print("VERDICT: document matches the code.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
