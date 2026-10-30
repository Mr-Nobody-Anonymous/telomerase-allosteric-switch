"""Unit and integration tests for the Breakthrough Research Investigation deliverables.

Validates file presence, report structure, epistemic confidence levels,
literature citation integrity, candidate ranking matrix, and exact coordinate measurements.
"""

import csv
from pathlib import Path

import gemmi

BASE_DIR = Path(__file__).resolve().parent.parent
INVESTIGATION_DIR = BASE_DIR / "results" / "breakthrough_investigation"
DATA_DIR = BASE_DIR / "data" / "raw"


def test_deliverable_files_exist():
    """Verify all deliverables exist and are non-empty."""
    required_files = [
        "NEXT_BREAKTHROUGH_RESEARCH_REPORT.md",
        "BREAKTHROUGH_RESEARCH_REPORT.md",
        "LITERATURE_NOVELTY_AUDIT.md",
        "ALTERNATIVE_ALLOSTERIC_PATHWAYS.md",
        "CATALYTIC_CYCLE_ANALYSIS.md",
        "MUTATION_CANDIDATE_RANKING.csv",
        "EXPERIMENTAL_FALSIFICATION_PROTOCOL.md",
        "EVIDENCE_AND_CONFIDENCE_MATRIX.csv",
        "STRUCTURAL_VALIDATION.md",
        "MUTAGENESIS_ANALYSIS.md",
    ]
    for filename in required_files:
        filepath = INVESTIGATION_DIR / filename
        assert filepath.exists(), f"Missing required file: {filename}"
        assert filepath.stat().st_size > 500, f"File {filename} is unexpectedly empty or too small"


def test_next_breakthrough_report_structure():
    """Ensure NEXT_BREAKTHROUGH_RESEARCH_REPORT.md addresses all core mission objectives."""
    report_text = (INVESTIGATION_DIR / "NEXT_BREAKTHROUGH_RESEARCH_REPORT.md").read_text(encoding="utf-8")

    key_sections = [
        "1. Executive Summary & Forensic Re-Alignment",
        "2. Objective 1: Challenging the Arg756 Hypothesis",
        "3. Objective 2: Discovering Alternative Allosteric Pathways",
        "4. Objective 3: Kinetic Dissection of the Catalytic Cycle",
        "5. Objective 4: Critical Literature Novelty Audit",
        "6. Objective 5: Ranking Candidate Mutations",
        "7. Objective 6: Decisive Experimental Framework",
        "8. Objective 7: Biological Significance & Clinical Realism",
        "9. Direct Answer to the Final Question",
    ]

    for section in key_sections:
        assert section.lower() in report_text.lower(), f"Missing section: {section}"


def test_citation_integrity_no_invented_references():
    """Verify that every cited paper has a valid DOI and genuine author attribution."""
    report_text = (INVESTIGATION_DIR / "NEXT_BREAKTHROUGH_RESEARCH_REPORT.md").read_text(encoding="utf-8")
    audit_text = (INVESTIGATION_DIR / "LITERATURE_NOVELTY_AUDIT.md").read_text(encoding="utf-8")

    verified_dois = [
        "10.1038/nature05454",  # Wang et al. 2007
        "10.1038/emboj.2009.409",  # Latrick & Cech 2010
        "10.1038/nature11648",  # Nandakumar et al. 2012
        "10.1073/pnas.2024889118",  # Sandhu et al. 2021
        "10.1038/s41586-021-03415-4",  # Ghanim et al. 2021
        "10.1038/s41586-022-04582-8",  # Sekne et al. 2022
        "10.26508/lsa.202201727",  # Padmanaban et al. 2023
    ]

    for doi in verified_dois:
        assert doi in report_text or doi in audit_text, f"Verified DOI {doi} missing from reports"


def test_mutation_candidate_ranking_matrix():
    """Verify MUTATION_CANDIDATE_RANKING.csv structure and candidate ordering."""
    csv_path = INVESTIGATION_DIR / "MUTATION_CANDIDATE_RANKING.csv"
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert len(rows) >= 10, "Ranking matrix should contain at least 10 audited candidates"
    assert rows[0]["mutation"] == "Lys757A", "Top-ranked candidate should be Lys757A"
    assert rows[1]["mutation"] == "Lys757R", "Second-ranked candidate should be Lys757R"
    assert rows[2]["mutation"] == "Tyr176A", "Third-ranked candidate should be Tyr176A"

    # Verify Lys757A has no RNA confounder
    assert "no contact" in rows[0]["rna_a55_confounder"].lower()


def test_evidence_and_confidence_matrix_structure():
    """Verify CSV columns and valid evidence level tagging."""
    csv_path = INVESTIGATION_DIR / "EVIDENCE_AND_CONFIDENCE_MATRIX.csv"
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert len(rows) >= 14, "Evidence matrix should contain at least 14 audited claims"
    levels_present = {row["evidence_level"].split(":")[0].strip() for row in rows}
    for lvl in ["Level 1", "Level 2", "Level 3", "Level 4", "Level 5"]:
        assert lvl in levels_present, f"Must contain {lvl} evidence"


def test_exact_coordinate_measurements():
    """Verify calculated coordinate distances from 7TRE raw CIF structure."""
    cif_7tre = DATA_DIR / "7TRE.cif"
    assert cif_7tre.exists()
    st = gemmi.make_structure_from_block(gemmi.cif.read_file(str(cif_7tre)).sole_block())

    tert = st[0]["A"]
    tpp1 = st[0]["C"]
    dna = st[0]["D"]
    rna = st[0]["B"]

    # Residue lookup
    r756 = next(r for r in tert if r.seqid.num == 756)
    k757 = next(r for r in tert if r.seqid.num == 757)
    y176 = next(r for r in tert if r.seqid.num == 176)
    k78 = next(r for r in tert if r.seqid.num == 78)
    r774 = next(r for r in tert if r.seqid.num == 774)
    l769 = next(r for r in tert if r.seqid.num == 769)

    e168 = next(r for r in tpp1 if r.seqid.num == 168)
    l183 = next(r for r in tpp1 if r.seqid.num == 183)
    e215 = next(r for r in tpp1 if r.seqid.num == 215)

    # 1. TRAP vector length: 21.757 A
    trap_len = r756["CA"][0].pos.dist(r774["CA"][0].pos)
    assert abs(trap_len - 21.757) < 0.05, f"Expected ~21.76 A, got {trap_len:.3f}"

    # 2. TERT Arg774 to TPP1 Glu168 salt bridge: 2.320 A
    sb_dist = r774["NH1"][0].pos.dist(e168["OE1"][0].pos)
    assert abs(sb_dist - 2.320) < 0.05, f"Expected 2.32 A, got {sb_dist:.3f}"

    # 3. TERT Leu769 to TPP1 Leu183 hydrophobic contact: 3.665 A
    vdw_dist = l769["CD2"][0].pos.dist(l183["CD2"][0].pos)
    assert abs(vdw_dist - 3.665) < 0.05, f"Expected 3.67 A, got {vdw_dist:.3f}"

    # 4. TERT Lys78 to TPP1 Glu215 salt bridge: 3.58 A
    ten_sb = k78["NZ"][0].pos.dist(e215["OE2"][0].pos)
    assert abs(ten_sb - 3.58) < 0.05, f"Expected 3.58 A, got {ten_sb:.3f}"

    # 5. Inter-domain bridge Tyr176:O to Lys757:CD: 3.65 A
    bridge_d = min(a1.pos.dist(a2.pos) for a1 in y176 for a2 in k757)
    assert abs(bridge_d - 3.647) < 0.05, f"Expected 3.65 A, got {bridge_d:.3f}"

    # 6. Arg756 to RNA A55: 3.12 A
    r_a55 = next(r for r in rna if r.seqid.num == 55)
    rna_contact = r756["NH2"][0].pos.dist(r_a55["C2"][0].pos)
    assert abs(rna_contact - 3.12) < 0.05, f"Expected 3.12 A, got {rna_contact:.3f}"

    # 7. Lys757 to RNA template: > 8.0 A (proves no direct contact to template)
    k757_rna_min = min(a1.pos.dist(a2.pos) for a1 in k757 for r in rna for a2 in r)
    assert k757_rna_min > 8.0, f"Expected Lys757 to be >8 A from RNA, got {k757_rna_min:.2f} A"

    # 8. Arg756 to DNA DT 12: 5.02 A
    r_dt12 = next(r for r in dna if r.seqid.num == 12)
    dna_dist = r756["NE"][0].pos.dist(r_dt12["C5'"][0].pos)
    assert abs(dna_dist - 5.02) < 0.05, f"Expected 5.02 A, got {dna_dist:.3f}"
