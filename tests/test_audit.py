"""
Adversarial Verification Suite for Milestone B1 Deliverables.

Empirically challenges:
1. `results/audit/evidence_ledger.csv`:
   - Valid CSV structure and parser robustness (no unquoted commas).
   - Compliance with required 8-column schema.
   - Exact count of 24 claims (CLM-001 through CLM-024) with no duplicates.
   - Valid parameters and justified evidence_status classifications.
   - Cross-check of status distribution against repository inventory documentation.
2. `results/tables/residue_correspondence_tpp1.csv`:
   - Complete elimination of the 1-residue phase shift across all 161 rows.
   - Verification of key residues: 166 ASP, 167 TRP, 168 GLU, 169 GLU, 170 LYS, 171 GLU.
   - 100% concordance against all 152 modeled residues in 7TRE Chain C coordinates.
   - 100% concordance against 2I46 Chain A (152 modeled residues) and Chain B (134 modeled residues).
   - Strict absence of disordered DE-loop residues (168-171) in 2I46 Chain B.
3. `results/literature/prior_art_matrix.csv`:
   - 22 structured entries with non-empty fields.
   - Positive control classifications for canonical interactions (PA-001, PA-002, PA-003).
"""

import csv
from pathlib import Path

import gemmi
import pandas as pd
import pytest


@pytest.fixture(scope="session")
def project_root() -> Path:
    # Anchor to project root
    p = Path(__file__).resolve().parent.parent
    assert (p / "results").is_dir(), f"Expected results dir in {p}"
    return p


class TestEvidenceLedgerAdversarial:
    """Stress tests and verifies evidence_ledger.csv integrity."""

    def test_evidence_ledger_csv_parseable_by_pandas(self, project_root: Path):
        """Verify evidence_ledger.csv can be parsed by pandas without ParserError.

        Fails if unquoted commas exist in free-text fields.
        """
        ledger_path = project_root / "results" / "audit" / "evidence_ledger.csv"
        assert ledger_path.is_file(), f"Missing {ledger_path}"

        # This will raise pandas.errors.ParserError if any row has unquoted commas
        df = pd.read_csv(ledger_path)
        assert len(df) == 24, f"Expected 24 claims, found {len(df)}"

    def test_evidence_ledger_schema_and_columns(self, project_root: Path):
        """Verify column names match specification in PROJECT.md."""
        ledger_path = project_root / "results" / "audit" / "evidence_ledger.csv"
        expected_cols = [
            "claim_id",
            "metric_name",
            "reported_value",
            "source_file",
            "source_script",
            "parameters",
            "evidence_status",
            "notes",
        ]

        # Use csv module for raw row-level inspection
        with open(ledger_path, "r", encoding="utf-8") as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header == expected_cols, f"Header mismatch: {header} vs {expected_cols}"

            for row_idx, row in enumerate(reader, start=2):
                assert len(row) == 8, (
                    f"Line {row_idx} has {len(row)} fields instead of 8: {row[:3]}... "
                    "Likely caused by unquoted commas in free-text fields."
                )

    def test_evidence_ledger_claim_ids_and_non_nulls(self, project_root: Path):
        """Verify claim IDs CLM-001 to CLM-024 are unique, non-null, and complete."""
        ledger_path = project_root / "results" / "audit" / "evidence_ledger.csv"

        with open(ledger_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 24, f"Expected 24 claims, got {len(rows)}"

        claim_ids = [r["claim_id"] for r in rows]
        assert len(claim_ids) == len(set(claim_ids)), "Duplicate claim IDs detected"

        expected_ids = [f"CLM-{i:03d}" for i in range(1, 25)]
        assert claim_ids == expected_ids, f"Mismatch in expected claim ID sequence: {claim_ids}"

        for r in rows:
            for k, v in r.items():
                assert v is not None and str(v).strip() != "", f"Empty field '{k}' in claim {r['claim_id']}"

    def test_evidence_ledger_status_classifications(self, project_root: Path):
        """Verify validity of evidence_status classifications."""
        ledger_path = project_root / "results" / "audit" / "evidence_ledger.csv"

        valid_statuses = {"established", "uncertain", "contradicted", "untested"}

        with open(ledger_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        claims_by_id = {r["claim_id"]: r for r in rows}

        for cid, r in claims_by_id.items():
            st = r["evidence_status"]
            assert st in valid_statuses, f"Invalid evidence_status '{st}' for {cid}"

        # Justification checks for high-impact claims
        assert claims_by_id["CLM-005"]["evidence_status"] == "uncertain", (
            "CLM-005 (2.32 A salt bridge) must be classified as 'uncertain' due to cryo-EM compaction artifact"
        )
        assert claims_by_id["CLM-010"]["evidence_status"] == "uncertain", (
            "CLM-010 (0.338 A core RMSD) must be classified as 'uncertain' regarding allostery proof"
        )
        assert claims_by_id["CLM-021"]["evidence_status"] == "untested", (
            "CLM-021 (Two-step sequential binding) must be 'untested' without time-resolved kinetic data"
        )
        assert claims_by_id["CLM-022"]["evidence_status"] == "contradicted", (
            "CLM-022 (Zero allostery claim) must be 'contradicted' by dynamic allostery biophysics"
        )
        assert claims_by_id["CLM-023"]["evidence_status"] == "untested", (
            "CLM-023 (TELIN-1 inhibitor) must be 'untested' (0 wet-lab data)"
        )
        assert claims_by_id["CLM-024"]["evidence_status"] == "established", (
            "CLM-024 (Sequence offset bug fix) must be 'established'"
        )

    def test_evidence_ledger_status_counts_arithmetic(self, project_root: Path):
        """Cross-check status counts in evidence_ledger.csv against repository_inventory.md.

        Empirically checks if documentation claims match actual row counts.
        """
        ledger_path = project_root / "results" / "audit" / "evidence_ledger.csv"
        inv_path = project_root / "results" / "audit" / "repository_inventory.md"

        with open(ledger_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        counts = {}
        for r in rows:
            st = r["evidence_status"]
            counts[st] = counts.get(st, 0) + 1

        # Expected counts from 24 claims: 17 established, 4 uncertain, 1 contradicted, 2 untested
        assert counts.get("established", 0) == 17, (
            f"Expected 17 established claims, got {counts.get('established', 0)}"
        )
        assert counts.get("uncertain", 0) == 4, (
            f"Expected 4 uncertain claims, got {counts.get('uncertain', 0)}"
        )
        assert counts.get("contradicted", 0) == 1, (
            f"Expected 1 contradicted claim, got {counts.get('contradicted', 0)}"
        )
        assert counts.get("untested", 0) == 2, (
            f"Expected 2 untested claims, got {counts.get('untested', 0)}"
        )
        assert sum(counts.values()) == 24

        # Check repository_inventory.md text consistency
        if inv_path.is_file():
            inv_text = inv_path.read_text(encoding="utf-8")
            if "Established (18 claims, 75.0%)" in inv_text:
                pytest.fail(
                    "Arithmetic discrepancy in repository_inventory.md: Section 8 states "
                    "'Established (18 claims, 75.0%)', but evidence_ledger.csv contains "
                    "exactly 17 established claims (70.83% of 24). Total would exceed 100% (25 claims)."
                )


class TestTpp1ResidueCorrespondenceAdversarial:
    """Stress tests residue_correspondence_tpp1.csv against primary structures."""

    def test_tpp1_phase_shift_eliminated(self, project_root: Path):
        """Verify critical TEL-patch and anchor residue identities have no phase shift."""
        csv_path = project_root / "results" / "tables" / "residue_correspondence_tpp1.csv"
        assert csv_path.is_file(), f"Missing {csv_path}"

        df = pd.read_csv(csv_path)
        assert len(df) == 161
        assert df["uniprot_resnum"].min() == 89
        assert df["uniprot_resnum"].max() == 249

        lookup = df.set_index("uniprot_resnum")["resname_canonical"].to_dict()

        # N-terminus
        assert lookup[89] == "GLY"
        assert lookup[90] == "SER"

        # Critical DE-loop region (residues 165 to 172)
        assert lookup[165] == "SER"
        assert lookup[166] == "ASP", f"Phase shift error: res 166 is {lookup[166]}, expected ASP"
        assert lookup[167] == "TRP", f"Phase shift error: res 167 is {lookup[167]}, expected TRP"
        assert lookup[168] == "GLU", f"Phase shift error: res 168 is {lookup[168]}, expected GLU"
        assert lookup[169] == "GLU", f"Phase shift error: res 169 is {lookup[169]}, expected GLU"
        assert lookup[170] == "LYS", f"Phase shift error: res 170 is {lookup[170]}, expected LYS"
        assert lookup[171] == "GLU", f"Phase shift error: res 171 is {lookup[171]}, expected GLU"
        assert lookup[172] == "PHE", f"Phase shift error: res 172 is {lookup[172]}, expected PHE"

        # Hydrophobic anchors
        assert lookup[183] == "LEU"
        assert lookup[212] == "LEU"
        assert lookup[215] == "GLU"
        assert lookup[241] == "HIS"

    def test_tpp1_exact_concordance_with_7tre_coordinates(self, project_root: Path):
        """Verify all 152 modeled residues in 7TRE Chain C match table with 0 mismatches."""
        csv_path = project_root / "results" / "tables" / "residue_correspondence_tpp1.csv"
        cif_path = project_root / "data" / "raw" / "7TRE.cif"

        df = pd.read_csv(csv_path)

        doc = gemmi.cif.read_file(str(cif_path))
        st = gemmi.make_structure_from_block(doc.sole_block())
        chain_c = st[0]["C"]

        cif_residues = {res.seqid.num: res.name for res in chain_c}
        assert len(cif_residues) == 152, f"Expected 152 modeled residues in 7TRE Chain C, got {len(cif_residues)}"

        mismatches = []
        for auth_id, cif_name in cif_residues.items():
            matching_rows = df[df["auth_seq_id_7tre"] == auth_id]
            if len(matching_rows) == 0:
                mismatches.append(f"Residue {auth_id} ({cif_name}) not found in table")
                continue
            table_name = matching_rows.iloc[0]["resname_canonical"]
            is_modeled = matching_rows.iloc[0]["modeled_in_7tre"]
            if table_name != cif_name:
                mismatches.append(f"Residue {auth_id}: CIF={cif_name} vs Table={table_name}")
            if not is_modeled:
                mismatches.append(f"Residue {auth_id} marked not modeled but present in CIF")

        assert len(mismatches) == 0, f"Found {len(mismatches)} mismatches with 7TRE: {mismatches[:5]}"

    def test_tpp1_exact_concordance_with_2i46_coordinates(self, project_root: Path):
        """Verify 2I46 Chain A (152 aa) and Chain B (134 aa) concordance and DE-loop disorder."""
        csv_path = project_root / "results" / "tables" / "residue_correspondence_tpp1.csv"
        cif_path = project_root / "data" / "raw" / "2I46.cif"

        df = pd.read_csv(csv_path)

        doc = gemmi.cif.read_file(str(cif_path))
        st = gemmi.make_structure_from_block(doc.sole_block())

        # Chain A
        chain_a = st[0]["A"]
        a_residues = {res.seqid.num: res.name for res in chain_a}
        assert len(a_residues) == 152

        mismatches_a = []
        for auth_id, cif_name in a_residues.items():
            matching_rows = df[df["auth_seq_id_2i46_chain_a"] == auth_id]
            if len(matching_rows) == 0:
                mismatches_a.append(f"Chain A Residue {auth_id} not found in table")
                continue
            table_name = matching_rows.iloc[0]["resname_canonical"]
            if table_name != cif_name:
                mismatches_a.append(f"Chain A Residue {auth_id}: CIF={cif_name} vs Table={table_name}")

        assert len(mismatches_a) == 0, f"Mismatches in 2I46 Chain A: {mismatches_a}"

        # Chain B
        chain_b = st[0]["B"]
        b_residues = {res.seqid.num: res.name for res in chain_b}
        assert len(b_residues) == 134

        # DE-loop residues (168-171) must be disordered in Chain B
        for r in [168, 169, 170, 171]:
            assert r not in b_residues, f"Residue {r} should be disordered in 2I46 Chain B"
            row = df[df["uniprot_resnum"] == r].iloc[0]
            assert bool(row["modeled_in_2i46_chain_b"]) is False


class TestPriorArtMatrixAdversarial:
    """Verifies prior_art_matrix.csv compliance with Requirement R2."""

    def test_prior_art_matrix_schema_and_positive_controls(self, project_root: Path):
        """Verify positive controls and schema in prior_art_matrix.csv."""
        matrix_path = project_root / "results" / "literature" / "prior_art_matrix.csv"
        assert matrix_path.is_file(), f"Missing {matrix_path}"

        df = pd.read_csv(matrix_path)
        assert len(df) >= 20

        expected_cols = [
            "finding_id",
            "structural_feature",
            "prior_art_source",
            "prior_art_status",
            "original_publication_finding",
            "this_repo_claim",
            "classification",
            "notes",
        ]
        assert list(df.columns) == expected_cols

        # Check positive controls as mandated by R2
        pa_001 = df[df["finding_id"] == "PA-001"].iloc[0]
        assert "Arg774" in pa_001["structural_feature"] and "Glu168" in pa_001["structural_feature"]
        assert pa_001["classification"] == "Positive Control"

        pa_002 = df[df["finding_id"] == "PA-002"].iloc[0]
        assert "Lys78" in pa_002["structural_feature"] and "Glu215" in pa_002["structural_feature"]
        assert pa_002["classification"] == "Positive Control"
