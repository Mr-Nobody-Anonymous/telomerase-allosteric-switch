"""
Unit and integration test suite for telotrade sequence mapping and residue correspondence (Milestone 2).

Covers:
- mmCIF polymer and coordinate sequence extraction via Gemmi and Biopython.
- Chain identifier disambiguation (label_asym_id vs auth_asym_id) across 7TRE, 7TRD, 7TRF, 7TRC, 2I46.
- Construct tag offset reconcilers:
    - TERT: auth_seq_num = label_seq_id - 35 (canonical UniProt O14746).
    - TPP1 in 7TRE: auth_seq_num = label_seq_id + 86 (canonical UniProt Q96AP0).
    - TPP1 in 2I46: auth_seq_num = label_seq_id + 89.
- 2I46 dimer asymmetry handling:
    - Chain A: resolved TEL-patch DE loop (residues 168-171 modeled).
    - Chain B: disordered DE loop (residues 168-171 missing coordinates).
    - Safe chain redirection and disorder detection warnings.
- Sequence-aware pairwise homology alignment:
    - 7TRE TERT <-> 7TRD TERT (934 equivalent residues).
    - 7TRE TERT <-> 7TRF TERT (934 equivalent residues).
    - 7TRE TPP1 <-> 2I46 TPP1 Chain A (152 equivalent residues).
    - 7TRE TPP1 <-> 2I46 TPP1 Chain B (134 mapped residues, DE-loop unmapped).
- UniProt residue index mapping and domain partitioning.
- Standardized CSV table generation:
    - residue_correspondence_tert.csv (1132 rows).
    - residue_correspondence_tpp1.csv (161 rows).
    - chain_entity_mapping.csv (19 rows).
"""

from pathlib import Path

import pandas as pd
import pytest

from telotrade.sequence_mapping import (
    TEL_PATCH_DE_LOOP,
    TEL_PATCH_RESIDUES,
    TERT_AUTH_OFFSET,
    TERT_DOMAINS,
    TERT_LABEL_OFFSET,
    TERT_UNIPROT_ACCESSION,
    TPP1_2I46_AUTH_OFFSET,
    TPP1_2I46_LABEL_OFFSET,
    TPP1_7TRE_AUTH_OFFSET,
    TPP1_7TRE_LABEL_OFFSET,
    TPP1_UNIPROT_ACCESSION,
    align_homologous_chains,
    align_sequences,
    check_tel_patch_integrity,
    find_project_root,
    generate_chain_entity_mapping,
    generate_tert_correspondence_table,
    generate_tpp1_correspondence_table,
    get_chain_sequence,
    get_modeled_residues,
    get_polymer_entity_sequence,
    get_tert_domain,
    get_uniprot_mapping,
    map_residue_correspondence,
    map_tert_correspondence,
    map_tpp1_correspondence,
    reconcile_label_from_auth,
    reconcile_tert_offset,
    reconcile_tpp1_offset,
    select_valid_tpp1_chain,
)


@pytest.fixture(scope="session")
def project_root() -> Path:
    return find_project_root()


@pytest.fixture(scope="session")
def raw_dir(project_root: Path) -> Path:
    path = project_root / "data" / "raw"
    assert path.is_dir(), f"Raw data directory does not exist: {path}"
    return path


# =============================================================================
# 1. Sequence Extraction & Parser Dual-Engine Verification
# =============================================================================


class TestSequenceExtraction:
    """Tests for mmCIF sequence and residue extraction via Gemmi and Biopython."""

    def test_7tre_tert_extraction_gemmi_and_biopython(self, raw_dir: Path):
        """Verify 7TRE TERT sequence extraction produces identical 934 residues across Gemmi and Biopython."""
        cif_path = raw_dir / "7TRE.cif"

        seq_g, nums_g = get_chain_sequence(cif_path, "A", use_auth=True, parser="gemmi")
        seq_b, nums_b = get_chain_sequence(cif_path, "A", use_auth=True, parser="biopython")

        assert len(seq_g) == 934
        assert len(nums_g) == 934
        assert seq_g == seq_b
        assert nums_g == nums_b
        assert seq_g.startswith(("CRVR", "C"))  # First modeled residue is Cys7
        assert nums_g[0] == 7
        assert nums_g[-1] == 1132

    def test_7tre_tpp1_extraction_gemmi_and_biopython(self, raw_dir: Path):
        """Verify 7TRE TPP1 sequence extraction produces identical 152 residues across Gemmi and Biopython."""
        cif_path = raw_dir / "7TRE.cif"

        seq_g, nums_g = get_chain_sequence(cif_path, "C", use_auth=True, parser="gemmi")
        seq_b, nums_b = get_chain_sequence(cif_path, "C", use_auth=True, parser="biopython")

        assert len(seq_g) == 152
        assert len(nums_g) == 152
        assert seq_g == seq_b
        assert nums_g == nums_b
        assert seq_g.startswith("SGRL")  # Ser90
        assert nums_g[0] == 90
        assert nums_g[-1] == 241

    def test_7tre_chain_disambiguation_author_vs_label(self, raw_dir: Path):
        """
        Verify chain disambiguation in 7TRE:
        auth 'A' is TERT (protein, 934 aa), label 'A' is hTR (RNA, 217 nt).
        auth 'B' is RNA, label 'B' is TERT.
        """
        cif_path = raw_dir / "7TRE.cif"

        # Using author chain convention
        seq_auth_a, nums_auth_a = get_chain_sequence(cif_path, "A", use_auth=True)
        assert len(seq_auth_a) == 934  # TERT protein
        assert nums_auth_a[0] == 7

        seq_auth_b, _ = get_chain_sequence(cif_path, "B", use_auth=True)
        assert len(seq_auth_b) == 217  # hTR RNA

        # Using CIF label convention
        seq_label_b, nums_label_b = get_chain_sequence(cif_path, "B", use_auth=False)
        assert len(seq_label_b) == 934  # TERT protein
        assert nums_label_b[0] == 42  # label_seq_id = 7 + 35

        seq_label_a, _ = get_chain_sequence(cif_path, "A", use_auth=False)
        assert len(seq_label_a) == 217  # hTR RNA

    def test_7tre_tpp1_chain_disambiguation_author_vs_label(self, raw_dir: Path):
        """
        Verify TPP1 chain disambiguation in 7TRE:
        auth 'C' is TPP1 (protein, 152 aa), label 'C' is DNA substrate (7 nt).
        auth 'D' is DNA, label 'D' is TPP1.
        """
        cif_path = raw_dir / "7TRE.cif"

        seq_auth_c, nums_auth_c = get_chain_sequence(cif_path, "C", use_auth=True)
        assert len(seq_auth_c) == 152  # TPP1
        assert nums_auth_c[0] == 90

        seq_auth_d, _ = get_chain_sequence(cif_path, "D", use_auth=True)
        assert len(seq_auth_d) == 7  # DNA substrate

        seq_label_d, nums_label_d = get_chain_sequence(cif_path, "D", use_auth=False)
        assert len(seq_label_d) == 152  # TPP1
        assert nums_label_d[0] == 4  # label_seq_id = 90 - 86

    def test_7trd_and_7trf_tert_lengths(self, raw_dir: Path):
        """Verify 7TRD and 7TRF TERT chains both have exactly 934 modeled residues."""
        cif_trd = raw_dir / "7TRD.cif"
        cif_trf = raw_dir / "7TRF.cif"

        seq_trd, nums_trd = get_chain_sequence(cif_trd, "A", use_auth=True)
        seq_trf, nums_trf = get_chain_sequence(cif_trf, "A", use_auth=True)

        assert len(seq_trd) == 934
        assert len(seq_trf) == 934
        assert nums_trd == nums_trf
        assert seq_trd == seq_trf

    def test_get_polymer_entity_sequence(self, raw_dir: Path):
        """Verify extraction of full unmodeled construct sequences from mmCIF entity table."""
        cif_path = raw_dir / "7TRE.cif"

        # Entity 2 is TERT (construct length 1167, containing 35-aa Twin-Strep tag)
        tert_poly = get_polymer_entity_sequence(cif_path, entity_id_or_chain=2)
        assert len(tert_poly) == 1167
        assert tert_poly.startswith("GHMSAWSHPQFEKGGGSGGGSGGSAWSHPQFEKSA")  # 35-aa Twin-Strep tag
        assert tert_poly[35:].startswith("MPRAPRCRAV")  # Canonical UniProt O14746

        # Entity 4 is TPP1 (construct length 163, containing 3-aa GAG tag)
        tpp1_poly = get_polymer_entity_sequence(cif_path, entity_id_or_chain=4)
        assert len(tpp1_poly) == 163
        assert tpp1_poly.startswith("GAGSGRLVL")

    def test_get_modeled_residues(self, raw_dir: Path):
        """Verify get_modeled_residues returns full residue objects with CA flags."""
        cif_path = raw_dir / "7TRE.cif"
        res_list = get_modeled_residues(cif_path, "A", use_auth=True)

        assert len(res_list) == 934
        first_res = res_list[0]
        assert first_res.resname == "CYS"
        assert first_res.auth_seq_num == 7
        assert first_res.label_seq_id == 42
        assert first_res.has_ca is True

    def test_nonexistent_file_raises_filenotfound(self, tmp_path: Path):
        """Verify get_chain_sequence raises FileNotFoundError for missing path."""
        with pytest.raises(FileNotFoundError):
            get_chain_sequence(tmp_path / "missing.cif", "A")

    def test_nonexistent_chain_raises_keyerror(self, raw_dir: Path):
        """Verify get_chain_sequence raises KeyError for invalid chain identifier."""
        with pytest.raises(KeyError):
            get_chain_sequence(raw_dir / "7TRE.cif", "Z")


# =============================================================================
# 2. Construct Tag Offset Reconcilers
# =============================================================================


class TestConstructTagOffsetReconcilers:
    """Tests for construct tag offset calculation and UniProt reconciliation."""

    def test_tert_offset_reconciler_constants(self):
        """Verify TERT offset constants: auth = label - 35."""
        assert TERT_AUTH_OFFSET == -35
        assert TERT_LABEL_OFFSET == 35
        assert TERT_UNIPROT_ACCESSION == "O14746"

    def test_reconcile_tert_offset(self):
        """Verify reconcile_tert_offset maps label_seq_id to canonical UniProt residue."""
        assert reconcile_tert_offset(42) == 7  # Cys7
        assert reconcile_tert_offset(1167) == 1132  # Asp1132
        assert reconcile_tert_offset(747) == 712  # Catalytic triad Asp712
        assert reconcile_tert_offset(903) == 868  # Catalytic triad Asp868
        assert reconcile_tert_offset(904) == 869  # Catalytic triad Asp869

    def test_tpp1_offset_reconciler_constants(self):
        """Verify TPP1 offset constants for 7TRE and 2I46."""
        assert TPP1_7TRE_AUTH_OFFSET == 86
        assert TPP1_7TRE_LABEL_OFFSET == -86
        assert TPP1_2I46_AUTH_OFFSET == 89
        assert TPP1_2I46_LABEL_OFFSET == -89
        assert TPP1_UNIPROT_ACCESSION == "Q96AP0"

    def test_reconcile_tpp1_offset_7tre(self):
        """Verify reconcile_tpp1_offset for 7TRE: auth = label + 86."""
        assert reconcile_tpp1_offset(4, structure_id="7TRE") == 90  # Ser90
        assert reconcile_tpp1_offset(82, structure_id="7TRE") == 168  # Glu168 (TEL-patch)
        assert reconcile_tpp1_offset(85, structure_id="7TRE") == 171  # Glu171 (TEL-patch)
        assert reconcile_tpp1_offset(97, structure_id="7TRE") == 183  # Leu183 (TEL-patch)
        assert reconcile_tpp1_offset(126, structure_id="7TRE") == 212  # Leu212 (TEL-patch)
        assert reconcile_tpp1_offset(155, structure_id="7TRE") == 241  # His241

    def test_reconcile_tpp1_offset_2i46(self):
        """Verify reconcile_tpp1_offset for 2I46: auth = label + 89."""
        assert reconcile_tpp1_offset(1, structure_id="2I46") == 90  # Ser90
        assert reconcile_tpp1_offset(79, structure_id="2I46") == 168  # Glu168
        assert reconcile_tpp1_offset(82, structure_id="2I46") == 171  # Glu171
        assert reconcile_tpp1_offset(94, structure_id="2I46") == 183  # Leu183
        assert reconcile_tpp1_offset(123, structure_id="2I46") == 212  # Leu212
        assert reconcile_tpp1_offset(152, structure_id="2I46") == 241  # His241

    def test_reconcile_label_from_auth(self):
        """Verify reverse offset conversion from author/UniProt to CIF label_seq_id."""
        assert reconcile_label_from_auth(7, protein="TERT") == 42
        assert reconcile_label_from_auth(1132, protein="TERT") == 1167
        assert reconcile_label_from_auth(168, protein="TPP1", structure_id="7TRE") == 82
        assert reconcile_label_from_auth(168, protein="TPP1", structure_id="2I46") == 79

    def test_get_uniprot_mapping_tert(self, raw_dir: Path):
        """Verify get_uniprot_mapping for TERT across both use_auth=True and use_auth=False."""
        cif_path = raw_dir / "7TRE.cif"

        map_auth = get_uniprot_mapping(cif_path, "A", use_auth=True)
        assert len(map_auth) == 934
        assert map_auth[7] == 7
        assert map_auth[1132] == 1132

        map_label = get_uniprot_mapping(cif_path, "B", use_auth=False)
        assert len(map_label) == 934
        assert map_label[42] == 7
        assert map_label[1167] == 1132

    def test_get_uniprot_mapping_tpp1(self, raw_dir: Path):
        """Verify get_uniprot_mapping for TPP1 in 7TRE and 2I46."""
        cif_7tre = raw_dir / "7TRE.cif"
        cif_2i46 = raw_dir / "2I46.cif"

        # 7TRE
        map_7tre_auth = get_uniprot_mapping(cif_7tre, "C", use_auth=True)
        assert len(map_7tre_auth) == 152
        assert map_7tre_auth[168] == 168

        map_7tre_label = get_uniprot_mapping(cif_7tre, "D", use_auth=False)
        assert len(map_7tre_label) == 152
        assert map_7tre_label[82] == 168

        # 2I46
        map_2i46_auth = get_uniprot_mapping(cif_2i46, "A", use_auth=True)
        assert len(map_2i46_auth) == 152
        assert map_2i46_auth[168] == 168

        map_2i46_label = get_uniprot_mapping(cif_2i46, "A", use_auth=False)
        assert len(map_2i46_label) == 152
        assert map_2i46_label[79] == 168


# =============================================================================
# 3. 2I46 Dimer Asymmetry & TEL-Patch Disorder
# =============================================================================


class Test2I46DimerAsymmetryAndDisorder:
    """Tests for 2I46 dimer asymmetry (Chain A resolved DE loop vs Chain B disordered)."""

    def test_2i46_chain_a_has_tel_patch_resolved(self, raw_dir: Path):
        """Verify 2I46 Chain A has all TEL-patch residues (168-171, 183, 212) resolved."""
        cif_path = raw_dir / "2I46.cif"
        res = check_tel_patch_integrity(cif_path, chain_id="A")

        assert res["chain_id"] == "A"
        assert res["resolved"] is True
        assert res["de_loop_resolved"] is True
        assert res["missing_residues"] == []
        assert set(TEL_PATCH_RESIDUES).issubset(set(res["modeled_residues"]))
        assert res["warning"] is None

    def test_2i46_chain_b_has_tel_patch_de_loop_disordered(self, raw_dir: Path):
        """Verify 2I46 Chain B has the TEL-patch DE loop (168-171) completely disordered/unmodeled."""
        cif_path = raw_dir / "2I46.cif"

        with pytest.warns(UserWarning, match="TEL-patch residue disorder detected"):
            res = check_tel_patch_integrity(cif_path, chain_id="B")

        assert res["chain_id"] == "B"
        assert res["resolved"] is False
        assert res["de_loop_resolved"] is False
        assert set(TEL_PATCH_DE_LOOP).issubset(set(res["missing_residues"]))
        # Invariant core anchors Leu183 and Leu212 ARE modeled in Chain B
        assert 183 in res["modeled_residues"]
        assert 212 in res["modeled_residues"]

    def test_7tre_tpp1_has_full_tel_patch_resolved(self, raw_dir: Path):
        """Verify 7TRE TPP1 (auth chain C) has full TEL patch resolved."""
        cif_path = raw_dir / "7TRE.cif"
        res = check_tel_patch_integrity(cif_path, chain_id="C")

        assert res["resolved"] is True
        assert res["de_loop_resolved"] is True
        assert res["missing_residues"] == []

    def test_select_valid_tpp1_chain_redirects_2i46_chain_b(self, raw_dir: Path):
        """Verify select_valid_tpp1_chain redirects Chain B to Chain A with a warning."""
        cif_path = raw_dir / "2I46.cif"

        with pytest.warns(UserWarning, match="Redirecting to 2I46 Chain A"):
            selected = select_valid_tpp1_chain(cif_path, preferred_chain="B", require_tel_patch=True)

        assert selected == "A"

    def test_select_valid_tpp1_chain_keeps_chain_a(self, raw_dir: Path):
        """Verify select_valid_tpp1_chain retains Chain A when requested."""
        cif_path = raw_dir / "2I46.cif"
        selected = select_valid_tpp1_chain(cif_path, preferred_chain="A", require_tel_patch=True)
        assert selected == "A"


# =============================================================================
# 4. Sequence-Aware Homology Alignment
# =============================================================================


class TestSequenceAwareHomologyAlignment:
    """Tests for pairwise sequence alignment and residue correspondence mapping."""

    def test_align_identical_sequences(self):
        """Verify align_sequences on identical strings yields 100% identity and 1:1 mapping."""
        s = "ACDEFGHIKLMNPQRSTVWY"
        res = align_sequences(s, s)
        assert res.identity == 1.0
        assert len(res.mapping) == len(s)
        for i in range(len(s)):
            assert res.mapping[i] == i

    def test_align_sequences_with_gap(self):
        """Verify align_sequences correctly handles insertions/deletions."""
        s1 = "ACDEFGHIK"
        s2 = "ACDFGHIK"  # Missing 'E' at index 3 in s1
        res = align_sequences(s1, s2)
        assert res.mapping[0] == 0  # A
        assert res.mapping[1] == 1  # C
        assert res.mapping[2] == 2  # D
        assert 3 not in res.mapping  # E deleted (index 3 in s1)
        assert res.mapping[4] == 3  # F
        assert res.unaligned_res1 == [3]

    def test_map_tert_7tre_vs_7trd(self, raw_dir: Path):
        """Verify map_residue_correspondence between 7TRE and 7TRD TERT maps all 934 residues."""
        cif_tre = raw_dir / "7TRE.cif"
        cif_trd = raw_dir / "7TRD.cif"

        mapping = map_residue_correspondence(cif_tre, "A", cif_trd, "A", use_auth=True)
        assert len(mapping) == 934
        assert mapping[7] == 7
        assert mapping[1132] == 1132
        assert mapping[712] == 712  # Catalytic Asp712

    def test_map_tert_7tre_vs_7trf(self, raw_dir: Path):
        """Verify map_residue_correspondence between 7TRE and 7TRF TERT maps all 934 residues."""
        cif_tre = raw_dir / "7TRE.cif"
        cif_trf = raw_dir / "7TRF.cif"

        mapping = map_residue_correspondence(cif_tre, "A", cif_trf, "A", use_auth=True)
        assert len(mapping) == 934
        assert mapping[7] == 7
        assert mapping[1132] == 1132

    def test_map_tpp1_7tre_vs_2i46_chain_a(self, raw_dir: Path):
        """
        Verify map_residue_correspondence between 7TRE TPP1 and 2I46 Chain A maps all 152 residues,
        including all TEL-patch positions (168, 169, 170, 171, 183, 212).
        """
        cif_tre = raw_dir / "7TRE.cif"
        cif_i46 = raw_dir / "2I46.cif"

        mapping = map_residue_correspondence(cif_tre, "C", cif_i46, "A", use_auth=True)
        assert len(mapping) == 152
        assert mapping[90] == 90
        assert mapping[241] == 241
        # Check all TEL-patch residues
        for tel_res in TEL_PATCH_RESIDUES:
            assert mapping[tel_res] == tel_res

    def test_map_tpp1_7tre_vs_2i46_chain_a_label_seq(self, raw_dir: Path):
        """
        Verify map_residue_correspondence under label_seq numbering reconciles the 3-aa GAG tag:
        label 82 in 7TRE (Glu168) maps to label 79 in 2I46 (Glu168).
        """
        cif_tre = raw_dir / "7TRE.cif"
        cif_i46 = raw_dir / "2I46.cif"

        mapping = map_residue_correspondence(cif_tre, "D", cif_i46, "A", use_auth=False)
        assert len(mapping) == 152
        # Ser90: 7TRE label 4 -> 2I46 label 1
        assert mapping[4] == 1
        # Glu168: 7TRE label 82 -> 2I46 label 79
        assert mapping[82] == 79
        # His241: 7TRE label 155 -> 2I46 label 152
        assert mapping[155] == 152

    def test_map_tpp1_7tre_vs_2i46_chain_b_detects_disorder(self, raw_dir: Path):
        """
        Verify map_residue_correspondence against 2I46 Chain B correctly maps only 134 residues,
        leaving TEL-patch DE-loop residues (168-171) unmapped.
        """
        cif_tre = raw_dir / "7TRE.cif"
        cif_i46 = raw_dir / "2I46.cif"

        aln_res = align_homologous_chains(cif_tre, "C", cif_i46, "B", use_auth=True)
        mapping = aln_res.mapping

        assert len(mapping) == 134
        # DE-loop residues are unmapped in Chain B
        for de_res in TEL_PATCH_DE_LOOP:
            assert de_res not in mapping
            assert de_res in aln_res.unaligned_res1
        # Invariant core hydrophobic anchors ARE mapped
        assert mapping[183] == 183
        assert mapping[212] == 212

    def test_map_tert_and_tpp1_convenience_wrappers(self, raw_dir: Path):
        """Verify convenience wrappers map_tert_correspondence and map_tpp1_correspondence."""
        cif_tre = raw_dir / "7TRE.cif"
        cif_trd = raw_dir / "7TRD.cif"
        cif_i46 = raw_dir / "2I46.cif"

        m_tert = map_tert_correspondence(cif_tre, cif_trd)
        assert len(m_tert) == 934

        m_tpp1 = map_tpp1_correspondence(cif_tre, cif_i46)
        assert len(m_tpp1) == 152


# =============================================================================
# 5. Domain Partitioning & Standardized Table Generators
# =============================================================================


class TestStandardizedTableGenerators:
    """Tests for table generators matching the project specification."""

    def test_tert_domain_definitions(self):
        """Verify get_tert_domain assigns correct domains according to PROJECT.md."""
        assert TERT_DOMAINS["TEN"] == (1, 196)
        assert TERT_DOMAINS["LINKER"] == (197, 325)
        assert TERT_DOMAINS["TRBD"] == (326, 540)
        assert TERT_DOMAINS["IFD_TRAP"] == (755, 785)
        assert TERT_DOMAINS["RT"] == (541, 935)
        assert TERT_DOMAINS["CTE"] == (936, 1132)

        assert get_tert_domain(10) == "TEN"
        assert get_tert_domain(200) == "LINKER"
        assert get_tert_domain(400) == "TRBD"
        assert get_tert_domain(600) == "RT"
        assert get_tert_domain(774) == "IFD_TRAP"  # Essential Arg774
        assert get_tert_domain(1000) == "CTE"

    def test_generate_tert_correspondence_table_schema(self, raw_dir: Path):
        """Verify generate_tert_correspondence_table generates 1132-row table with required schema."""
        df = generate_tert_correspondence_table(raw_dir)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1132
        expected_cols = [
            "uniprot_resnum",
            "resname_canonical",
            "label_seq_id_7tre",
            "auth_seq_id_7tre",
            "modeled_in_7tre",
            "label_seq_id_7trd",
            "auth_seq_id_7trd",
            "modeled_in_7trd",
            "label_seq_id_7trf",
            "auth_seq_id_7trf",
            "modeled_in_7trf",
            "domain_assignment",
        ]
        for col in expected_cols:
            assert col in df.columns, f"Missing column {col} in TERT correspondence table"

        # Check modeled counts
        assert df["modeled_in_7tre"].sum() == 934
        assert df["modeled_in_7trd"].sum() == 934
        assert df["modeled_in_7trf"].sum() == 934

        # Check unmodeled linker (179-321)
        linker_rows = df[(df["uniprot_resnum"] >= 179) & (df["uniprot_resnum"] <= 321)]
        assert not linker_rows["modeled_in_7tre"].any()

        # Check modeled anchor Cys7
        cys7 = df[df["uniprot_resnum"] == 7].iloc[0]
        assert cys7["resname_canonical"] == "CYS"
        assert cys7["auth_seq_id_7tre"] == 7
        assert cys7["label_seq_id_7tre"] == 42
        assert cys7["modeled_in_7tre"] is True or cys7["modeled_in_7tre"] == 1

    def test_generate_tpp1_correspondence_table_schema(self, raw_dir: Path):
        """Verify generate_tpp1_correspondence_table generates 161-row table with required schema."""
        df = generate_tpp1_correspondence_table(raw_dir)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 161  # residues 89 to 249
        expected_cols = [
            "uniprot_resnum",
            "resname_canonical",
            "label_seq_id_7tre",
            "auth_seq_id_7tre",
            "modeled_in_7tre",
            "label_seq_id_2i46_chain_a",
            "auth_seq_id_2i46_chain_a",
            "modeled_in_2i46_chain_a",
            "modeled_in_2i46_chain_b",
            "is_tel_patch_residue",
        ]
        for col in expected_cols:
            assert col in df.columns, f"Missing column {col} in TPP1 correspondence table"

        # Check modeled counts
        assert df["modeled_in_7tre"].sum() == 152
        assert df["modeled_in_2i46_chain_a"].sum() == 152
        assert df["modeled_in_2i46_chain_b"].sum() == 134

        # Check TEL patch rows
        tel_rows = df[df["is_tel_patch_residue"]]
        assert len(tel_rows) == 6
        assert set(tel_rows["uniprot_resnum"]) == set(TEL_PATCH_RESIDUES)

        # Check Glu168
        glu168 = df[df["uniprot_resnum"] == 168].iloc[0]
        assert glu168["resname_canonical"] == "GLU"
        assert glu168["auth_seq_id_7tre"] == 168
        assert glu168["label_seq_id_7tre"] == 82
        assert glu168["auth_seq_id_2i46_chain_a"] == 168
        assert glu168["label_seq_id_2i46_chain_a"] == 79
        assert glu168["modeled_in_7tre"] is True or glu168["modeled_in_7tre"] == 1
        assert glu168["modeled_in_2i46_chain_a"] is True or glu168["modeled_in_2i46_chain_a"] == 1
        assert glu168["modeled_in_2i46_chain_b"] is False or glu168["modeled_in_2i46_chain_b"] == 0

    def test_generate_chain_entity_mapping(self, raw_dir: Path):
        """Verify generate_chain_entity_mapping covers all 5 structures."""
        df = generate_chain_entity_mapping(raw_dir)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 19
        assert set(df["pdb_id"]) == {"7TRE", "7TRD", "7TRF", "7TRC", "2I46"}

        # 7TRE TERT entity
        tert_7tre = df[(df["pdb_id"] == "7TRE") & (df["label_asym_id"] == "B")].iloc[0]
        assert tert_7tre["auth_asym_id"] == "A"
        assert tert_7tre["uniprot_accession"] == "O14746"
        assert tert_7tre["modeled_residues_count"] == 934
