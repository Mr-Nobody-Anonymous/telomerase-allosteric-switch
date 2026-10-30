"""
Adversarial stress-testing suite for Milestone 2:
Sequence-Aware Chain & Residue Correspondence Engine.
Authored by Challenger 1 (challenger_m2_1).

Empirically probes:
1. Sequence alignment resilience to empty sequences, gaps, insertions, terminal truncations,
   alanine scanning mutations, and completely non-homologous sequences.
2. Construct tag offset reconcilers under negative, zero, and out-of-range residue indices.
3. Chain disambiguation and UniProt mapping logic across heterogeneous macromolecular entities
   (TERT, TPP1, hTR RNA, telomeric DNA, nucleosomal histones).
4. Gemmi vs Biopython parser parity and failure modes under label_asym_id lookups.
5. Strict preservation of 2I46 Chain A (resolved DE-loop) vs Chain B (disordered DE-loop) asymmetry.
"""

from pathlib import Path

import pytest

from telotrade.sequence_mapping import (
    TEL_PATCH_DE_LOOP,
    TEL_PATCH_RESIDUES,
    align_homologous_chains,
    align_sequences,
    check_tel_patch_integrity,
    find_project_root,
    get_chain_sequence,
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


@pytest.fixture(scope="module")
def project_root() -> Path:
    return find_project_root()


@pytest.fixture(scope="module")
def raw_dir(project_root: Path) -> Path:
    p = project_root / "data" / "raw"
    assert p.is_dir()
    return p


# =============================================================================
# Dimension 1: Sequence Alignment with Gaps, Insertions, Mutations, & Non-Homology
# =============================================================================


class TestAdversarialAlignmentRobustness:
    """Empirically stress-test pairwise alignment and coordinate mapping."""

    def test_empty_sequence_handled_gracefully_with_empty_result(self):
        """
        Verify that align_sequences returns an empty AlignmentResult with score=0.0,
        identity=0.0, and mapping={} when an empty string is provided, rather than
        raising ValueError.
        """
        res1 = align_sequences("", "ACDEFGHIKL")
        assert res1.score == 0.0
        assert res1.identity == 0.0
        assert res1.mapping == {}
        assert res1.reverse_mapping == {}
        assert res1.aligned_pairs == []
        assert res1.unaligned_res2 == list(range(10))

        res2 = align_sequences("ACDEFGHIKL", "")
        assert res2.score == 0.0
        assert res2.identity == 0.0
        assert res2.mapping == {}
        assert res2.reverse_mapping == {}
        assert res2.aligned_pairs == []
        assert res2.unaligned_res1 == list(range(10))

        res3 = align_sequences("", "")
        assert res3.score == 0.0
        assert res3.identity == 0.0
        assert res3.mapping == {}

    def test_single_and_multiple_internal_insertions_preserve_indices(self):
        """
        Test that insertions (gaps in seq1 or seq2) do not induce index drift
        in flanking or downstream aligned blocks.
        """
        # seq1 has 20 aa; seq2 has a 5-aa insertion in the middle
        s1 = "ACDEFGHIKLMNPQRSTVWY"
        s2 = "ACDEFXXXXXGHIKLMNPQRSTVWY"

        res = align_sequences(s1, s2)
        assert res.score > 0
        # Residues 0..4 in s1 map to 0..4 in s2
        for i in range(5):
            assert res.mapping[i] == i

        # Residues 5..19 in s1 map to 10..24 in s2 (exact 5-residue shift matching insertion)
        for i in range(5, 20):
            assert res.mapping[i] == i + 5

        # Unaligned residues in s2 must be exactly the inserted block [5, 6, 7, 8, 9]
        assert res.unaligned_res2 == [5, 6, 7, 8, 9]
        assert res.unaligned_res1 == []

    def test_double_internal_gaps_maintain_downstream_fidelity(self):
        """Test two distinct deletion events separated by structured regions."""
        s1 = "AAAAABBBBBCCCCCDDDDDEEEEE"
        # s2 has deletions of BBBBB and DDDDD
        s2 = "AAAAACCCCCEEEEE"

        res = align_sequences(s1, s2)
        # Block A: 0..4 -> 0..4
        for i in range(5):
            assert res.mapping[i] == i
        # Block C: 10..14 in s1 -> 5..9 in s2
        for i in range(10, 15):
            assert res.mapping[i] == i - 5
        # Block E: 20..24 in s1 -> 10..14 in s2
        for i in range(20, 25):
            assert res.mapping[i] == i - 10

        # Unaligned residues in s1 must be the two deleted blocks
        assert sorted(res.unaligned_res1) == [5, 6, 7, 8, 9, 15, 16, 17, 18, 19]

    def test_terminal_truncations_do_not_shift_remaining_blocks(self):
        """Test severe N-terminal and C-terminal truncations without terminal character collisions."""
        s_full = "MSTNPKPQRKWTHRLVRPWIEEVRAL"
        s_trunc_n = s_full[10:]  # 10 aa N-term truncated: KWTHRLVRPWIEEVRAL
        s_trunc_c = s_full[:14]  # 12 aa C-term truncated: MSTNPKPQRKWTHR

        res_n = align_sequences(s_full, s_trunc_n)
        for i in range(10, len(s_full)):
            assert res_n.mapping[i] == i - 10

        res_c = align_sequences(s_full, s_trunc_c)
        for i in range(14):
            assert res_c.mapping[i] == i

    def test_global_alignment_terminal_gap_stretching_artifact(self):
        """
        FINDING 5 EMPIRICAL PROBE:
        In global alignment mode (Align.PairwiseAligner mode='global' without free end gaps),
        if a truncated query ends with a residue that matches the terminal residue of the target,
        the global aligner stretches an internal gap across the middle rather than leaving
        a terminal gap, because matching the final residue scores higher than end-gap penalization.
        """
        s_full = "MSTNPKPQRKWTHRLVRPWIEEVRAL"  # Ends with 'L' at index 25; index 14 is also 'L'
        s_trunc_c = s_full[:15]               # Ends with 'L' at index 14

        res = align_sequences(s_full, s_trunc_c)
        # Because of global end-gap penalty, query index 14 ('L') is aligned to target index 25 ('L')!
        assert res.mapping[25] == 14
        # Target index 14 is left unaligned in the internal gap
        assert 14 in res.unaligned_res1

    def test_alanine_scanning_mutations_pair_homologous_positions(self):
        """
        Verify that multi-site point mutations (e.g. TEL-patch alanine scan)
        maintain 1:1 residue correspondence without opening spurious gaps.
        """
        wt_patch = "TFTGSERRDEEEELVDRN"
        # Alanine scan on DEEEEL -> AAAAAA
        mut_patch = "TFTGSERRAAAAAAVDRN"

        res = align_sequences(wt_patch, mut_patch)
        assert len(res.mapping) == len(wt_patch)
        # All positions should map 1:1 without gaps
        for i in range(len(wt_patch)):
            assert res.mapping[i] == i
        # Identity should reflect 6 mutated residues out of 18 (66.7% match)
        assert pytest.approx(res.identity, 0.01) == 12 / 18

    def test_non_homologous_chains_produce_low_identity_negative_score(self, raw_dir: Path):
        """
        Empirically challenge map_residue_correspondence on completely non-homologous
        chains (TERT vs Histone H2B).
        """
        cif_7tre = raw_dir / "7TRE.cif"
        cif_7trf = raw_dir / "7TRF.cif"

        res = align_homologous_chains(
            cif_path_1=cif_7tre,
            chain_1="A",  # TERT
            cif_path_2=cif_7trf,
            chain_2="F",  # Histone H2B (non-homologous)
            use_auth=True,
        )
        # Sequence identity between TERT and Histone is minimal (< 5%)
        assert res.identity < 0.05
        # Alignment score is negative
        assert res.score < 0

        # Finding: map_residue_correspondence silently returns spurious mapping without error
        mapping = map_residue_correspondence(cif_7tre, "A", cif_7trf, "F", use_auth=True)
        assert isinstance(mapping, dict)
        assert len(mapping) > 0  # Spurious mapping returned for non-homologous proteins


# =============================================================================
# Dimension 2: Construct Tag Offset Reconcilers & Out-of-Range Handling
# =============================================================================


class TestAdversarialOffsetReconcilers:
    """Stress-test offset translations with boundary, negative, and invalid indices."""

    def test_reconcile_tert_offset_boundary_and_negative_indices(self):
        """
        Test TERT offset reconciler under out-of-range and negative values.
        Construct: label 1..35 is Twin-Strep tag; label 42..1167 is modeled TERT (auth 7..1132).
        """
        # Exactly at tag boundary
        assert reconcile_tert_offset(35) == 0  # UniProt has no residue 0!
        assert reconcile_tert_offset(36) == 1  # First UniProt residue
        assert reconcile_tert_offset(42) == 7  # First modeled residue in 7TRE

        # Tag residues (label 1..35) produce negative UniProt sequence numbers
        assert reconcile_tert_offset(1) == -34
        assert reconcile_tert_offset(10) == -25

        # Negative and zero inputs produce unphysical outputs
        assert reconcile_tert_offset(0) == -35
        assert reconcile_tert_offset(-10) == -45

        # Upper bound overflow beyond full-length 1132
        assert reconcile_tert_offset(1167) == 1132
        assert reconcile_tert_offset(2000) == 1965

    def test_reconcile_tpp1_offset_boundary_and_negative_indices(self):
        """
        Test TPP1 offset reconciler under boundary and negative values.
        In 7TRE: auth = label + 86.
        In 2I46: auth = label + 89.
        """
        # 7TRE construct tag: label 1..3 is GAG linker (auth 87..89)
        assert reconcile_tpp1_offset(1, "7TRE") == 87
        assert reconcile_tpp1_offset(3, "7TRE") == 89
        assert reconcile_tpp1_offset(4, "7TRE") == 90  # First modeled residue

        # 2I46 construct starts at Ser90 (label 1)
        assert reconcile_tpp1_offset(1, "2I46") == 90
        assert reconcile_tpp1_offset(79, "2I46") == 168  # Glu168

        # Negative and zero inputs
        assert reconcile_tpp1_offset(0, "7TRE") == 86
        assert reconcile_tpp1_offset(-50, "7TRE") == 36
        assert reconcile_tpp1_offset(-100, "2I46") == -11

    def test_reconcile_label_from_auth_produces_negative_labels(self):
        """
        FINDING 4 EMPIRICAL PROBE:
        reconcile_label_from_auth produces negative CIF label_seq_id values
        when given valid author residue numbers outside the construct range.
        mmCIF requires label_seq_id >= 1.
        """
        # TPP1 in 7TRE: auth - 86 -> for residue 50: 50 - 86 = -36 (illegal negative CIF label)
        neg_label_7tre = reconcile_label_from_auth(50, protein="TPP1", structure_id="7TRE")
        assert neg_label_7tre == -36

        # TPP1 in 2I46: auth - 89 -> for residue 10: 10 - 89 = -79
        neg_label_2i46 = reconcile_label_from_auth(10, protein="TPP1", structure_id="2I46")
        assert neg_label_2i46 == -79

        # TERT: auth + 35 -> for negative auth index: -10 + 35 = 25
        assert reconcile_label_from_auth(-10, protein="TERT") == 25

    def test_get_tert_domain_out_of_bounds_handling(self):
        """Verify get_tert_domain gracefully handles out-of-bounds inputs."""
        assert get_tert_domain(0) == "UNKNOWN"
        assert get_tert_domain(-100) == "UNKNOWN"
        assert get_tert_domain(1133) == "UNKNOWN"
        assert get_tert_domain(9999) == "UNKNOWN"

        # Structural domains inside valid range
        assert get_tert_domain(1) == "TEN"
        assert get_tert_domain(196) == "TEN"
        assert get_tert_domain(197) == "LINKER"
        assert get_tert_domain(325) == "LINKER"
        assert get_tert_domain(326) == "TRBD"
        assert get_tert_domain(540) == "TRBD"
        assert get_tert_domain(755) == "IFD_TRAP"
        assert get_tert_domain(785) == "IFD_TRAP"
        assert get_tert_domain(754) == "RT"
        assert get_tert_domain(786) == "RT"
        assert get_tert_domain(936) == "CTE"
        assert get_tert_domain(1132) == "CTE"


# =============================================================================
# Dimension 3: UniProt Mapping & Cross-Entity Chain Conflation
# =============================================================================


class TestAdversarialUniProtMapping:
    """
    Empirically probe get_uniprot_mapping for entity conflation bugs.
    """

    def test_7tre_rna_chain_returns_empty_mapping_not_negative_uniprot(self, raw_dir: Path):
        """
        Verify that in 7TRE, Label Chain A (hTR RNA) is not misclassified as TERT,
        and returns an empty mapping instead of generating negative UniProt residue numbers.
        """
        cif_7tre = raw_dir / "7TRE.cif"

        # Query label chain A (which is RNA) with use_auth=False
        mapping_label = get_uniprot_mapping(cif_7tre, chain_id="A", use_auth=False)
        assert mapping_label == {}

        # Query auth chain B (which is RNA) with use_auth=True
        mapping_auth = get_uniprot_mapping(cif_7tre, chain_id="B", use_auth=True)
        assert mapping_auth == {}

    def test_7tre_dna_chain_returns_empty_mapping_not_tpp1_uniprot(self, raw_dir: Path):
        """
        Verify that in 7TRE, Label Chain C (DNA) is not misclassified as TPP1,
        and returns an empty mapping instead of corrupting DNA residue numbering.
        """
        cif_7tre = raw_dir / "7TRE.cif"

        # Query label chain C (which is DNA) with use_auth=False
        mapping_label = get_uniprot_mapping(cif_7tre, chain_id="C", use_auth=False)
        assert mapping_label == {}

        # Query auth chain D (which is DNA) with use_auth=True
        mapping_auth = get_uniprot_mapping(cif_7tre, chain_id="D", use_auth=True)
        assert mapping_auth == {}

    def test_7trf_histone_returns_empty_mapping_not_tert_uniprot(self, raw_dir: Path):
        """
        Verify that in 7TRF, Histone chains are not misclassified as TERT,
        and return an empty mapping.
        """
        cif_7trf = raw_dir / "7TRF.cif"

        # Query label chain A (Histone H2B) with use_auth=False
        mapping_label = get_uniprot_mapping(cif_7trf, chain_id="A", use_auth=False)
        assert mapping_label == {}

        # Query auth chain F (Histone H2B) with use_auth=True
        mapping_auth = get_uniprot_mapping(cif_7trf, chain_id="F", use_auth=True)
        assert mapping_auth == {}


# =============================================================================
# Dimension 4: Gemmi vs Biopython Parity and Label Discrepancy
# =============================================================================


class TestAdversarialParserDiscrepancies:
    """Stress-test Gemmi vs Biopython parser parity."""

    def test_biopython_parser_honors_label_asym_id_parity_with_gemmi(self, raw_dir: Path):
        """
        Verify that Biopython MMCIFParser honors use_auth=False and achieves complete
        parity with Gemmi on label_asym_id lookup and label_seq_id extraction.
        """
        cif_7tre = raw_dir / "7TRE.cif"

        seq_gemmi, nums_gemmi = get_chain_sequence(cif_7tre, "B", use_auth=False, parser="gemmi")
        seq_bio, nums_bio = get_chain_sequence(cif_7tre, "B", use_auth=False, parser="biopython")

        # Gemmi and Biopython both return TERT protein
        assert len(seq_gemmi) == 934
        assert len(seq_bio) == 934
        assert seq_gemmi.startswith("CRAVRSLLRS")
        assert seq_bio.startswith("CRAVRSLLRS")
        assert nums_gemmi[0] == 42  # label_seq_id
        assert nums_bio[0] == 42

        # Proves total parity between parsers when use_auth=False
        assert seq_gemmi == seq_bio
        assert nums_gemmi == nums_bio


# =============================================================================
# Dimension 5: 2I46 Dimer Asymmetry & TEL-Patch Preservation
# =============================================================================


class TestAdversarial2I46DimerAsymmetry:
    """Empirically challenge 2I46 Chain A vs Chain B asymmetry."""

    def test_2i46_chain_a_and_b_disorder_strictly_preserved(self, raw_dir: Path):
        """
        Verify that Chain A has 152 residues and fully modeled DE-loop (168-171),
        while Chain B has 134 residues and completely missing DE-loop (168-171)
        and missing beta-strand loop (126-139).
        """
        cif_2i46 = raw_dir / "2I46.cif"

        _, nums_a = get_chain_sequence(cif_2i46, "A", use_auth=True)
        _, nums_b = get_chain_sequence(cif_2i46, "B", use_auth=True)

        assert len(nums_a) == 152
        assert len(nums_b) == 134

        # Chain A contains all TEL-patch residues
        for r in TEL_PATCH_RESIDUES:
            assert r in nums_a

        # Chain B is missing exactly 168-171
        for r in TEL_PATCH_DE_LOOP:
            assert r not in nums_b

        # Missing residues in Chain B are exactly [126..139] and [168..171]
        missing_in_b = sorted(set(nums_a) - set(nums_b))
        expected_missing = list(range(126, 140)) + list(range(168, 172))
        assert missing_in_b == expected_missing

    def test_2i46_a_vs_b_pairwise_alignment_has_zero_index_shift(self, raw_dir: Path):
        """
        Verify that pairwise alignment of Chain A vs Chain B handles the two gaps
        (126-139 and 168-171) without introducing ANY off-by-one or index shifts
        in the remaining 134 aligned residues.
        """
        cif_2i46 = raw_dir / "2I46.cif"

        res = align_homologous_chains(cif_2i46, "A", cif_2i46, "B", use_auth=True)
        assert len(res.mapping) == 134

        # Every single mapped residue must match its identical counterpart: r1 == r2
        for r1, r2 in res.mapping.items():
            assert r1 == r2, f"Index shift detected at residue {r1} -> {r2}!"

        # Unaligned residues in Chain A must be exactly the 18 missing residues
        assert sorted(res.unaligned_res1) == list(range(126, 140)) + list(range(168, 172))
        assert res.unaligned_res2 == []

    def test_select_valid_tpp1_chain_redirection_with_warning(self, raw_dir: Path):
        """Verify select_valid_tpp1_chain redirects 2I46 Chain B to A and emits UserWarning."""
        cif_2i46 = raw_dir / "2I46.cif"

        with pytest.warns(UserWarning, match="2I46 Chain B has disordered TEL-patch DE-loop"):
            ch = select_valid_tpp1_chain(cif_2i46, preferred_chain="B", require_tel_patch=True)
        assert ch == "A"

        # When require_tel_patch=False, allows Chain B
        ch_b = select_valid_tpp1_chain(cif_2i46, preferred_chain="B", require_tel_patch=False)
        assert ch_b == "B"

    def test_check_tel_patch_integrity_distinguishes_a_and_b(self, raw_dir: Path):
        """Verify check_tel_patch_integrity correctly evaluates resolved vs disordered status."""
        cif_2i46 = raw_dir / "2I46.cif"

        info_a = check_tel_patch_integrity(cif_2i46, "A")
        assert info_a["resolved"] is True
        assert info_a["de_loop_resolved"] is True
        assert info_a["missing_residues"] == []

        with pytest.warns(UserWarning, match="TEL-patch residue disorder detected"):
            info_b = check_tel_patch_integrity(cif_2i46, "B")
        assert info_b["resolved"] is False
        assert info_b["de_loop_resolved"] is False
        assert info_b["missing_residues"] == [168, 169, 170, 171]


# =============================================================================
# Dimension 6: Full Multi-Chain Homology Correspondence Across Targets
# =============================================================================


class TestAdversarialMultiTargetCorrespondence:
    """Verify correspondence stability across 7TRE, 7TRD, 7TRF, and 2I46."""

    def test_tert_correspondence_100_percent_identity_across_all_934_residues(self, raw_dir: Path):
        """
        Verify that TERT modeled residues in 7TRE, 7TRD, and 7TRF align with 100% identity
        and zero index drift across all 934 positions.
        """
        cif_7tre = raw_dir / "7TRE.cif"
        cif_7trd = raw_dir / "7TRD.cif"
        cif_7trf = raw_dir / "7TRF.cif"

        map_trd = map_tert_correspondence(cif_7tre, cif_7trd)
        map_trf = map_tert_correspondence(cif_7tre, cif_7trf)

        assert len(map_trd) == 934
        assert len(map_trf) == 934

        for r_tre, r_trd in map_trd.items():
            assert r_tre == r_trd

        for r_tre, r_trf in map_trf.items():
            assert r_tre == r_trf

    def test_tpp1_7tre_to_2i46_a_maps_all_152_residues_with_100_percent_identity(self, raw_dir: Path):
        """
        Verify that TPP1 in 7TRE (152 modeled residues) maps completely to 2I46 Chain A
        (152 modeled residues) with 1:1 identity across all positions 90-241.
        """
        cif_7tre = raw_dir / "7TRE.cif"
        cif_2i46 = raw_dir / "2I46.cif"

        mapping = map_tpp1_correspondence(cif_7tre, cif_2i46)
        assert len(mapping) == 152

        for r_tre, r_2i46 in mapping.items():
            assert r_tre == r_2i46
            assert 90 <= r_tre <= 241
