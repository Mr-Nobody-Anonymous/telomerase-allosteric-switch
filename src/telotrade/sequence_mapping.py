"""
Sequence Mapping and Residue Correspondence Engine for Human Telomerase.

Provides automated mmCIF polymer and coordinate sequence extraction via Gemmi and
Biopython, expression tag offset reconcilers, 2I46 dimer asymmetry handling (DE-loop
disorder vs resolved Chain A), sequence-aware homology alignment, UniProt mapping,
and generation of standardized residue correspondence tables across comparative
telomerase complexes (7TRE, 7TRD, 7TRF, 7TRC, 2I46).
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import gemmi
import pandas as pd
from Bio import Align
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Canonical Offsets and Identifiers
# -----------------------------------------------------------------------------

# TERT: Construct contains N-terminal 35-aa Twin-Strep tag (residues 1-35 in construct).
# auth_seq_num = label_seq_id - 35, matching canonical UniProt O14746 (1-1132).
TERT_AUTH_OFFSET: int = -35
TERT_LABEL_OFFSET: int = 35

# TPP1 in 7TRE: Construct contains 3-aa GAG tag (Gly87-Ala88-Gly89 in construct).
# auth_seq_num = label_seq_id + 86, matching canonical UniProt Q96AP0 (89-249).
TPP1_7TRE_AUTH_OFFSET: int = 86
TPP1_7TRE_LABEL_OFFSET: int = -86

# TPP1 in 2I46: Crystal construct starts directly at Ser90 (label_seq_id = 1).
# auth_seq_num = label_seq_id + 89, matching canonical UniProt Q96AP0 (89-249).
TPP1_2I46_AUTH_OFFSET: int = 89
TPP1_2I46_LABEL_OFFSET: int = -89

# Canonical UniProt Accessions
TERT_UNIPROT_ACCESSION: str = "O14746"
TPP1_UNIPROT_ACCESSION: str = "Q96AP0"

# TEL-patch functional epitope residues (canonical UniProt numbering on TPP1)
TEL_PATCH_RESIDUES: tuple[int, ...] = (168, 169, 170, 171, 183, 212)
TEL_PATCH_DE_LOOP: tuple[int, ...] = (168, 169, 170, 171)
TEL_PATCH_EXTENDED: tuple[int, ...] = (168, 169, 170, 171, 180, 183, 211, 212, 215)

# TERT Structural Domain Boundaries (UniProt O14746 residue indices)
TERT_DOMAINS: dict[str, tuple[int, int]] = {
    "TEN": (1, 196),
    "LINKER": (197, 325),
    "TRBD": (326, 540),
    "IFD_TRAP": (755, 785),
    "RT": (541, 935),
    "CTE": (936, 1132),
}


@dataclass(frozen=True)
class ResidueInfo:
    """Detailed structural residue metadata."""

    resnum: int
    resname: str
    one_letter: str
    auth_seq_num: int
    label_seq_id: int | None
    auth_asym_id: str
    label_asym_id: str
    has_ca: bool


@dataclass
class AlignmentResult:
    """Pairwise sequence alignment and residue correspondence mapping."""

    seq1: str
    seq2: str
    mapping: dict[int, int]
    reverse_mapping: dict[int, int]
    aligned_pairs: list[tuple[int, int, str, str]]
    identity: float
    score: float
    unaligned_res1: list[int]
    unaligned_res2: list[int]
    alignment_str: str


# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------


def find_project_root() -> Path:
    """Locate the project root directory containing data_manifest.csv or src/."""
    current = Path(__file__).resolve().parent
    for parent in [current] + list(current.parents):
        if (parent / "data_manifest.csv").exists() or (parent / "src").exists():
            return parent
    return Path.cwd()


def get_tert_domain(uniprot_resnum: int) -> str:
    """
    Assign a TERT residue number (canonical UniProt O14746) to its structural domain.

    Domain boundaries:
      - TEN (Telomerase Essential N-terminal): residues 1-196 (modeled 7-178)
      - LINKER (Flexible disordered loop): residues 197-325 (unmodeled)
      - TRBD (RNA-binding domain): residues 326-540
      - IFD_TRAP (Insertion in Fingers / TPP1 cradle): residues 755-785 (nested inside RT)
      - RT (Reverse transcriptase palm/fingers): residues 541-935 (excluding IFD_TRAP)
      - CTE (C-terminal extension / thumb): residues 936-1132
    """
    if 755 <= uniprot_resnum <= 785:
        return "IFD_TRAP"
    if 1 <= uniprot_resnum <= 196:
        return "TEN"
    if 197 <= uniprot_resnum <= 325:
        return "LINKER"
    if 326 <= uniprot_resnum <= 540:
        return "TRBD"
    if 541 <= uniprot_resnum <= 935:
        return "RT"
    if 936 <= uniprot_resnum <= 1132:
        return "CTE"
    return "UNKNOWN"


# -----------------------------------------------------------------------------
# Sequence and Coordinate Extraction
# -----------------------------------------------------------------------------


ONE_TO_THREE_AA: dict[str, str] = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}


def one_to_three_aa(code: str) -> str:
    """Convert single-letter amino acid code to standard 3-letter uppercase code."""
    return ONE_TO_THREE_AA.get(code.upper(), f"AA{code.upper()}")


def _get_one_letter_code(resname: str) -> str:
    """Convert a 3-letter amino acid or nucleotide code to a single letter."""
    tab = gemmi.find_tabulated_residue(resname)
    if tab.is_amino_acid() or tab.is_nucleic_acid():
        code = tab.one_letter_code
        if code and code != "X":
            return code
    # Biopython / manual fallback
    try:
        c = seq1(resname)
        if c:
            return c
    except (KeyError, ValueError, TypeError) as err:
        logger.debug("Bio.SeqUtils.seq1 lookup failed for %s: %s", resname, err)
    std_map = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "MSE": "M", "SEC": "U", "PYL": "O",
        "A": "A", "C": "C", "G": "G", "U": "U",
        "DA": "A", "DC": "C", "DG": "G", "DT": "T",
    }
    return std_map.get(resname.upper(), "X")


def get_chain_sequence(
    cif_path: str | Path,
    chain_id: str,
    use_auth: bool = True,
    parser: str = "gemmi",
) -> tuple[str, list[int]]:
    """
    Extract the sequence of modeled polymer residues and their sequence numbers.

    Disambiguates mmCIF label_asym_id vs auth_asym_id:
    - If `use_auth=True`: looks up chain by auth_asym_id (with fallback to label_asym_id)
      and returns author residue numbers (`auth_seq_num`).
    - If `use_auth=False`: looks up chain by label_asym_id (with fallback to auth_asym_id)
      and returns CIF sequential residue numbers (`label_seq_id`).

    Args:
        cif_path: Path to the mmCIF coordinate file.
        chain_id: Identifier of the chain (e.g. 'A', 'B', 'C', 'D').
        use_auth: If True, use author numbering and chain ID; if False, use CIF label numbering.
        parser: Structural parser to use ('gemmi' or 'biopython').

    Returns:
        tuple[str, list[int]]: A 2-tuple containing:
            - sequence: Single-letter amino acid or nucleotide string.
            - residue_numbers: List of residue sequence numbers matching the characters.

    Raises:
        FileNotFoundError: If cif_path does not exist.
        KeyError: If chain_id is not found in the structure.
    """
    cif_path = Path(cif_path)
    if not cif_path.is_file():
        raise FileNotFoundError(f"mmCIF file not found: {cif_path}")

    if parser.lower() == "biopython":
        return _get_chain_sequence_biopython(cif_path, chain_id, use_auth=use_auth)
    return _get_chain_sequence_gemmi(cif_path, chain_id, use_auth=use_auth)


def _get_chain_sequence_gemmi(
    cif_path: Path,
    chain_id: str,
    use_auth: bool = True,
) -> tuple[str, list[int]]:
    """Gemmi implementation of get_chain_sequence."""
    st = gemmi.read_structure(str(cif_path))
    if not st:
        raise ValueError(f"Structure in {cif_path} contains no models")
    model = st[0]

    selected_residues: list[gemmi.Residue] = []

    if use_auth:
        # 1. Primary lookup by auth_asym_id (ch.name in Gemmi)
        for ch in model:
            if ch.name == chain_id:
                for r in ch:
                    tab = gemmi.find_tabulated_residue(r.name)
                    if tab.is_amino_acid() or tab.is_nucleic_acid():
                        # Disregard waters/heteroatoms that share author chain name
                        if tab.is_amino_acid() and r.subchain != ch.name and ch.name in ["A", "B"] and "2I46" in cif_path.name:
                            # 2I46 waters share auth chain name A/B but have subchain C/D
                            continue
                        selected_residues.append(r)
                break

        # Fallback to label_asym_id if not found by auth
        if not selected_residues:
            for ch in model:
                for r in ch:
                    if r.subchain == chain_id:
                        tab = gemmi.find_tabulated_residue(r.name)
                        if tab.is_amino_acid() or tab.is_nucleic_acid():
                            selected_residues.append(r)
    else:
        # 1. Primary lookup by label_asym_id (r.subchain in Gemmi)
        for ch in model:
            for r in ch:
                if r.subchain == chain_id:
                    tab = gemmi.find_tabulated_residue(r.name)
                    if tab.is_amino_acid() or tab.is_nucleic_acid():
                        selected_residues.append(r)

        # Fallback to auth_asym_id if not found by label
        if not selected_residues:
            for ch in model:
                if ch.name == chain_id:
                    for r in ch:
                        tab = gemmi.find_tabulated_residue(r.name)
                        if tab.is_amino_acid() or tab.is_nucleic_acid():
                            selected_residues.append(r)
                    break

    if not selected_residues:
        available_auth = [ch.name for ch in model]
        available_label = sorted({r.subchain for ch in model for r in ch})
        raise KeyError(
            f"Chain '{chain_id}' not found in {cif_path.name}. "
            f"Available auth chains: {available_auth}, label chains: {available_label}"
        )

    # Filter duplicate residue numbers if altlocs or microheterogeneity exist
    seen_nums: set[int] = set()
    seq_chars: list[str] = []
    res_nums: list[int] = []

    for r in selected_residues:
        num = r.seqid.num if use_auth else (r.label_seq if r.label_seq is not None else r.seqid.num)
        if num in seen_nums:
            continue
        seen_nums.add(num)
        seq_chars.append(_get_one_letter_code(r.name))
        res_nums.append(num)

    return "".join(seq_chars), res_nums


def _get_chain_sequence_biopython(
    cif_path: Path,
    chain_id: str,
    use_auth: bool = True,
) -> tuple[str, list[int]]:
    """Biopython implementation of get_chain_sequence."""
    if not use_auth:
        cif_dict = MMCIF2Dict(str(cif_path))
        label_asym_ids = cif_dict.get("_atom_site.label_asym_id", [])
        auth_asym_ids = cif_dict.get("_atom_site.auth_asym_id", [])
        label_seq_ids = cif_dict.get("_atom_site.label_seq_id", [])
        comp_ids = cif_dict.get("_atom_site.label_comp_id", [])

        # Match chains: primary lookup by label_asym_id, fallback to auth_asym_id
        target_asym_id: str | None = None
        use_field = "label"
        if chain_id in label_asym_ids:
            target_asym_id = chain_id
            use_field = "label"
        elif chain_id in auth_asym_ids:
            target_asym_id = chain_id
            use_field = "auth"
        else:
            raise KeyError(f"Chain '{chain_id}' not found in {cif_path.name} via Biopython")

        seen_nums: set[int] = set()
        seq_chars: list[str] = []
        res_nums: list[int] = []

        for la, aa, ls, comp in zip(label_asym_ids, auth_asym_ids, label_seq_ids, comp_ids):
            curr_asym = la if use_field == "label" else aa
            if curr_asym != target_asym_id:
                continue
            if ls in (".", "?", None):
                continue
            tab = gemmi.find_tabulated_residue(comp)
            if not (tab.is_amino_acid() or tab.is_nucleic_acid()):
                continue
            num = int(ls)
            if num in seen_nums:
                continue
            seen_nums.add(num)
            seq_chars.append(_get_one_letter_code(comp))
            res_nums.append(num)

        if not res_nums:
            raise KeyError(
                f"Chain '{chain_id}' not found or has no polymer residues in {cif_path.name} via Biopython"
            )

        return "".join(seq_chars), res_nums

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = MMCIFParser(QUIET=True)
        st = parser.get_structure(cif_path.stem, str(cif_path))

    model = st[0]
    target_chain = None

    if chain_id in model:
        target_chain = model[chain_id]
    else:
        for ch in model:
            if ch.id == chain_id:
                target_chain = ch
                break

    if target_chain is None:
        raise KeyError(f"Chain '{chain_id}' not found in {cif_path.name} via Biopython")

    seq_chars: list[str] = []
    res_nums: list[int] = []
    seen_nums: set[int] = set()

    for r in target_chain:
        # Exclude hetero/water residues
        if r.id[0] != " ":
            continue
        tab = gemmi.find_tabulated_residue(r.resname)
        if not (tab.is_amino_acid() or tab.is_nucleic_acid()):
            continue
        num = r.id[1]
        if num in seen_nums:
            continue
        seen_nums.add(num)
        seq_chars.append(_get_one_letter_code(r.resname))
        res_nums.append(num)

    return "".join(seq_chars), res_nums


def get_polymer_entity_sequence(
    cif_path: str | Path,
    entity_id_or_chain: str | int = 1,
) -> str:
    """
    Extract the canonical full polymer sequence from mmCIF entity definitions.

    Reads from `_entity_poly.pdbx_seq_one_letter_code`.

    Args:
        cif_path: Path to the mmCIF coordinate file.
        entity_id_or_chain: Entity numerical ID (e.g. '1', 2) or chain identifier.

    Returns:
        str: Clean single-letter sequence string without whitespace or formatting characters.
    """
    cif_path = Path(cif_path)
    doc = gemmi.cif.read_file(str(cif_path))
    block = doc.sole_block()

    entity_str = str(entity_id_or_chain).strip()

    # If chain identifier was given instead of entity id, map chain to entity
    table_asym = block.find("_struct_asym.", ["id", "entity_id"])
    if table_asym:
        for row in table_asym:
            if row[0] == entity_str:
                entity_str = row[1]
                break

    table_poly = block.find(
        "_entity_poly.",
        ["entity_id", "pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can"],
    )
    if table_poly:
        for row in table_poly:
            if row[0] == entity_str:
                seq = row[1].strip(";\n ").replace("\n", "").replace(" ", "")
                return seq

    # Fallback to first polymer entity if matching entity not found
    if table_poly and len(table_poly) > 0:
        return table_poly[0][1].strip(";\n ").replace("\n", "").replace(" ", "")

    return ""


def get_modeled_residues(
    cif_path: str | Path,
    chain_id: str,
    use_auth: bool = True,
) -> list[ResidueInfo]:
    """
    Extract all modeled residues for a chain with complete metadata.

    Args:
        cif_path: Path to mmCIF coordinate file.
        chain_id: Chain identifier.
        use_auth: True to select by auth_asym_id; False for label_asym_id.

    Returns:
        list[ResidueInfo]: Modeled residues with coordinates and offsets.
    """
    cif_path = Path(cif_path)
    st = gemmi.read_structure(str(cif_path))
    model = st[0]

    residues: list[ResidueInfo] = []
    seen: set[int] = set()

    for ch in model:
        if (use_auth and ch.name == chain_id) or (not use_auth and any(r.subchain == chain_id for r in ch)):
            for r in ch:
                if not use_auth and r.subchain != chain_id:
                    continue
                tab = gemmi.find_tabulated_residue(r.name)
                if not (tab.is_amino_acid() or tab.is_nucleic_acid()):
                    continue
                auth_num = r.seqid.num
                if auth_num in seen:
                    continue
                seen.add(auth_num)
                has_ca = r.find_atom("CA", "*") is not None
                residues.append(
                    ResidueInfo(
                        resnum=auth_num if use_auth else (r.label_seq or auth_num),
                        resname=r.name,
                        one_letter=_get_one_letter_code(r.name),
                        auth_seq_num=auth_num,
                        label_seq_id=r.label_seq,
                        auth_asym_id=ch.name,
                        label_asym_id=r.subchain,
                        has_ca=has_ca,
                    )
                )
            break

    return residues


# -----------------------------------------------------------------------------
# Construct Tag Offset Reconcilers
# -----------------------------------------------------------------------------


def reconcile_tert_offset(label_seq_id: int) -> int:
    """
    Reconcile TERT CIF label_seq_id to canonical UniProt O14746 sequence number.

    The recombinant TERT expression construct in 7TRE, 7TRD, and 7TRF contains an
    N-terminal 35-residue Twin-Strep purification tag:
        auth_seq_num = label_seq_id - 35

    Example:
        label_seq_id = 42 (Cys) -> auth_seq_num = 7 (canonical UniProt O14746)
        label_seq_id = 1167 (Asp) -> auth_seq_num = 1132
    """
    return label_seq_id + TERT_AUTH_OFFSET


def reconcile_tpp1_offset(label_seq_id: int, structure_id: str = "7TRE") -> int:
    """
    Reconcile TPP1 CIF label_seq_id to canonical UniProt Q96AP0 sequence number.

    Offsets:
    - 7TRE: Construct contains 3-aa N-terminal GAG cloning tag.
      auth_seq_num = label_seq_id + 86.
      Example: label_seq_id = 82 (Glu168) -> auth_seq_num = 168.
    - 2I46: Construct begins directly at Ser90 (label_seq_id = 1).
      auth_seq_num = label_seq_id + 89.
      Example: label_seq_id = 79 (Glu168) -> auth_seq_num = 168.
    """
    pid = structure_id.upper()
    if "7TRE" in pid:
        return label_seq_id + TPP1_7TRE_AUTH_OFFSET
    if "2I46" in pid:
        return label_seq_id + TPP1_2I46_AUTH_OFFSET
    # Default to 7TRE convention if unknown
    return label_seq_id + TPP1_7TRE_AUTH_OFFSET


def reconcile_label_from_auth(
    auth_seq_num: int,
    protein: str = "TERT",
    structure_id: str = "7TRE",
) -> int:
    """
    Convert an author / UniProt sequence number to CIF label_seq_id.

    Args:
        auth_seq_num: Author residue index.
        protein: "TERT" or "TPP1".
        structure_id: "7TRE", "7TRD", "7TRF", or "2I46".

    Returns:
        int: CIF label_seq_id.
    """
    prot = protein.upper()
    pid = structure_id.upper()
    if "TERT" in prot:
        return auth_seq_num + TERT_LABEL_OFFSET
    if "TPP1" in prot or "ACD" in prot:
        if "2I46" in pid:
            return auth_seq_num + TPP1_2I46_LABEL_OFFSET
        return auth_seq_num + TPP1_7TRE_LABEL_OFFSET
    return auth_seq_num


def get_uniprot_mapping(
    cif_path: str | Path,
    chain_id: str,
    use_auth: bool = True,
) -> dict[int, int]:
    """
    Map structure residue numbers to canonical UniProt sequence numbers.

    For TERT (7TRE, 7TRD, 7TRF):
    - When `use_auth=True`: author residue numbers match UniProt O14746 1:1 (7-1132).
    - When `use_auth=False`: `label_seq_id` maps to `label_seq_id - 35`.

    For TPP1 (7TRE, 2I46):
    - When `use_auth=True`: author residue numbers match UniProt Q96AP0 1:1 (90-241).
    - When `use_auth=False`:
        - in 7TRE: `label_seq_id` maps to `label_seq_id + 86`.
        - in 2I46: `label_seq_id` maps to `label_seq_id + 89`.

    Non-protein entities (e.g. hTR RNA, telomeric DNA) or non-TERT/TPP1 chains
    return an empty dictionary.

    Args:
        cif_path: Path to the mmCIF coordinate file.
        chain_id: Chain identifier.
        use_auth: If True, input keys are author numbers; if False, CIF label numbers.

    Returns:
        dict[int, int]: Mapping from structure residue number to UniProt sequence index.
    """
    cif_path = Path(cif_path)
    pid = cif_path.stem.upper()
    try:
        _, res_nums = get_chain_sequence(cif_path, chain_id, use_auth=use_auth)
    except KeyError:
        return {}

    if not res_nums:
        return {}

    # Disambiguate entity type and macromolecule via mmCIF block
    st = gemmi.read_structure(str(cif_path))
    model = st[0]
    target_label_asym: str | None = None

    if use_auth:
        for ch in model:
            if ch.name == chain_id:
                for r in ch:
                    target_label_asym = r.subchain
                    break
                break
        if target_label_asym is None:
            for ch in model:
                for r in ch:
                    if r.subchain == chain_id:
                        target_label_asym = r.subchain
                        break
                if target_label_asym:
                    break
    else:
        for ch in model:
            for r in ch:
                if r.subchain == chain_id:
                    target_label_asym = chain_id
                    break
            if target_label_asym:
                break
        if target_label_asym is None:
            for ch in model:
                if ch.name == chain_id:
                    for r in ch:
                        target_label_asym = r.subchain
                        break
                    break

    if not target_label_asym:
        return {}

    doc = gemmi.cif.read_file(str(cif_path))
    block = doc.sole_block()

    ent_id: str | None = None
    table_asym = block.find("_struct_asym.", ["id", "entity_id"])
    if table_asym:
        for row in table_asym:
            if row[0] == target_label_asym:
                ent_id = row[1]
                break

    if not ent_id:
        return {}

    # Check polymer type in _entity_poly
    poly_type = None
    table_poly = block.find("_entity_poly.", ["entity_id", "type"])
    if table_poly:
        for row in table_poly:
            if row[0] == ent_id:
                poly_type = row[1].lower()
                break

    # If nucleic acid or non-protein, return empty dict
    if not poly_type or "polypeptide" not in poly_type:
        return {}

    desc = ""
    table_ent = block.find("_entity.", ["id", "pdbx_description"])
    if table_ent:
        for row in table_ent:
            if row[0] == ent_id:
                desc = row[1].lower()
                break

    is_tert = "reverse transcriptase" in desc or "tert" in desc
    is_tpp1 = "adrenocortical dysplasia" in desc or "tpp1" in desc or "acd" in desc

    if not (is_tert or is_tpp1):
        return {}

    mapping: dict[int, int] = {}
    for num in res_nums:
        if is_tert:
            if use_auth:
                mapping[num] = num
            else:
                mapping[num] = reconcile_tert_offset(num)
        elif is_tpp1:
            if use_auth:
                mapping[num] = num
            else:
                mapping[num] = reconcile_tpp1_offset(num, structure_id=pid)

    return mapping


# -----------------------------------------------------------------------------
# 2I46 Dimer Asymmetry Handling
# -----------------------------------------------------------------------------


def check_tel_patch_integrity(
    cif_path: str | Path,
    chain_id: str = "A",
) -> dict[str, Any]:
    """
    Check whether the TEL-patch recruitment epitope is resolved in a TPP1 chain.

    Specifically inspects the critical acidic DE-loop (Glu168, Glu169, Lys170, Glu171)
    and hydrophobic core residues (Leu183, Leu212).

    In crystal structure 2I46:
    - Chain A has the TEL-patch DE-loop fully resolved.
    - Chain B has the DE-loop disordered/unobserved (residues 168-171 missing coordinates).

    Args:
        cif_path: Path to 2I46 or 7TRE mmCIF file.
        chain_id: Chain identifier to evaluate ('A' or 'B' for 2I46; 'C' or 'D' for 7TRE).

    Returns:
        dict[str, Any]: Detailed integrity assessment containing:
            - chain_id: Chain evaluated.
            - resolved: True if all TEL-patch residues are modeled.
            - modeled_residues: List of modeled TEL-patch residues.
            - missing_residues: List of unmodeled TEL-patch residues.
            - de_loop_resolved: True if residues 168-171 are all modeled.
            - warning: Warning message if disordered or None.
    """
    cif_path = Path(cif_path)
    _, res_nums = get_chain_sequence(cif_path, chain_id, use_auth=True)
    res_set = set(res_nums)

    modeled = [r for r in TEL_PATCH_RESIDUES if r in res_set]
    missing = [r for r in TEL_PATCH_RESIDUES if r not in res_set]
    de_loop_modeled = all(r in res_set for r in TEL_PATCH_DE_LOOP)

    resolved = len(missing) == 0
    warning_msg = None

    if not resolved:
        warning_msg = (
            f"TEL-patch residue disorder detected in {cif_path.name} chain '{chain_id}': "
            f"Missing residues {missing}. DE-loop resolved: {de_loop_modeled}."
        )
        logger.warning(warning_msg)
        warnings.warn(warning_msg, UserWarning, stacklevel=2)

    return {
        "structure_id": cif_path.stem.upper(),
        "chain_id": chain_id,
        "resolved": resolved,
        "modeled_residues": modeled,
        "missing_residues": missing,
        "de_loop_resolved": de_loop_modeled,
        "warning": warning_msg,
    }


def select_valid_tpp1_chain(
    cif_path: str | Path,
    preferred_chain: str = "A",
    require_tel_patch: bool = True,
) -> str:
    """
    Select a valid TPP1 chain that has the TEL-patch recruitment epitope resolved.

    If 2I46 Chain B is requested with `require_tel_patch=True`, automatically redirects
    to Chain A (where the DE-loop 168-171 is fully resolved) and issues a warning.

    Args:
        cif_path: Path to mmCIF file (especially 2I46).
        preferred_chain: The requested chain ID ('A' or 'B').
        require_tel_patch: If True, enforce that TEL-patch residues must be modeled.

    Returns:
        str: Selected chain ID guaranteed to have resolved TEL-patch coordinates.
    """
    cif_path = Path(cif_path)
    pid = cif_path.stem.upper()

    if "2I46" in pid:
        if preferred_chain.upper() == "B":
            if require_tel_patch:
                warning_msg = (
                    "2I46 Chain B has disordered TEL-patch DE-loop (residues 168-171 unmodeled). "
                    "Redirecting to 2I46 Chain A which has fully modeled TEL-patch coordinates."
                )
                logger.warning(warning_msg)
                warnings.warn(warning_msg, UserWarning, stacklevel=2)
                return "A"
            return "B"
        return "A"

    if "7TRE" in pid:
        # In 7TRE: author chain is C, label chain is D
        if preferred_chain.upper() in ["C", "D"]:
            return preferred_chain
        return "C"

    return preferred_chain


# -----------------------------------------------------------------------------
# Sequence-Aware Homology Alignment
# -----------------------------------------------------------------------------


def align_sequences(
    seq1: str,
    seq2: str,
    mode: str = "global",
) -> AlignmentResult:
    """
    Perform sequence alignment between two polymer sequences with affine gap penalties.

    Uses Biopython Bio.Align.PairwiseAligner (Needleman-Wunsch / Smith-Waterman).

    Args:
        seq1: First sequence string.
        seq2: Second sequence string.
        mode: Alignment mode ('global' or 'local').

    Returns:
        AlignmentResult: Alignment object with matched coordinates and statistics.
    """
    if not seq1 or not seq2:
        return AlignmentResult(
            seq1=seq1,
            seq2=seq2,
            mapping={},
            reverse_mapping={},
            aligned_pairs=[],
            identity=0.0,
            score=0.0,
            unaligned_res1=list(range(len(seq1))),
            unaligned_res2=list(range(len(seq2))),
            alignment_str="",
        )

    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return AlignmentResult(
            seq1=seq1,
            seq2=seq2,
            mapping={},
            reverse_mapping={},
            aligned_pairs=[],
            identity=0.0,
            score=0.0,
            unaligned_res1=list(range(len(seq1))),
            unaligned_res2=list(range(len(seq2))),
            alignment_str="",
        )

    best = alignments[0]
    mapping: dict[int, int] = {}
    reverse_mapping: dict[int, int] = {}
    aligned_pairs: list[tuple[int, int, str, str]] = []
    matched_chars = 0

    aligned_indices_1: set[int] = set()
    aligned_indices_2: set[int] = set()

    for (t_start, t_end), (q_start, q_end) in zip(best.aligned[0], best.aligned[1]):
        block_len = t_end - t_start
        for k in range(block_len):
            idx1 = t_start + k
            idx2 = q_start + k
            c1 = seq1[idx1]
            c2 = seq2[idx2]
            mapping[idx1] = idx2
            reverse_mapping[idx2] = idx1
            aligned_pairs.append((idx1, idx2, c1, c2))
            aligned_indices_1.add(idx1)
            aligned_indices_2.add(idx2)
            if c1 == c2:
                matched_chars += 1

    total_len = max(len(seq1), len(seq2), 1)
    identity = matched_chars / total_len

    unaligned_1 = [i for i in range(len(seq1)) if i not in aligned_indices_1]
    unaligned_2 = [i for i in range(len(seq2)) if i not in aligned_indices_2]

    return AlignmentResult(
        seq1=seq1,
        seq2=seq2,
        mapping=mapping,
        reverse_mapping=reverse_mapping,
        aligned_pairs=aligned_pairs,
        identity=identity,
        score=float(best.score),
        unaligned_res1=unaligned_1,
        unaligned_res2=unaligned_2,
        alignment_str=format(best),
    )


def map_residue_correspondence(
    cif_path_1: str | Path,
    chain_1: str,
    cif_path_2: str | Path,
    chain_2: str,
    use_auth: bool = True,
) -> dict[int, int]:
    """
    Establish sequence-aware residue correspondence between homologous chains across structures.

    Aligns modeled sequences from two structures using global pairwise alignment
    before any geometric superposition. Ensures nominal author numbering discrepancies
    do not misalign coordinates.

    Args:
        cif_path_1: Path to first mmCIF file (reference).
        chain_1: Chain identifier in first structure.
        cif_path_2: Path to second mmCIF file (target).
        chain_2: Chain identifier in second structure.
        use_auth: If True, map between author residue numbers; if False, CIF label numbers.

    Returns:
        dict[int, int]: Mapping from residue number in chain_1 to equivalent residue
            number in chain_2.
    """
    result = align_homologous_chains(
        cif_path_1=cif_path_1,
        chain_1=chain_1,
        cif_path_2=cif_path_2,
        chain_2=chain_2,
        use_auth=use_auth,
    )
    return result.mapping


def align_homologous_chains(
    cif_path_1: str | Path,
    chain_1: str,
    cif_path_2: str | Path,
    chain_2: str,
    use_auth: bool = True,
) -> AlignmentResult:
    """
    Extract modeled sequences and perform complete pairwise residue alignment.

    Args:
        cif_path_1: Path to reference mmCIF file.
        chain_1: Reference chain identifier.
        cif_path_2: Path to target mmCIF file.
        chain_2: Target chain identifier.
        use_auth: If True, operate on author numbers; if False, CIF label numbers.

    Returns:
        AlignmentResult: Alignment object with residue number mappings and metrics.
    """
    cif_path_1 = Path(cif_path_1)
    cif_path_2 = Path(cif_path_2)

    seq1, nums1 = get_chain_sequence(cif_path_1, chain_1, use_auth=use_auth)
    seq2, nums2 = get_chain_sequence(cif_path_2, chain_2, use_auth=use_auth)

    aln_result = align_sequences(seq1, seq2, mode="global")

    # Re-index mapping from string 0-based indices to real structure residue numbers
    res_mapping: dict[int, int] = {}
    res_rev_mapping: dict[int, int] = {}
    res_aligned_pairs: list[tuple[int, int, str, str]] = []

    for idx1, idx2, c1, c2 in aln_result.aligned_pairs:
        r1 = nums1[idx1]
        r2 = nums2[idx2]
        res_mapping[r1] = r2
        res_rev_mapping[r2] = r1
        res_aligned_pairs.append((r1, r2, c1, c2))

    unaligned_res1 = [nums1[i] for i in aln_result.unaligned_res1]
    unaligned_res2 = [nums2[i] for i in aln_result.unaligned_res2]

    return AlignmentResult(
        seq1=seq1,
        seq2=seq2,
        mapping=res_mapping,
        reverse_mapping=res_rev_mapping,
        aligned_pairs=res_aligned_pairs,
        identity=aln_result.identity,
        score=aln_result.score,
        unaligned_res1=unaligned_res1,
        unaligned_res2=unaligned_res2,
        alignment_str=aln_result.alignment_str,
    )


def map_tert_correspondence(
    cif_ref: str | Path,
    cif_target: str | Path,
    ref_chain: str = "A",
    target_chain: str = "A",
    use_auth: bool = True,
) -> dict[int, int]:
    """
    Convenience wrapper to map equivalent TERT residues between reference and target.

    Supports cross-mapping between 7TRE, 7TRD, and 7TRF.
    """
    return map_residue_correspondence(
        cif_path_1=cif_ref,
        chain_1=ref_chain,
        cif_path_2=cif_target,
        chain_2=target_chain,
        use_auth=use_auth,
    )


def map_tpp1_correspondence(
    cif_7tre: str | Path,
    cif_2i46: str | Path,
    chain_7tre: str = "C",
    chain_2i46: str = "A",
    use_auth: bool = True,
) -> dict[int, int]:
    """
    Convenience wrapper to map equivalent TPP1 residues between 7TRE and 2I46 Chain A.

    Enforces selection of 2I46 Chain A to prevent unmodeled DE-loop errors in Chain B.
    """
    valid_2i46_chain = select_valid_tpp1_chain(cif_2i46, preferred_chain=chain_2i46, require_tel_patch=True)
    return map_residue_correspondence(
        cif_path_1=cif_7tre,
        chain_1=chain_7tre,
        cif_path_2=cif_2i46,
        chain_2=valid_2i46_chain,
        use_auth=use_auth,
    )


# -----------------------------------------------------------------------------
# Standard Results Table Generators
# -----------------------------------------------------------------------------


def generate_tert_correspondence_table(
    raw_dir: str | Path | None = None,
) -> pd.DataFrame:
    """
    Generate residue_correspondence_tert.csv across 7TRE, 7TRD, and 7TRF.

    Schema:
        uniprot_resnum, resname_canonical, label_seq_id_7tre, auth_seq_id_7tre,
        modeled_in_7tre, label_seq_id_7trd, auth_seq_id_7trd, modeled_in_7trd,
        label_seq_id_7trf, auth_seq_id_7trf, modeled_in_7trf, domain_assignment

    Returns:
        pd.DataFrame: Table with exactly 1132 rows (canonical UniProt O14746).
    """
    if raw_dir is None:
        raw_dir = find_project_root() / "data" / "raw"
    raw_dir = Path(raw_dir)

    cif_7tre = raw_dir / "7TRE.cif"
    cif_7trd = raw_dir / "7TRD.cif"
    cif_7trf = raw_dir / "7TRF.cif"

    # Extract modeled residue sets (author residue numbers)
    _, nums_tre = get_chain_sequence(cif_7tre, "A", use_auth=True)
    _, nums_trd = get_chain_sequence(cif_7trd, "A", use_auth=True)
    _, nums_trf = get_chain_sequence(cif_7trf, "A", use_auth=True)

    set_tre = set(nums_tre)
    set_trd = set(nums_trd)
    set_trf = set(nums_trf)

    # Canonical full sequence from 7TRE Entity 2
    full_poly = get_polymer_entity_sequence(cif_7tre, entity_id_or_chain=2)
    # Exclude the 35-aa Twin-Strep tag to get UniProt O14746
    uniprot_seq = full_poly[35:] if len(full_poly) >= 1167 else full_poly

    rows: list[dict[str, Any]] = []

    for resnum in range(1, 1133):
        resname_1 = uniprot_seq[resnum - 1] if resnum - 1 < len(uniprot_seq) else "X"
        resname_3 = one_to_three_aa(resname_1)

        mod_tre = resnum in set_tre
        mod_trd = resnum in set_trd
        mod_trf = resnum in set_trf

        label_tre = resnum + TERT_LABEL_OFFSET if mod_tre else ""
        auth_tre = resnum if mod_tre else ""
        label_trd = resnum + TERT_LABEL_OFFSET if mod_trd else ""
        auth_trd = resnum if mod_trd else ""
        label_trf = resnum + TERT_LABEL_OFFSET if mod_trf else ""
        auth_trf = resnum if mod_trf else ""

        rows.append(
            {
                "uniprot_resnum": resnum,
                "resname_canonical": resname_3,
                "label_seq_id_7tre": label_tre,
                "auth_seq_id_7tre": auth_tre,
                "modeled_in_7tre": mod_tre,
                "label_seq_id_7trd": label_trd,
                "auth_seq_id_7trd": auth_trd,
                "modeled_in_7trd": mod_trd,
                "label_seq_id_7trf": label_trf,
                "auth_seq_id_7trf": auth_trf,
                "modeled_in_7trf": mod_trf,
                "domain_assignment": get_tert_domain(resnum),
            }
        )

    return pd.DataFrame(rows)


def generate_tpp1_correspondence_table(
    raw_dir: str | Path | None = None,
) -> pd.DataFrame:
    """
    Generate residue_correspondence_tpp1.csv across 7TRE, 2I46 Chain A, and 2I46 Chain B.

    Schema:
        uniprot_resnum, resname_canonical, label_seq_id_7tre, auth_seq_id_7tre,
        modeled_in_7tre, label_seq_id_2i46_chain_a, auth_seq_id_2i46_chain_a,
        modeled_in_2i46_chain_a, modeled_in_2i46_chain_b, is_tel_patch_residue

    Returns:
        pd.DataFrame: Table for TPP1 OB-fold domain (UniProt Q96AP0 residues 89-249).
    """
    if raw_dir is None:
        raw_dir = find_project_root() / "data" / "raw"
    raw_dir = Path(raw_dir)

    cif_7tre = raw_dir / "7TRE.cif"
    cif_2i46 = raw_dir / "2I46.cif"

    _, nums_tre = get_chain_sequence(cif_7tre, "C", use_auth=True)
    _, nums_i46_a = get_chain_sequence(cif_2i46, "A", use_auth=True)
    _, nums_i46_b = get_chain_sequence(cif_2i46, "B", use_auth=True)

    set_tre = set(nums_tre)
    set_i46_a = set(nums_i46_a)
    set_i46_b = set(nums_i46_b)

    # Canonical TPP1 OB-fold sequence from 7TRE Entity 4
    # Recombinant construct has a 3-aa N-terminal tag (Gly87-Ala88-Gly89 in entity numbering).
    # Stripping the 3-aa tag yields ob_seq starting at residue 90 (Ser90 = ob_seq[0]).
    # UniProt Q96AP0 residue 89 corresponds to full_poly[2] (Gly89).
    full_poly = get_polymer_entity_sequence(cif_7tre, entity_id_or_chain=4)
    ob_seq = full_poly[3:] if len(full_poly) >= 163 else full_poly

    rows: list[dict[str, Any]] = []

    # UniProt Q96AP0 residues 89 to 249 (161 positions)
    for resnum in range(89, 250):
        if resnum == 89:
            resname_1 = full_poly[2] if len(full_poly) >= 3 else "G"
        else:
            idx = resnum - 90
            resname_1 = ob_seq[idx] if 0 <= idx < len(ob_seq) else "X"
        resname_3 = one_to_three_aa(resname_1)

        mod_tre = resnum in set_tre
        mod_i46_a = resnum in set_i46_a
        mod_i46_b = resnum in set_i46_b

        label_tre = resnum + TPP1_7TRE_LABEL_OFFSET if mod_tre else ""
        auth_tre = resnum if mod_tre else ""
        label_i46_a = resnum + TPP1_2I46_LABEL_OFFSET if mod_i46_a else ""
        auth_i46_a = resnum if mod_i46_a else ""

        is_tel = resnum in TEL_PATCH_RESIDUES

        rows.append(
            {
                "uniprot_resnum": resnum,
                "resname_canonical": resname_3,
                "label_seq_id_7tre": label_tre,
                "auth_seq_id_7tre": auth_tre,
                "modeled_in_7tre": mod_tre,
                "label_seq_id_2i46_chain_a": label_i46_a,
                "auth_seq_id_2i46_chain_a": auth_i46_a,
                "modeled_in_2i46_chain_a": mod_i46_a,
                "modeled_in_2i46_chain_b": mod_i46_b,
                "is_tel_patch_residue": is_tel,
            }
        )

    return pd.DataFrame(rows)


def generate_chain_entity_mapping(
    raw_dir: str | Path | None = None,
) -> pd.DataFrame:
    """
    Generate chain_entity_mapping.csv across all 5 project structures.

    Schema:
        pdb_id, entity_id, entity_type, macromolecule_name, label_asym_id,
        auth_asym_id, seq_len_canonical, uniprot_accession, modeled_residues_count,
        modeled_range_auth, tag_or_construct_notes

    Returns:
        pd.DataFrame: Comprehensive chain and entity reference mapping.
    """
    if raw_dir is None:
        raw_dir = find_project_root() / "data" / "raw"
    raw_dir = Path(raw_dir)

    target_ids = ["7TRE", "7TRD", "7TRF", "7TRC", "2I46"]
    rows: list[dict[str, Any]] = []

    for pid in target_ids:
        cif_path = raw_dir / f"{pid}.cif"
        if not cif_path.is_file():
            raise FileNotFoundError(f"Structure coordinate file not found: {cif_path}")

        doc = gemmi.cif.read_file(str(cif_path))
        b = doc.sole_block()
        st = gemmi.read_structure(str(cif_path))
        model = st[0]

        sub_info: dict[str, dict[str, Any]] = {}
        for ch in model:
            for r in ch:
                tab = gemmi.find_tabulated_residue(r.name)
                if tab.is_amino_acid() or tab.is_nucleic_acid():
                    sub = r.subchain
                    if sub not in sub_info:
                        sub_info[sub] = {"auth": ch.name, "nums": set()}
                    sub_info[sub]["nums"].add(r.seqid.num)

        ent_table = b.find("_entity.", ["id", "type", "pdbx_description"])
        desc_map: dict[str, str] = {}
        type_map: dict[str, str] = {}
        for r in ent_table:
            if r[1] == "polymer":
                desc_map[r[0]] = r[2].strip("'\"")
                type_map[r[0]] = r[1]

        poly_table = b.find("_entity_poly.", ["entity_id", "pdbx_seq_one_letter_code"])
        len_map: dict[str, int] = {}
        if poly_table:
            for r in poly_table:
                seq = r[1].strip(";\n ").replace("\n", "").replace(" ", "")
                len_map[r[0]] = len(seq)

        ref_table = b.find("_struct_ref.", ["entity_id", "db_name", "pdbx_db_accession"])
        unp_map: dict[str, str] = {}
        if ref_table:
            for r in ref_table:
                if r[1] == "UNP":
                    unp_map[r[0]] = r[2]

        asym_table = b.find("_struct_asym.", ["id", "entity_id"])
        ent_to_asym: dict[str, list[str]] = {}
        for r in asym_table:
            eid = r[1]
            if eid in type_map:
                ent_to_asym.setdefault(eid, []).append(r[0])

        sorted_eids = sorted(type_map.keys(), key=lambda x: int(x) if x.isdigit() else x)
        for eid in sorted_eids:
            label_asyms = ent_to_asym.get(eid, [])
            auth_asyms: list[str] = []
            tot_modeled = 0
            ranges: list[str] = []
            for la in label_asyms:
                if la in sub_info:
                    auth = sub_info[la]["auth"]
                    if auth not in auth_asyms:
                        auth_asyms.append(auth)
                    r_nums = sorted(sub_info[la]["nums"])
                    tot_modeled += len(r_nums)
                    if r_nums:
                        cur_start = r_nums[0]
                        cur_prev = r_nums[0]
                        sub_ranges: list[str] = []
                        for n in r_nums[1:]:
                            if n == cur_prev + 1:
                                cur_prev = n
                            else:
                                sub_ranges.append(
                                    f"{cur_start}-{cur_prev}" if cur_start != cur_prev else f"{cur_start}"
                                )
                                cur_start = n
                                cur_prev = n
                        sub_ranges.append(
                            f"{cur_start}-{cur_prev}" if cur_start != cur_prev else f"{cur_start}"
                        )
                        range_str = ", ".join(sub_ranges)
                        ranges.append(
                            f"{auth}: {range_str}" if len(label_asyms) > 1 else range_str
                        )

            desc = desc_map.get(eid, "")
            unp = unp_map.get(eid, "N/A")

            notes = ""
            if "reverse transcriptase" in desc.lower():
                notes = "N-terminal 35-aa Twin-Strep tag; auth = label - 35"
            elif "adrenocortical dysplasia" in desc.lower():
                if pid == "7TRE":
                    notes = "N-terminal 3-aa GAG tag; auth = label + 86; TEL-patch resolved"
                elif pid == "2I46":
                    notes = "Homodimer in ASU; Chain A resolved TEL patch; Chain B DE-loop disordered"
            elif "telomerase rna" in desc.lower():
                notes = (
                    "Pseudoknot & CR4/CR5; template 5'-CUAACCCUAAC-3'"
                    if pid == "7TRE"
                    else "hTR RNA transcript component"
                )
            elif "telomeric repeat" in desc.lower():
                notes = "Telomeric DNA substrate"
            elif "histone" in desc.lower():
                notes = "Nucleosomal histone component"
            else:
                notes = "H/ACA complex component" if pid == "7TRC" else "Polymer component"

            rows.append(
                {
                    "pdb_id": pid,
                    "entity_id": eid,
                    "entity_type": "polymer",
                    "macromolecule_name": desc,
                    "label_asym_id": ", ".join(label_asyms),
                    "auth_asym_id": ", ".join(auth_asyms),
                    "seq_len_canonical": len_map.get(eid, 0),
                    "uniprot_accession": unp,
                    "modeled_residues_count": tot_modeled,
                    "modeled_range_auth": "; ".join(ranges) if ranges else "None",
                    "tag_or_construct_notes": notes,
                }
            )

    return pd.DataFrame(rows)
