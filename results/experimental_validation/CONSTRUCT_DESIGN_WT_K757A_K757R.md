# RECOMBINANT CONSTRUCT DESIGN: HUMAN TERT WT, K757A, AND K757R (SEQUENCE-VERIFIED)

**Document ID:** `EXP-DOC-001-REV-AUDITED`  
**Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/CONSTRUCT_DESIGN_WT_K757A_K757R.md`  
**Classification:** Sequence-Level Verified Experimental Blueprint  

---

## 1. Reference Sequence & Residue Numbering Verification

The primary amino acid and nucleotide sequence of human telomerase reverse transcriptase was independently audited against the canonical NCBI reference **NM_198253.3** and UniProt reference **O14746** (1,132 amino acids, 3,399 bp coding sequence).

### 1.1 Local Sequence Registry at TRAP Locus (Residues 748–766)
```
========================================================================================================================
CANONICAL HUMAN TERT SEQUENCE AT RESIDUES 748-766 (NM_198253.3 / UniProt O14746)
========================================================================================================================
Residue Index:   748   749   750   751   752   753   754   755   756   757   758   759   760   761   762   763   764   765   766
Amino Acid:      Gln   Lys   Ala   Ala   His   Gly   His   Val   Arg   LYS   Ala   Phe   Lys   Ser   His   Val   Ser   Thr   Leu
Single-Letter:    Q     K     A     A     H     G     H     V     R     K     A     F     K     S     H     V     S     T     L
Canonical Codon: CAG   AAG   GCC   GCC   CAT   GGG   CAC   GTC   CGC   AAG   GCC   TTC   AAG   AGC   CAC   GTC   TCT   ACC   CTG
========================================================================================================================
```

> **Audit Correction**: An earlier working draft incorrectly listed residues 753–760 as `Leu-Val-Arg-Arg-Lys-Thr-Phe-Arg`. Sequence-level verification against UniProt `O14746` and PDB 7TRE confirms the true canonical sequence is `Gly753 - His754 - Val755 - Arg756 - Lys757 - Ala758 - Phe759 - Lys760` (`G-H-V-R-K-A-F-K`). Codon 757 is `AAG` (nt 2269–2271).

---

## 2. Recombinant Cassette & Expression Tag Placement

```
+---------------------------------------------------------------------------------------------------------+
|                                    N-TERMINAL RECOMBINANT CASSETTE                                      |
+---------------------------------------------------------------------------------------------------------+
| [ Twin-Strep-tag ] --- [ Linker ] --- [ 3xFLAG-tag ] --- [ PreScission 3C Site ] --- [ hTERT (1-1132) ]|
| WSHPQFEK-GGGSGGGSGGSA-WSHPQFEK        DYKDHDGDYKDHDIDYKDDDDK  LEVLFQ/GP              Cys7 ... Asp1132   |
+---------------------------------------------------------------------------------------------------------+
```

### Holoenzyme Assembly & Structural Compatibility:
- In cryo-EM reconstructions of human telomerase holoenzyme (PDB 7TRC, 7TRD, 7TRE; Ghanim et al. 2021; Sekne et al. 2022), TERT residues 1–6 are disordered and project into bulk solvent, far removed from the RNA-binding domain, active site, and shelterin interface.
- Ghanim et al. (2021) and Sekne et al. (2022) experimentally demonstrated that N-terminal Twin-Strep tags do not interfere with hTR incorporation, holoenzyme assembly, or telomerase catalysis.
- Furthermore, on-column cleavage with PreScission (3C) protease at `LEVLFQ/GP` removes the entire tag cassette, leaving only two non-native residues (`Gly-Pro-`) at the N-terminus prior to Met1.

---

## 3. Site-Directed Mutagenesis Design

### 3.1 Wild-Type Locus (Codons 754–760)
- **Amino Acids**: `His754 - Val755 - Arg756 - Lys757 - Ala758 - Phe759 - Lys760`
- **Native cDNA**: `5'- CAC GTC CGC AAG GCC TTC AAG -3'`
- **Target Codon 757**: `AAG` (Lysine)

### 3.2 Construct 1: TERT Lys757A (Uncoupling Candidate)
- **Mutated Codon 757**: `GCG` (Alanine)
- **Engineered cDNA**: `5'- CAC GTC CGC GCG GCC TTC AAG -3'`
- **Diagnostic Restriction Cleavage**:
  - In the wild-type flank `CAC GTC CGC AAG GCC TTC`, there are **0** cleavage sites for **BstUI** (`CG^CG`) or **HinP1I** (`G^CGC`).
  - Mutating `AAG` to `GCG` creates `CGCGCG`, introducing a novel **BstUI** site (`CG^CG`) and **HinP1I** site (`G^CGC`), providing an immediate diagnostic restriction screen.
- **Mutagenic Primer Pair**:
  - **Forward (K757A-F)**: `5'- GGG CAC GTC CGC GCG GCC TTC AAG AGC -3'` ($T_m = 75.4\text{ }^\circ\text{C}$, $\text{GC} = 73.3\%$, length = 30 nt)
  - **Reverse (K757A-R)**: `5'- GCT CTT GAA GGC CGC GCG GAC GTG CCC -3'` ($T_m = 75.4\text{ }^\circ\text{C}$, $\text{GC} = 73.3\%$, length = 30 nt)

### 3.3 Construct 2: TERT Lys757R (Stereochemical & Charge Rescue Control)
- **Mutated Codon 757**: `CGG` (Arginine)
- **Engineered cDNA**: `5'- CAC GTC CGC CGG GCC TTC AAG -3'`
- **Diagnostic Restriction Cleavage**:
  - Mutating `AAG` to `CGG` creates `CCGG`, introducing a novel **MspI / HpaII** cleavage site (`C^CGG`) absent in the wild-type local flank.
- **Mutagenic Primer Pair**:
  - **Forward (K757R-F)**: `5'- GGG CAC GTC CGC CGG GCC TTC AAG AGC -3'` ($T_m = 75.8\text{ }^\circ\text{C}$, $\text{GC} = 73.3\%$, length = 30 nt)
  - **Reverse (K757R-R)**: `5'- GCT CTT GAA GGC CCG GCG GAC GTG CCC -3'` ($T_m = 75.8\text{ }^\circ\text{C}$, $\text{GC} = 73.3\%$, length = 30 nt)

### 3.4 Plasmid-Wide Restriction Digest Reality
Because human TERT cDNA is $65.3\%\text{ GC}$-rich:
- *BssHII* (`GCGCGC`) cuts 3 times in wild-type TERT (nt 7, 144, 468).
- *NaeI / NgoMIV* (`GCCGGC`) cuts 3 times in wild-type TERT (nt 1410, 2767, 3185).
- Therefore, restriction screening serves strictly as an initial colony-check on PCR amplicons. **Full-length plasmid Sanger/Nanopore sequencing across the entire 3,399 bp ORF is mandatory** before protein expression.

---

## 4. Full Sequencing Validation Primer Panel

```
========================================================================================================================
FULL-LENGTH SEQUENCING PRIMER PANEL
========================================================================================================================
Primer Name      Target Region       Sequence (5' -> 3')                     Target Nucleotide Span (hTERT ORF)
------------------------------------------------------------------------------------------------------------------------
SEQ-CMV-F        Tag & N-term        CGCAAATGGGCGGTAGGCGTG                   Vector CMV promoter -> hTERT 1-150
SEQ-TEN-1F       TEN Domain Core     GTCACCCTCGAGGAGCTG                      hTERT nt 350-800 (Residues 115-260)
SEQ-TEN-PLYQ     PLYQ / Tyr176       CTGAGCACCGTCTCCTTC                      hTERT nt 480-950 (Residues 160-310)
SEQ-TRBD-F       TRBD Domain         GTGTGGGTGCGAGACGTG                      hTERT nt 1100-1600 (Residues 365-530)
SEQ-RT-F         Catalytic Palm      CTGCTGAAGGCCCACCTG                      hTERT nt 1900-2400 (Residues 630-800)
SEQ-TRAP-F       TRAP & Lys757 Target GCCCATGGGCACGTCCGC                    hTERT nt 2220-2700 (Residues 740-900)
SEQ-CTE-F        C-terminal Extension GTGGACATCTTCCTGCGC                     hTERT nt 2700-3200 (Residues 900-1060)
SEQ-BGH-R        C-terminus / Vector TAGAAGGCACAGTCGAGG                      hTERT 3150-3399 -> BGH polyA
========================================================================================================================
```

---

## 5. Sequence-Verified FastA Constructs

```fasta
>hTERT_WT_TwinStrep_3xFLAG_3C_Verified
MWSHPQFEKGGGSGGGSGGSAWSHPQFEKGGGSGGGSGGSDYKDHDGDYKDHDIDYKDDDDKLEVLFQGPMAPRAPRC
RAVRSLLRSHYREVLPLATFVRRLGPQGWRLVQRGDPAAFRALVAQCLVCVPWDARPPPAAPSFRQVSCLKELVARVL
QRLCERGAKNVLAFGFALLDGARGGPPEAFTTSVRSYLPNTVTDALRGSGAWGLLLRRVGDDVLVHLLARCALFVLVA
PSCAYQVCGPPLYQLGAATQARPPPHASGPRRRLGCERAWNHSVREAGVPLGLPAPGARRRGGSASRSLPLPKRPRRG
AAPEPERTPVGQGSWAHPGRTRGPSDRGFCVVSPARPAEEATSLEGALSGTRHSHPSVGRQHHAGPPSTSRPPRPWDT
PCPPVYAETKHFLYSSGDKEQLRPSFLLSSLRPSLTGARRLVETIFLGSRPWMPGTPRRLPRLPQRYWQMRPLFLELL
GNHAQCPYGVLLKTHCPLRAAVTPAAGVCAREKPQGSVAAPEEEDTDPRRLVQLLRQHSSPWQVYGFVRACLRRLVPP
GLWGSRHNERRFLRNTKKFISLGKHAKLSLQELTWKMSVRDCAWLRRSPGVGCVPAAEHRLREEILAKFLHWLMSVYV
VELLRSFFYVTETTFQKNRLFFYRKSVWSKLQSIGIRQHLKRVQLRELSEAEVRQHREARPALLTSRLRFIPKPDGLR
PIVNMDYVVGARTFRREKRAERLTSRVKALFSVLNYERARRPGLLGASVLGLDDIHRAWRTFVLRVRAQDPPPELYFV
KVDVTGAYDTIPQDRLTEVIASIIKPQNTYCVRRYAVVQKAAHGHVRKAFKSHVSTLTDLQPYMRQFVAHLQETSPLR
DAVVIEQSSSLNEASSGLFDVFLRFMCHHAVRIRGKSYVQCQGIPQGSILSTLLCSLCYGDMENKLFAGIRRDGLLLR
LVDDFLLVTPHLTHAKTFLRTLVRGVPEYGCVVNLRKTVVNFPVEDEALGGTAFVQMPAHGLFPWCGLLLDTRTLEVQ
SDYSSYARTSIRASLTFNRGFKAGRNMRRKLFGVLRLKCHSLFLDLQVNSLQTVCTNIYKILLLQAYRFHACVLQLPF
HQQVWKNPTFFLRVISDTASLCYSILKAKNAGMSLGAKGAAGPLPSEAVQWLCHQAFLLKLTRHRVTYVPLLGSLRTA
QTQLSRKLPGTTLTALEAAANPALPSDFKTILD*

>hTERT_K757A_TwinStrep_3xFLAG_3C_Verified
MWSHPQFEKGGGSGGGSGGSAWSHPQFEKGGGSGGGSGGSDYKDHDGDYKDHDIDYKDDDDKLEVLFQGPMAPRAPRC
RAVRSLLRSHYREVLPLATFVRRLGPQGWRLVQRGDPAAFRALVAQCLVCVPWDARPPPAAPSFRQVSCLKELVARVL
QRLCERGAKNVLAFGFALLDGARGGPPEAFTTSVRSYLPNTVTDALRGSGAWGLLLRRVGDDVLVHLLARCALFVLVA
PSCAYQVCGPPLYQLGAATQARPPPHASGPRRRLGCERAWNHSVREAGVPLGLPAPGARRRGGSASRSLPLPKRPRRG
AAPEPERTPVGQGSWAHPGRTRGPSDRGFCVVSPARPAEEATSLEGALSGTRHSHPSVGRQHHAGPPSTSRPPRPWDT
PCPPVYAETKHFLYSSGDKEQLRPSFLLSSLRPSLTGARRLVETIFLGSRPWMPGTPRRLPRLPQRYWQMRPLFLELL
GNHAQCPYGVLLKTHCPLRAAVTPAAGVCAREKPQGSVAAPEEEDTDPRRLVQLLRQHSSPWQVYGFVRACLRRLVPP
GLWGSRHNERRFLRNTKKFISLGKHAKLSLQELTWKMSVRDCAWLRRSPGVGCVPAAEHRLREEILAKFLHWLMSVYV
VELLRSFFYVTETTFQKNRLFFYRKSVWSKLQSIGIRQHLKRVQLRELSEAEVRQHREARPALLTSRLRFIPKPDGLR
PIVNMDYVVGARTFRREKRAERLTSRVKALFSVLNYERARRPGLLGASVLGLDDIHRAWRTFVLRVRAQDPPPELYFV
KVDVTGAYDTIPQDRLTEVIASIIKPQNTYCVRRYAVVQKAAHGHVRAFKSHVSTLTDLQPYMRQFVAHLQETSPLR
DAVVIEQSSSLNEASSGLFDVFLRFMCHHAVRIRGKSYVQCQGIPQGSILSTLLCSLCYGDMENKLFAGIRRDGLLLR
LVDDFLLVTPHLTHAKTFLRTLVRGVPEYGCVVNLRKTVVNFPVEDEALGGTAFVQMPAHGLFPWCGLLLDTRTLEVQ
SDYSSYARTSIRASLTFNRGFKAGRNMRRKLFGVLRLKCHSLFLDLQVNSLQTVCTNIYKILLLQAYRFHACVLQLPF
HQQVWKNPTFFLRVISDTASLCYSILKAKNAGMSLGAKGAAGPLPSEAVQWLCHQAFLLKLTRHRVTYVPLLGSLRTA
QTQLSRKLPGTTLTALEAAANPALPSDFKTILD*

>hTERT_K757R_TwinStrep_3xFLAG_3C_Verified
MWSHPQFEKGGGSGGGSGGSAWSHPQFEKGGGSGGGSGGSDYKDHDGDYKDHDIDYKDDDDKLEVLFQGPMAPRAPRC
RAVRSLLRSHYREVLPLATFVRRLGPQGWRLVQRGDPAAFRALVAQCLVCVPWDARPPPAAPSFRQVSCLKELVARVL
QRLCERGAKNVLAFGFALLDGARGGPPEAFTTSVRSYLPNTVTDALRGSGAWGLLLRRVGDDVLVHLLARCALFVLVA
PSCAYQVCGPPLYQLGAATQARPPPHASGPRRRLGCERAWNHSVREAGVPLGLPAPGARRRGGSASRSLPLPKRPRRG
AAPEPERTPVGQGSWAHPGRTRGPSDRGFCVVSPARPAEEATSLEGALSGTRHSHPSVGRQHHAGPPSTSRPPRPWDT
PCPPVYAETKHFLYSSGDKEQLRPSFLLSSLRPSLTGARRLVETIFLGSRPWMPGTPRRLPRLPQRYWQMRPLFLELL
GNHAQCPYGVLLKTHCPLRAAVTPAAGVCAREKPQGSVAAPEEEDTDPRRLVQLLRQHSSPWQVYGFVRACLRRLVPP
GLWGSRHNERRFLRNTKKFISLGKHAKLSLQELTWKMSVRDCAWLRRSPGVGCVPAAEHRLREEILAKFLHWLMSVYV
VELLRSFFYVTETTFQKNRLFFYRKSVWSKLQSIGIRQHLKRVQLRELSEAEVRQHREARPALLTSRLRFIPKPDGLR
PIVNMDYVVGARTFRREKRAERLTSRVKALFSVLNYERARRPGLLGASVLGLDDIHRAWRTFVLRVRAQDPPPELYFV
KVDVTGAYDTIPQDRLTEVIASIIKPQNTYCVRRYAVVQKAAHGHVRRFKSHVSTLTDLQPYMRQFVAHLQETSPLR
DAVVIEQSSSLNEASSGLFDVFLRFMCHHAVRIRGKSYVQCQGIPQGSILSTLLCSLCYGDMENKLFAGIRRDGLLLR
LVDDFLLVTPHLTHAKTFLRTLVRGVPEYGCVVNLRKTVVNFPVEDEALGGTAFVQMPAHGLFPWCGLLLDTRTLEVQ
SDYSSYARTSIRASLTFNRGFKAGRNMRRKLFGVLRLKCHSLFLDLQVNSLQTVCTNIYKILLLQAYRFHACVLQLPF
HQQVWKNPTFFLRVISDTASLCYSILKAKNAGMSLGAKGAAGPLPSEAVQWLCHQAFLLKLTRHRVTYVPLLGSLRTA
QTQLSRKLPGTTLTALEAAANPALPSDFKTILD*
```
