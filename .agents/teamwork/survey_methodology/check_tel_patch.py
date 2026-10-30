import json

with open(r'C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\survey_raw.json', 'r') as f:
    data = json.load(f)

# Inspect UniProt sequence for human TPP1 (Q96AP0) if available or query it
import urllib.request

url = "https://rest.uniprot.org/uniprotkb/Q96AP0.fasta"
try:
    with urllib.request.urlopen(url) as resp:
        fasta = resp.read().decode('utf-8')
    uniprot_seq = "".join(fasta.splitlines()[1:])
    print(f"UniProt Q96AP0 length: {len(uniprot_seq)}")
    # Find TEL patch residues around 168-171, 183, 212
    # In UniProt Q96AP0:
    print(f"UniProt 165-175: {uniprot_seq[164:175]} (indices 165 to 175, 1-based)")
    print(f"UniProt 168: {uniprot_seq[167]}, 169: {uniprot_seq[168]}, 170: {uniprot_seq[169]}, 171: {uniprot_seq[170]}")
    print(f"UniProt 180: {uniprot_seq[179]}, 183: {uniprot_seq[182]}, 212: {uniprot_seq[211]}")
    
    # Now find where 7TRE TPP1 aligns with UniProt Q96AP0
    tpp1_7tre = data['7TRE']['entities']['4']['entity_poly']['pdbx_seq_one_letter_code_can']
    print("\nAligning 7TRE TPP1 to UniProt:")
    # search substring
    sub = tpp1_7tre[4:24] # skip GAGS
    pos = uniprot_seq.find(sub)
    print(f"Substring '{sub}' found at UniProt 0-based {pos} (1-based {pos+1})")
    print(f"So 7TRE residue at index 4 (1-based 5, 'G') corresponds to UniProt {pos+1}")
except Exception as e:
    print("Error querying UniProt:", e)
