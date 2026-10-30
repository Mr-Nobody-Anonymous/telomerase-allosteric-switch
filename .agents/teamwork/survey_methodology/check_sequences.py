import json

with open(r'C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\survey_raw.json', 'r') as f:
    data = json.load(f)

# 7TRE TPP1 sequence
tpp1_7tre_seq = data['7TRE']['entities']['4']['entity_poly']['pdbx_seq_one_letter_code_can']
print("7TRE TPP1 canonical sequence (len", len(tpp1_7tre_seq), "):")
print(tpp1_7tre_seq)

# 2I46 TPP1 sequence
tpp1_2i46_seq = data['2I46']['entities']['1']['entity_poly']['pdbx_seq_one_letter_code_can']
print("\n2I46 TPP1 canonical sequence (len", len(tpp1_2i46_seq), "):")
print(tpp1_2i46_seq)

# Check TERT canonical sequences across 7TRE, 7TRD, 7TRF
tert_7tre_seq = data['7TRE']['entities']['2']['entity_poly']['pdbx_seq_one_letter_code_can']
tert_7trd_seq = data['7TRD']['entities']['2']['entity_poly']['pdbx_seq_one_letter_code_can']
tert_7trf_seq = data['7TRF']['entities']['4']['entity_poly']['pdbx_seq_one_letter_code_can']

print(f"\nTERT sequence lengths: 7TRE={len(tert_7tre_seq)}, 7TRD={len(tert_7trd_seq)}, 7TRF={len(tert_7trf_seq)}")
print("7TRE == 7TRD?", tert_7tre_seq == tert_7trd_seq)
print("7TRE == 7TRF?", tert_7tre_seq == tert_7trf_seq)
