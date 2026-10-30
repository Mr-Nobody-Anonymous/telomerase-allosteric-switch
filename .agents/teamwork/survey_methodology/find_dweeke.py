import urllib.request

url = "https://rest.uniprot.org/uniprotkb/Q96AP0.fasta"
with urllib.request.urlopen(url) as resp:
    fasta = resp.read().decode('utf-8')
seq = "".join(fasta.splitlines()[1:])

print("UniProt length:", len(seq))
print("First 100:", seq[:100])

# Search for DWEEKE
pos_dwee = seq.find("DWEEKE")
print(f"'DWEEKE' found at 0-based {pos_dwee} (1-based {pos_dwee+1})")
if pos_dwee != -1:
    print(f"Residues at {pos_dwee-2} to {pos_dwee+10}: {seq[pos_dwee-2:pos_dwee+10]}")
    # In 'DWEEKE':
    # D is at pos_dwee+1
    # W is at pos_dwee+2
    # E is at pos_dwee+3
    # E is at pos_dwee+4
    # K is at pos_dwee+5
    # E is at pos_dwee+6
    print(f"Residue {pos_dwee+1} = {seq[pos_dwee]}")
    print(f"Residue {pos_dwee+2} = {seq[pos_dwee+1]}")
    print(f"Residue {pos_dwee+3} = {seq[pos_dwee+2]}")
    print(f"Residue {pos_dwee+4} = {seq[pos_dwee+3]}")
    print(f"Residue {pos_dwee+5} = {seq[pos_dwee+4]}")
    print(f"Residue {pos_dwee+6} = {seq[pos_dwee+5]}")
