import urllib.request

from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq

# Fetch canonical NM_198253.3 CDS
url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_198253.3&rettype=fasta_cds_na&retmode=text'
with urllib.request.urlopen(url, timeout=10) as response:
    text = response.read().decode('utf-8')
    cds = ''.join([line.strip() for line in text.split('\n') if line and not line.startswith('>')])

wt_seq = Seq(cds)
# Mutate 757 AAG to GCG (K757A) and CGG (K757R)
start = (757 - 1) * 3
k757a_seq = Seq(cds[:start] + 'GCG' + cds[start+3:])
k757r_seq = Seq(cds[:start] + 'CGG' + cds[start+3:])

print("1. Restriction sites in WT hTERT CDS:")
rb = RestrictionBatch(['BssHII', 'NgoMIV', 'NaeI', 'NotI', 'AscI', 'FseI', 'PacI', 'PmeI', 'SbfI', 'AgeI', 'XhoI', 'EcoRI', 'BamHI', 'HindIII'])
analysis_wt = rb.search(wt_seq)
for enz, sites in analysis_wt.items():
    print(f"  {enz}: {len(sites)} sites at {sites}")

print("\n2. Comparing WT vs K757A (AAG -> GCG):")
analysis_k757a = rb.search(k757a_seq)
for enz in rb:
    diff = len(analysis_k757a[enz]) - len(analysis_wt[enz])
    if diff != 0:
        print(f"  {enz}: {len(analysis_wt[enz])} -> {len(analysis_k757a[enz])} sites (delta = {diff:+d}) at new site {analysis_k757a[enz][-1]}")

print("\n3. Comparing WT vs K757R (AAG -> CGG):")
analysis_k757r = rb.search(k757r_seq)
for enz in rb:
    diff = len(analysis_k757r[enz]) - len(analysis_wt[enz])
    if diff != 0:
        print(f"  {enz}: {len(analysis_wt[enz])} -> {len(analysis_k757r[enz])} sites (delta = {diff:+d}) at new site {analysis_k757r[enz][-1]}")
