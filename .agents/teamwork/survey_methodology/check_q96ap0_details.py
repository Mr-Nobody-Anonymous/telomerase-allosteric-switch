import urllib.request
import json

# Fetch UniProt Q96AP0 json entry
url = "https://rest.uniprot.org/uniprotkb/Q96AP0.json"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req) as resp:
    data = json.load(resp)

print("Canonical sequence length:", len(data['sequence']['value']))

# Let's inspect comments for isoforms
for c in data.get('comments', []):
    if c.get('commentType') == 'ALTERNATIVE PRODUCTS':
        for iso in c.get('isoforms', []):
            print("Isoform name:", iso.get('name', {}).get('value'))
            print("Isoform ids:", iso.get('isoformIds'))
            print("Isoform sequence status:", iso.get('isoformSequenceStatus'))

# Also check features for OB-fold / domain
for f in data.get('features', []):
    if 'domain' in f.get('type', '').lower() or 'region' in f.get('type', '').lower():
        print("Feature:", f.get('type'), f.get('description'), f.get('location'))
