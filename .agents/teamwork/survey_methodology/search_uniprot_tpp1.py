import urllib.request
import json

url = "https://rest.uniprot.org/uniprotkb/search?query=gene:ACD+AND+organism_id:9606&format=json"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req) as resp:
    res = json.load(resp)
    for r in res.get('results', []):
        acc = r.get('primaryAccession')
        entry_name = r.get('uniProtkbId')
        seq = r.get('sequence', {}).get('value')
        print(f"UniProt {acc} ({entry_name}) len={len(seq)}")
        # Check isoforms
        comments = r.get('comments', [])
        for c in comments:
            if c.get('commentType') == 'ALTERNATIVE PRODUCTS':
                for iso in c.get('isoforms', []):
                    print("  Isoform:", iso.get('name', {}).get('value'), iso.get('isoformIds'))
