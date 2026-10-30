import requests

pdb_list = ['7TRE', '7TRD', '7TRF', '7TRC', '2I46']
print("Testing official download endpoints:")
for p in pdb_list:
    mid = p[1:3].lower()
    pid = p.lower()
    pdf_url = f'https://files.rcsb.org/pub/pdb/validation_reports/{mid}/{pid}/{pid}_full_validation.pdf'
    xml_url = f'https://files.rcsb.org/pub/pdb/validation_reports/{mid}/{pid}/{pid}_validation.xml.gz'
    cif_url = f'https://files.rcsb.org/download/{p}.cif'
    cif_gz_url = f'https://files.rcsb.org/download/{p}.cif.gz'
    
    r_pdf = requests.head(pdf_url, headers={'User-Agent': 'Mozilla/5.0'})
    r_xml = requests.head(xml_url, headers={'User-Agent': 'Mozilla/5.0'})
    r_cif = requests.head(cif_url, headers={'User-Agent': 'Mozilla/5.0'})
    r_cif_gz = requests.head(cif_gz_url, headers={'User-Agent': 'Mozilla/5.0'})
    
    print(f"[{p}]")
    print(f"  mmCIF:     HTTP {r_cif.status_code} ({r_cif.headers.get('Content-Length', '?')} bytes) -> {cif_url}")
    print(f"  mmCIF.gz:  HTTP {r_cif_gz.status_code} ({r_cif_gz.headers.get('Content-Length', '?')} bytes) -> {cif_gz_url}")
    print(f"  Val PDF:   HTTP {r_pdf.status_code} ({r_pdf.headers.get('Content-Length', '?')} bytes) -> {pdf_url}")
    print(f"  Val XML.gz:HTTP {r_xml.status_code} ({r_xml.headers.get('Content-Length', '?')} bytes) -> {xml_url}")

# Also check EMDB map URLs
emdb_maps = {
    '7TRE': 'EMD-26087',
    '7TRD': 'EMD-26086',
    '7TRF': 'EMD-26088',
    '7TRC': 'EMD-26085'
}
print("\nTesting EMDB Map Endpoints:")
for p, emd in emdb_maps.items():
    emd_num = emd.split('-')[1]
    emdb_url = f'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{emd}/map/emd_{emd_num}.map.gz'
    r_emd = requests.head(emdb_url, headers={'User-Agent': 'Mozilla/5.0'})
    print(f"[{p} / {emd}] Map: HTTP {r_emd.status_code} ({r_emd.headers.get('Content-Length', '?')} bytes) -> {emdb_url}")
