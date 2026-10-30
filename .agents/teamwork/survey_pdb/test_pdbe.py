import requests

for p in ['7TRE', '7TRD', '7TRF', '7TRC', '2I46']:
    pid = p.lower()
    pdf_url = f'https://www.ebi.ac.uk/pdbe/entry-files/download/{pid}_full_validation.pdf'
    xml_url = f'https://www.ebi.ac.uk/pdbe/entry-files/download/{pid}_validation.xml.gz'
    r_pdf = requests.head(pdf_url, headers={'User-Agent': 'Mozilla/5.0'})
    r_xml = requests.head(xml_url, headers={'User-Agent': 'Mozilla/5.0'})
    print(f"{p}: PDF={r_pdf.status_code} ({r_pdf.headers.get('content-length', '?')} bytes), XML={r_xml.status_code} ({r_xml.headers.get('content-length', '?')} bytes)")
