"""Pin the source prediction behind the leading CHI sequence; no transcript claim."""
import hashlib
import json
import re
import urllib.request
from pathlib import Path


def parse_record(payload):
    text = payload.decode()
    version = re.search(r'^VERSION\s+(\S+)', text, re.M).group(1)
    sequence = re.sub('[^a-zA-Z]', '', text.split('ORIGIN', 1)[1].split('//', 1)[0]).upper()
    coded_by = re.search(r'/coded_by="([^"]+)"', text).group(1)
    coded_by = re.sub(r'\s+', '', coded_by)
    if not coded_by.startswith('complement(join(') or not coded_by.endswith('))'):
        raise ValueError('Unexpected coding-location grammar')
    pieces = coded_by[len('complement(join('):-2].split(',')
    exons = []
    for piece in pieces:
        match = re.fullmatch(r'([A-Z0-9]+\.\d+):(\d+)\.\.(\d+)', piece)
        if not match:
            raise ValueError('Unresolved source exon coordinates')
        accession, start, end = match.groups()
        start, end = int(start), int(end)
        if start > end:
            raise ValueError('Invalid genomic interval')
        exons.append({'accession': accession, 'start': start, 'end': end, 'strand': -1, 'length_nt': end-start+1})
    return {'protein_version': version, 'sequence': sequence, 'coded_by': coded_by,
        'genomic_order_exons': exons, 'transcript_order_exons': list(reversed(exons)),
        'coding_length_nt': sum(e['length_nt'] for e in exons)}


def run():
    annotation_path = Path('data/raw/chalcone-annotations/A0A7J6I409.json')
    site_path = Path('data/reports/phase1-chalcone-site-review.json')
    annotation = json.loads(annotation_path.read_text())
    source = Path('data/raw/chalcone-gene-model/KAF4401769.1.gb')
    source.parent.mkdir(parents=True, exist_ok=True)
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=KAF4401769.1&rettype=gb&retmode=text'
    if not source.exists():
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = response.read()
        parsed = parse_record(payload)
        if parsed['protein_version'] != 'KAF4401769.1' or parsed['sequence'] != annotation['sequence']['value']:
            raise ValueError('Source prediction differs from screened protein')
        source.write_bytes(payload)
    parsed = parse_record(source.read_bytes())
    if parsed['protein_version'] != 'KAF4401769.1' or parsed['sequence'] != annotation['sequence']['value']:
        raise ValueError('Cached prediction mismatch')
    sites = next(r for r in json.loads(site_path.read_text())['rows'] if r['accession'] == 'A0A7J6I409')
    domain = [f for f in annotation['features'] if f['type'] == 'Domain']
    report = {'schema': 'cannabis-chalcone-gene-model-review-v1', 'model_eligible': False,
        'accession': 'A0A7J6I409', 'source_url': url, 'source_record': parsed,
        'gene_annotation': annotation['genes'], 'domain_annotations': domain,
        'reference_alignment': sites,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, annotation_path, site_path)},
        'summary': {'protein_length_aa': len(parsed['sequence']), 'exons': len(parsed['genomic_order_exons']),
            'coding_length_nt': parsed['coding_length_nt'], 'new_exact_enzyme_assignments': 0},
        'next_test': 'Retrieve the exact genomic intervals and verify spliced translation, then seek independent transcript or proteomic support for the terminal extensions and splice junctions. Do not trim the sequence or treat annotation-domain endpoints as construct boundaries. Test stereoselective activity with chiral product analysis and spontaneous-conversion controls.',
        'claim_boundary': 'The pinned NCBI conceptual translation matches the screened UniProt sequence, but these records derive from the same gene prediction and are not independent expression evidence. Coding-location arithmetic is not a genomic translation check. The extra sequence is unresolved: no claim of artifact, transit peptide, functional fusion, mature protein boundary, or enzyme activity.'}
    Path('data/reports/phase1-chalcone-gene-model.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
