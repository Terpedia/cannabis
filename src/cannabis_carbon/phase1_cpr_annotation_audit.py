"""Name-annotated Cannabis reductase leads, not proven P450 partnerships."""
import hashlib
import json
import urllib.request
from pathlib import Path
from .genome import _fasta

NAME = 'NADPH--cytochrome P450 reductase'


def run():
    source = Path('data/reports/phase1-fnsii-search.json')
    parent = json.loads(source.read_text())
    proteome = Path(parent['proteome_path'])
    if hashlib.sha256(proteome.read_bytes()).hexdigest() != parent['proteome_sha256']:
        raise ValueError('Changed proteome')
    sequences = _fasta(proteome)
    headers = [line[1:] for line in proteome.read_text().splitlines() if line.startswith('>')]
    selected = [h for h in headers if NAME in h]
    raw = Path('data/raw/cpr-annotation-audit')
    raw.mkdir(parents=True, exist_ok=True)
    rows = []
    for header in selected:
        acc = header.split('|')[1]
        path = raw / (acc + '.json')
        url = 'https://rest.uniprot.org/uniprotkb/' + acc + '.json'
        if not path.exists():
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = response.read()
            annotation = json.loads(payload)
            if annotation['primaryAccession'] != acc or annotation['sequence']['value'] != sequences[acc]:
                raise ValueError('Annotation differs from pinned proteome')
            path.write_bytes(payload)
        annotation = json.loads(path.read_text())
        if annotation['primaryAccession'] != acc or annotation['sequence']['value'] != sequences[acc]:
            raise ValueError('Cached annotation differs from pinned proteome')
        rows.append({'accession': acc, 'source_header': header, 'sequence': sequences[acc],
            'annotation': annotation, 'snapshot': str(path), 'source_url': url,
            'model_eligible': False, 'evidence_class': 'name-annotation-only-reductase-lead',
            'compatible_fnsii_partners': [], 'partner_status': 'unverified'})
    if len({r['accession'] for r in rows}) != len(rows):
        raise ValueError('Duplicate candidate identity')
    paths = [source, proteome, *[Path(r['snapshot']) for r in rows]]
    report = {'schema': 'cannabis-cpr-annotation-audit-v1', 'model_eligible': False,
        'selection': {'header_contains': NAME, 'method': 'literal-name scan; not sequence-homology discovery'},
        'rows': rows, 'related_fnsii_candidate_accessions': sorted(p['accession'] for p in parent['cannabis_candidates']),
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'proteome_sequences': len(sequences), 'name_annotated_reductase_leads': len(rows),
            'verified_fnsii_partnerships': 0, 'new_exact_enzyme_assignments': 0},
        'next_steps': ['Use characterized plant CPR sequences for whole-proteome homology discovery, including proteins without CPR names.',
            'Review flavin/NADPH domains, membrane targeting, expression and paired P450 turnover.',
            'Reconcile carrier redox states and NADPH regeneration before a coupled balanced scenario.'],
        'claim_boundary': 'A name annotation and exact sequence match do not prove reductase activity or compatibility with any of the 95 FNS-II leads. Empty partner lists mean unverified, not incompatible. This scan cannot exclude unnamed or misannotated reductases. No organic seed, cofactor supply, reaction-model addition or atom-tracing claim.'}
    Path('data/reports/phase1-cpr-annotation-audit.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
