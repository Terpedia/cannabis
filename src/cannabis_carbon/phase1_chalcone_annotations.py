"""Sequence-verified full annotation context for every screened CHI reference."""
import hashlib
import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def run(name='chalcone'):
    if name not in ('chalcone', 'flavone'):
        raise ValueError('Unknown annotation batch')
    source = Path('data/reports/phase1-' + name + '-search.json')
    search = json.loads(source.read_text())
    sequences = {p['accession']: p['sequence'] for p in search['reference_sequences'] + search['cannabis_candidates']}
    refs = {p['accession'] for p in search['reference_sequences']}
    raw = Path('data/raw/' + name + '-annotations')
    raw.mkdir(parents=True, exist_ok=True)

    def fetch(acc):
        path = raw / (acc + '.json')
        url = 'https://rest.uniprot.org/uniprotkb/' + acc + '.json'
        if not path.exists():
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = response.read()
            data = json.loads(payload)
            if data['primaryAccession'] != acc or data['sequence']['value'] != sequences[acc]:
                raise ValueError('Annotation sequence does not match search')
            path.write_bytes(payload)
        data = json.loads(path.read_text())
        if data['primaryAccession'] != acc or data['sequence']['value'] != sequences[acc]:
            raise ValueError('Cached sequence mismatch')
        return {'accession': acc, 'role': 'reference' if acc in refs else 'Cannabis-lead',
            'annotation': data, 'source_url': url, 'snapshot': str(path),
            'snapshot_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'model_eligible': False,
            'passing_alignment_ids': [a['id'] for a in search['passing_alignments'] if acc in (a['reference_accession'], a['cannabis_accession'])]}

    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(fetch, sorted(sequences)))
    report = {'schema': 'cannabis-' + name + '-annotation-context-v1', 'model_eligible': False, 'rows': rows,
        'source_sha256': {str(source): hashlib.sha256(source.read_bytes()).hexdigest(), **{r['snapshot']: r['snapshot_sha256'] for r in rows}},
        'summary': {'sequence_verified_references': len(refs), 'sequence_verified_Cannabis_leads': len(sequences) - len(refs), 'new_exact_enzyme_assignments': 0},
        'claim_boundary': 'Every reference and passing Cannabis lead retained, including references with no passing hits. Original function, kinetic, caution, similarity and mutagenesis annotations and evidence codes are preserved without transferring activity. No individual reference is promoted merely because its EC annotation is reviewed. Full substrate, stereoselectivity, domain and sequence review remains required.'}
    Path('data/reports/phase1-' + name + '-annotations.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
