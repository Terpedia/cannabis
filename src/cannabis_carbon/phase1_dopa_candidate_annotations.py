"""Sequence-verified annotation context for DOPA-lyase experimental leads."""
import hashlib
import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def run():
    source = Path('data/reports/phase1-dopa-lyase-search.json')
    search = json.loads(source.read_text())
    sequences = {p['accession']: p['sequence'] for p in search['cannabis_candidates']}
    raw = Path('data/raw/dopa-candidate-annotations')
    raw.mkdir(parents=True, exist_ok=True)

    def fetch(acc):
        path = raw / (acc + '.json')
        url = 'https://rest.uniprot.org/uniprotkb/' + acc + '.json'
        if not path.exists():
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = response.read()
            data = json.loads(payload)
            if data['primaryAccession'] != acc or data['sequence']['value'] != sequences[acc]:
                raise ValueError('Annotation and searched sequence differ')
            path.write_bytes(payload)
        data = json.loads(path.read_text())
        if data['primaryAccession'] != acc or data['sequence']['value'] != sequences[acc]:
            raise ValueError('Cached annotation and searched sequence differ')
        return {'accession': acc, 'source_url': url, 'snapshot': str(path),
            'snapshot_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'annotation': data, 'model_eligible': False,
            'passing_alignment_ids': [a['id'] for a in search['passing_alignments'] if a['cannabis_accession'] == acc],
            'proposed_test': 'Compare L-DOPA, L-tyrosine and L-phenylalanine turnover and identify caffeate, coumarate and cinnamate independently. Include substrate-autoxidation and no-enzyme controls. Reference L-DOPA activity requires independent validation.',
            'claim_boundary': 'Protein/gene/domain annotations are retained with original evidence codes. Sequence identity to the searched proteome is verified; annotation is not a Cannabis biochemical assay or exact substrate assignment.'}

    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(fetch, sorted(sequences)))
    report = {'schema': 'cannabis-dopa-candidate-annotation-context-v1', 'rows': rows,
        'model_eligible': False, 'source_sha256': {str(source): hashlib.sha256(source.read_bytes()).hexdigest(),
            **{r['snapshot']: r['snapshot_sha256'] for r in rows}},
        'summary': {'sequence_verified_protein_leads': len(rows), 'new_exact_enzyme_assignments': 0},
        'claim_boundary': 'Experimental-lead context only. No model integration, reaction-direction resolution, net-route rescue or atom-tracing claim.'}
    Path('data/reports/phase1-dopa-candidate-annotations.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
