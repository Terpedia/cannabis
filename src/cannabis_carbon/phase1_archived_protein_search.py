"""Screen explicitly resolved original-source archived sequences, not replacements."""
import hashlib
import json
from pathlib import Path
from .phase1_archived_references import validate_entry
from .phase1_new_protein_search import screen, export_table


def run():
    source = Path('data/reports/phase1-archived-references.json')
    discovery = json.loads(source.read_text())
    for p, digest in discovery['source_sha256'].items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest() != digest:
            raise ValueError('Archive discovery source changed')
    references = {r['accession']: r for r in discovery['reference_sequences']}
    if len(references) != len(discovery['reference_sequences']):
        raise ValueError('Duplicate resolved archive sequence')
    for record in discovery['retrievals']:
        if record['status'] != 'retrieved':
            continue
        accession = record['requested_accession']
        if hashlib.sha256(record['response_text'].encode()).hexdigest() != record['response_sha256']:
            raise ValueError('Archive response changed')
        validated = validate_entry(accession, json.loads(record['response_text']))
        if accession in references and references[accession] != validated:
            raise ValueError('Resolved archive sequence differs from source response')
    output = Path('data/reports/phase1-archived-protein-search.json')
    screen(discovery, source, Path('data/raw/phase1-archived-protein-search'), output,
        references, discovery['retrievals'],
        evidence_class='original-MARTS-archived-sequence-homology-for-inferred-stoichiometry',
        additional_blockers=('archived-sequence-is-not-functional-annotation',
            'inferred-inorganic-stoichiometry-unverified', 'original-MARTS-exact-product-identity-unverified'),
        claim_boundary=discovery['claim_boundary'])
    count = export_table(output, Path('data/derived/phase1-archived-protein-search.ndjson'))
    print(json.dumps({'export_rows': count, 'sha256': hashlib.sha256(output.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    run()
