"""Test a competing substrate-class hypothesis without changing exact reaction assignments."""
import hashlib
import json
from pathlib import Path
from .phase1_family_search import parse_references
from .phase1_new_protein_search import screen


def run():
    raw = Path('data/raw/sphingolipid-alternative-search')
    references, matches, retrievals, hashes = {}, [], [], {}
    for accession in ('Q9ZRP7', 'Q3EBF7'):
        path = raw / (accession + '.json')
        protein = json.loads(path.read_text())
        if protein['primaryAccession'] != accession:
            raise ValueError('Reference identity mismatch')
        references.update(parse_references(f'>{accession}\n{protein["sequence"]["value"]}\n'.encode(), {accession}))
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        hashes[str(path)] = digest
        url = f'https://rest.uniprot.org/uniprotkb/{accession}.json'
        matches.append({'accession': accession, 'model_eligible': False,
            'match_type': 'alternative-substrate-class-reference-not-exact-delta6',
            'source_activity': [c for c in protein['comments'] if c['commentType'] == 'CATALYTIC ACTIVITY']})
        retrievals.append({'requested_accessions': [accession], 'status': 'retrieved',
                          'missing_accessions': [], 'snapshot': str(path), 'sha256': digest, 'url': url})
    boundary = ('Competing sphingolipid delta-8 substrate-class hypothesis, screened against the full Cannabis proteome. '
                'This is not an exact delta-6 acyl-CoA reaction join. Reference assay evidence remains source annotation; '
                'Cannabis specificity and product geometry require experiments. No pathway-model promotion.')
    discovery = {'schema': 'cannabis-sphingolipid-alternative-discovery-v1', 'source_sha256': hashes,
        'rows': [{'reaction_id': 'alternative-hypothesis:sphingolipid-delta8',
                  'target_ids': [], 'priority_target_ids': [], 'hypothesis_ids': [],
                  'reference_matches': matches}], 'claim_boundary': boundary}
    source = Path('data/reports/phase1-sphingolipid-alternative-references.json')
    source.write_text(json.dumps(discovery, separators=(',', ':')) + '\n')
    screen(discovery, source, raw, Path('data/reports/phase1-sphingolipid-alternative-search.json'),
           references, retrievals, evidence_class='alternative-substrate-class-homology-lead',
           additional_blockers=('not-an-exact-delta6-reaction-assignment',), claim_boundary=boundary)


if __name__ == '__main__':
    run()
