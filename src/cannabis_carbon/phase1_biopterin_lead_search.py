"""Screen a generic-substrate lead without promoting an exact reaction assignment."""
import hashlib
import json
from pathlib import Path

from .phase1_family_search import parse_references
from .phase1_new_protein_search import screen


def run():
    review_path = Path('data/curation/biopterin-reference-review.json')
    review = json.loads(review_path.read_text())
    snapshot = Path(review['source']['raw_path'])
    raw_bytes = snapshot.read_bytes()
    if hashlib.sha256(raw_bytes).hexdigest() != review['source']['sha256']:
        raise ValueError('Reference annotation snapshot changed')
    protein = json.loads(raw_bytes)
    accession = protein['primaryAccession']
    if accession != review['source']['accession']:
        raise ValueError('Reference accession mismatch')
    sequence = protein['sequence']['value']
    references = parse_references(f'>{accession}\n{sequence}\n'.encode(), {accession})
    parent_path = Path('data/reports/phase1-nonplant-reference-review.json')
    parent = json.loads(parent_path.read_text())
    gap = next(r for r in parent['rows'] if r['reaction_id'] == review['reaction_id'])
    boundary = ('Generic-substrate pteridine reductase reference screening across the full pinned Cannabis proteome. '
                'Matches are specificity-unverified experimental leads, NOT exact-Rhea enzyme candidates eligible '
                'for the balanced pathway model. Cofactor, stereochemistry, direction, compartment and all-input '
                'supply require independent evidence. Atom tracing remains deferred.')
    discovery = {
        'schema': 'cannabis-biopterin-generic-lead-discovery-v1',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (review_path, snapshot, parent_path)},
        'claim_boundary': boundary,
        'rows': [{
            'reaction_id': gap['reaction_id'], 'target_ids': gap['target_ids'],
            'priority_target_ids': gap.get('priority_target_ids', []),
            'hypothesis_ids': gap.get('hypothesis_ids', []),
            'reference_matches': [{
                'accession': accession,
                'match_type': 'generic-substrate-literature-lead-not-exact-rhea',
                'exact_reaction_annotation_match': False,
                'model_eligible': False,
                'annotation_observations': review['annotation_observations'],
                'source': review['source']
            }]
        }]
    }
    source = Path('data/reports/phase1-biopterin-lead-references.json')
    source.write_text(json.dumps(discovery, separators=(',', ':')) + '\n')
    retrievals = [{
        'requested_accessions': [accession], 'status': 'retrieved',
        'url': review['source']['url'], 'snapshot': str(snapshot),
        'sha256': review['source']['sha256'], 'missing_accessions': [],
        'sequence_source': 'UniProt JSON sequence.value; no independent FASTA retrieval'
    }]
    screen(discovery, source, Path('data/raw/phase1-biopterin-lead-search'),
           Path('data/reports/phase1-biopterin-lead-search.json'), references, retrievals,
           evidence_class='generic-substrate-specificity-unverified-homology-lead',
           additional_blockers=('not-eligible-for-exact-reaction-model',
                                'generic-reference-substrate-not-exact-biopterin-assignment'),
           claim_boundary=boundary)


if __name__ == '__main__':
    run()
