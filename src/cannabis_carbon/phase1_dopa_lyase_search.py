"""Full-proteome search of a literature-reported DOPA activity lead."""
import hashlib
import json
import urllib.request
from pathlib import Path
from .phase1_family_search import parse_references
from .phase1_new_protein_search import screen, export_table


def run():
    review_path = Path('data/curation/dopa-lyase-reference-review.json')
    parent_path = Path('data/reports/phase1-current-gap-priority.json')
    review, parent = [json.loads(p.read_text()) for p in (review_path, parent_path)]
    gap = next(r for r in parent['rows'] if r['reaction_id'] == review['reaction_id'])
    selected_uses = []
    for use in gap['selected_uses']:
        tids = [t for t in use['target_ids'] if t in gap['remaining_target_ids']]
        if tids:
            selected_uses.append({**use, 'target_ids': tids})
    accession = review['reference_accession']
    raw = Path('data/raw/phase1-dopa-lyase-search')
    raw.mkdir(parents=True, exist_ok=True)
    snapshot = raw / (accession + '.json')
    url = 'https://rest.uniprot.org/uniprotkb/' + accession + '.json'
    if not snapshot.exists():
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = response.read()
        if json.loads(payload)['primaryAccession'] != accession:
            raise ValueError('Reference accession mismatch')
        snapshot.write_bytes(payload)
    protein = json.loads(snapshot.read_text())
    if protein['primaryAccession'] != accession:
        raise ValueError('Reference accession mismatch')
    references = parse_references(('>' + accession + '\n' + protein['sequence']['value'] + '\n').encode(), {accession})
    boundary = ('Literature-reported L-DOPA activity with assay data not shown; product-binding '
        'structure is not turnover evidence. UniProt annotates a different substrate. Homology leads '
        'remain ineligible for exact-reaction integration pending substrate and reference-activity review. '
        'No Cannabis activity, physiological direction, complete pathway, or atom-tracing claim.')
    discovery = {'schema': 'cannabis-dopa-lyase-lead-discovery-v1', 'model_eligible': False,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (review_path, parent_path, snapshot)},
        'claim_boundary': boundary, 'rows': [{
            'reaction_id': gap['reaction_id'], 'left': gap['reaction']['left'], 'right': gap['reaction']['right'],
            'sources': gap['reaction']['sources'], 'target_ids': gap['remaining_target_ids'],
            'priority_target_ids': gap['remaining_target_ids'], 'hypothesis_ids': [],
            'selected_uses': selected_uses, 'historical_selected_uses': gap['selected_uses'],
            'prior_searches': gap['prior_searches'],
            'model_eligible': False, 'reference_matches': [{
                'accession': accession, 'model_eligible': False, 'exact_reaction_annotation_match': False,
                'match_type': 'literature-reported-activity; assay-data-not-shown',
                'review': review, 'annotation_comments': protein.get('comments', []),
                'source_url': url}]}]}
    source = Path('data/reports/phase1-dopa-lyase-references.json')
    source.write_text(json.dumps(discovery, separators=(',', ':')) + '\n')
    retrievals = [{'requested_accessions': [accession], 'status': 'retrieved', 'url': url,
        'snapshot': str(snapshot), 'sha256': hashlib.sha256(snapshot.read_bytes()).hexdigest(),
        'missing_accessions': [], 'sequence_source': 'UniProt JSON sequence.value'}]
    output = Path('data/reports/phase1-dopa-lyase-search.json')
    screen(discovery, source, raw, output, references, retrievals,
        evidence_class='reported-reference-activity-specificity-unverified-homology-lead',
        additional_blockers=('reference-L-DOPA-assay-data-not-shown', 'not-eligible-for-exact-reaction-model'),
        claim_boundary=boundary)
    result = json.loads(output.read_text())
    result['model_eligible'] = False
    for row in result['rows']:
        row['model_eligible'] = False
    output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    export_table(output, Path('data/derived/phase1-dopa-lyase-search.ndjson'))


if __name__ == '__main__':
    run()
