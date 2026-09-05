"""Whole-proteome search against a source-defined desaturase region."""
import hashlib
import json
from pathlib import Path
from .phase1_family_search import parse_references
from .phase1_new_protein_search import screen


def run():
    raw = Path('data/raw/desaturase-catalytic-domain')
    protein_path, domain_path = raw / 'O95864.json', raw / 'O95864-PF00487.json'
    protein, domain = [json.loads(p.read_text()) for p in (protein_path, domain_path)]
    match = domain['proteins'][0]
    if match['accession'].upper() != protein['primaryAccession'] or match['protein_length'] != protein['sequence']['length']:
        raise ValueError('Protein/domain identity mismatch')
    fragments = match['entry_protein_locations'][0]['fragments']
    if len(fragments) != 1 or fragments[0]['dc-status'] != 'CONTINUOUS':
        raise ValueError('Expected one continuous source-defined domain')
    start, end = fragments[0]['start'], fragments[0]['end']
    sequence = protein['sequence']['value'][start-1:end]
    identifier = f'O95864_PF00487_{start}_{end}'
    references = parse_references(f'>{identifier}\n{sequence}\n'.encode(), {identifier})
    activity = next(c for c in protein['comments'] if c['commentType'] == 'CATALYTIC ACTIVITY'
                    and any(x['id'] == 'RHEA:47144' for x in c['reaction']['reactionCrossReferences']))
    parent_path = Path('data/reports/phase1-remaining-gap-references.json')
    parent = json.loads(parent_path.read_text())
    gap = next(r for r in parent['rows'] if 'RHEA:47145' in r['source_reaction_ids'])
    boundary = ('Source-defined PF00487 catalytic-region homology search across the complete Cannabis proteome. '
                'The reference is a subsequence, not full-length O95864. Its exact Rhea activity is similarity-inferred. '
                'Passing hits are domain-level, specificity-unverified leads; not automatic exact-reaction model additions. '
                'Original full-protein screens are unchanged. Atom tracing deferred.')
    row = {k: gap[k] for k in ('reaction_id', 'target_ids', 'priority_target_ids', 'hypothesis_ids')}
    row['reference_matches'] = [{'accession': identifier, 'parent_accession': 'O95864',
        'subsequence_start_1based': start, 'subsequence_end_inclusive': end,
        'source_domain': domain['metadata']['accession'], 'source_activity': activity,
        'model_eligible': False, 'evidence_class': 'similarity-inferred-reference-domain-lead'}]
    discovery = {'schema': 'cannabis-desaturase-domain-discovery-v1', 'rows': [row],
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (protein_path, domain_path, parent_path)},
        'claim_boundary': boundary}
    source = Path('data/reports/phase1-desaturase-domain-references.json')
    source.write_text(json.dumps(discovery, separators=(',', ':')) + '\n')
    retrievals = [{'requested_accessions': [identifier], 'status': 'retrieved', 'missing_accessions': [],
        'url': 'https://rest.uniprot.org/uniprotkb/O95864.json', 'snapshot': str(protein_path),
        'sha256': hashlib.sha256(protein_path.read_bytes()).hexdigest(),
        'domain_url': 'https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00487/protein/uniprot/O95864/',
        'sequence_derivation': row['reference_matches'][0]}]
    screen(discovery, source, raw, Path('data/reports/phase1-desaturase-domain-search.json'),
           references, retrievals, evidence_class='desaturase-domain-specificity-unverified-lead',
           additional_blockers=('domain-only-homology-not-full-enzyme-evidence', 'reference-activity-inferred-by-similarity'),
           claim_boundary=boundary)


if __name__ == '__main__':
    run()
