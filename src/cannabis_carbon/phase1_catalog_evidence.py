"""Compact, exact-reaction evidence supplement; frozen net chemistry is unchanged."""
import hashlib
import json
from pathlib import Path

from .phase1_screened_overlay import build_overlay
from .phase1_scope import write_rows


def build(catalog, search):
    added = build_overlay({'reactions': catalog['reactions'], 'hypotheses': []}, search,
        'https://github.com/Terpedia/cannabis/blob/main/data/reports/phase1-catalog-protein-search.json')
    by_reaction = {e['reaction_id']: e for e in added['enzyme_evidence']}
    reactions = {r['id']: r for r in catalog['reactions']}
    existing = {e['id'] for e in catalog['enzyme_evidence']}
    for evidence in added['enzyme_evidence']:
        if evidence['id'] in existing or reactions[evidence['reaction_id']]['enzyme_evidence_ids']:
            raise ValueError('Supplement overlaps existing candidate evidence')
        evidence['evidence_class'] = 'catalog-net-gap-direction-unresolved-reference-homology'
    updates = []
    for cert in catalog['certificates']:
        before = cert['missing_candidate_reaction_ids']
        after = [rid for rid in before if rid not in by_reaction]
        if len(after) != len(before):
            updates.append({'compound_id': cert['compound_id'],
                'baseline_missing_candidate_reaction_ids': before,
                'missing_candidate_reaction_ids': after})
    by_compound = {c['compound_id']: c for c in updates}
    targets = [{ 'cannabisdb_id': t['cannabisdb_id'], **by_compound[t['compound_id']]}
        for t in catalog['targets'] if t['compound_id'] in by_compound]
    remaining = {g['id'] for g in catalog['gap_priorities']} - by_reaction.keys()
    newly_closed = [t['cannabisdb_id'] for t in targets if not t['missing_candidate_reaction_ids']]
    all_closed = sum(bool(t['certificate_compound_id']) and not
        by_compound.get(t['compound_id'], t)['missing_candidate_reaction_ids'] for t in catalog['targets'])
    summary = {'added_candidate_equations': len(by_reaction),
        'distinct_new_candidate_proteins': added['summary']['distinct_cannabis_proteins'],
        'baseline_missing_candidate_equations': len(catalog['gap_priorities']),
        'remaining_missing_candidate_equations': len(remaining),
        'target_records_with_added_evidence': len(targets),
        'newly_candidate_linked_selected_certificate_target_ids': newly_closed,
        'selected_certificate_targets_with_candidates_for_all_steps': all_closed,
        'selected_certificate_targets_with_remaining_gaps': catalog['summary']['target_status_counts']['exact-net-conversion-hypothesis'] - all_closed}
    return {'schema': 'cannabis-carbon.phase1-catalog-evidence.v1',
        'enzyme_evidence': added['enzyme_evidence'], 'certificate_updates': updates,
        'target_updates': targets, 'summary': summary,
        'view_boundary': f'Balanced-catalog diagnostic: 304 target records have exact net balances; {all_closed} selected certificates now have candidates for all steps, not confirmed activity. Red edges still lack candidates; blue edges show newly screened candidates. Net chemistry and baseline snapshots are unchanged.',
        'claim_boundary': 'Evidence-only supplement to frozen chemistry-only net certificates. No equation, extent, target identity, exchange, direction or startup result changes. Candidate proteins are not proven enzyme activity or specificity; internal-pool origin and physiological feasibility remain unresolved. No net feasibility is recomputed and no necessity is inferred. Atom tracing remains deferred.'}


def run():
    paths = [Path('data/reports', n + '.json') for n in (
        'phase1-catalog-net-gaps', 'phase1-catalog-protein-search', 'phase1-catalog-references')]
    catalog, search, references = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if search['source_discovery_sha256'] != hashes[str(paths[2])] or references['source_sha256'][str(paths[0])] != hashes[str(paths[0])]:
        raise ValueError('Evidence lineage mismatch')
    output = build(catalog, search); output['source_sha256'] = hashes
    payload = json.dumps(output, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-catalog-evidence.json').write_text(payload)
    groups = [('enzyme_evidence', 'enzyme_evidence', 'id'), ('certificate_update', 'certificate_updates', 'compound_id'),
        ('target_update', 'target_updates', 'cannabisdb_id')]
    records = [('metadata', 'supplement', {k:v for k,v in output.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        records.extend((kind, row[key], row) for row in output[collection])
    sha = hashlib.sha256(payload.encode()).hexdigest()
    count = write_rows(records, sha, Path('data/derived/phase1-catalog-evidence.ndjson'))
    print(json.dumps({'summary': output['summary'], 'bytes': len(payload.encode()), 'rows': count, 'sha256': sha}))


if __name__ == '__main__':
    run()
