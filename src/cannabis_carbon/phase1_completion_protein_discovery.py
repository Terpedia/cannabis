"""Exact existing-evidence joins and source-only screening queue for completions."""
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

BOUNDARY = ('Original-MARTS-reference homology leads for inferred stoichiometry, not enzyme confirmation. '
    'No proteins are inherited from composition-only Rhea templates. Exact product identity, full '
    'stoichiometry, direction and Cannabis activity remain unverified. Atom tracing is deferred.')


def build(completions, audit, scope, network, review):
    variants = {v['id']: v for v in completions['variants']}
    sources = {s['id']: s['source_record'] for s in audit['source_ledger']}
    baseline = {r['id']: r for r in network['reactions']}
    candidates = scope['candidate_reaction_evidence_ids']
    reviewed = {r['requested_accession']: r['status'] for r in review['references']}
    reviews = {r['id']: r for r in review['rows']}
    target_ids, priority = defaultdict(list), defaultdict(list)
    for t in completions['targets']:
        for hid in t['completion_ids']:
            target_ids[hid].append(t['cannabisdb_id'])
            if not t['baseline_exact_balanced_participation']:
                priority[hid].append(t['cannabisdb_id'])
    rows, existing, excluded = [], [], []
    for h in completions['completions']:
        rid = h['balanced_equation_id']
        if rid in candidates:
            r = baseline.get(rid)
            if not r or r['left'] != h['left'] or r['right'] != h['right']:
                raise ValueError('Exact equation evidence join failed side/stoichiometry validation')
            existing.append({'id': h['id'], 'reaction_id': rid, 'target_ids': target_ids[h['id']],
                'evidence_ids': candidates[rid], 'join_method': 'exact-full-canonical-equation-and-coefficients',
                'claim_boundary': 'Existing candidate evidence for the same full equation, not new activity evidence or confirmation of MARTS completion.'})
            continue
        refs = defaultdict(list); source_reviews = []
        for sid in variants[h['variant_id']]['source_record_ids']:
            s = sources[sid]; accession = s.get('source_uniprot_id')
            reason = 'missing-UniProt-reference' if not accession else 'UniParc-requires-explicit-resolution' if accession.startswith('UPI') else \
                'unsupported-accession-format' if not re.fullmatch(r'[A-Z0-9]{6,10}', accession) else \
                'known-inactive-reference-requires-resolution' if reviewed.get(accession) == 'inactive-UniProt-entry-requires-explicit-resolution' else None
            if reason:
                excluded.append({'id': h['id'] + ':' + sid, 'completion_id': h['id'], 'source_record_id': sid,
                    'source_uniprot_id': accession, 'source_genbank_id': s.get('source_genbank_id'), 'status': reason})
            else:
                refs[accession].append(sid)
            if sid in reviews:
                source_reviews.append({'source_record_id': sid, 'status': reviews[sid]['status'],
                    'reference_status': reviews[sid]['reference_status'],
                    'stereo_only_lead': any(a['stereo_only_diagnostic_lead'] for m in reviews[sid]['catalytic_activity_reviews'] for a in m['equation_assessments'])})
        rows.append({'reaction_id': rid, 'hypothesis_ids': [h['id']], 'target_ids': target_ids[h['id']],
            'priority_target_ids': priority[h['id']], 'left': h['left'], 'right': h['right'],
            'marts_forward_direction': h['marts_forward_direction'],
            'source_record_ids': variants[h['variant_id']]['source_record_ids'],
            'reference_matches': [{'accession': acc, 'source_record_ids': ids,
                'join_method': 'original-MARTS-transformation-source-reference; not a cofactor-template enzyme',
                'claim_boundary': BOUNDARY} for acc, ids in sorted(refs.items())],
            'prior_source_reviews': source_reviews, 'claim_boundary': BOUNDARY})
    if len({r['reaction_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate reaction in search queue')
    return {'schema': 'cannabis-carbon.phase1-completion-protein-discovery.v1', 'rows': rows,
        'existing_evidence_matches': existing, 'excluded_source_references': excluded,
        'summary': {'completions': len(completions['completions']), 'existing_exact_equation_evidence': len(existing),
            'queued_equations': len(rows), 'queued_UniProt_identifiers': len({m['accession'] for r in rows for m in r['reference_matches']}),
            'queued_target_ids': len({t for r in rows for t in r['target_ids']}),
            'queued_priority_target_ids': len({t for r in rows for t in r['priority_target_ids']}),
            'excluded_reference_status_counts': dict(Counter(r['status'] for r in excluded))},
        'claim_boundary': BOUNDARY}


def run():
    names = ['phase1-marts-completions', 'phase1-marts-audit', 'phase1-candidate-scope',
             'phase1-full-balanced-network', 'phase1-marts-gap-references']
    paths = [Path('data/reports', n + '.json') for n in names]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    output = Path('data/reports/phase1-completion-protein-discovery.json')
    output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
