"""Expose the separate completion sensitivity scenario in the existing net viewer."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from .phase1_net_view import build as evidence_bundle


def build(sensitivity, baseline, evidence_sources, completion_evidence, marts_audit):
    base = evidence_bundle(baseline, evidence_sources)
    certificates = base['certificates'] + sensitivity['additional_net_certificates']
    if len({c['compound_id'] for c in certificates}) != len(certificates):
        raise ValueError('Sensitivity and baseline certificates overlap')
    used = {s['reaction_id'] for c in certificates for s in c['steps']}
    reactions = {r['id']: r for r in base['reactions']}
    for r in sensitivity['certificate_baseline_reactions']:
        if r['id'] in reactions and reactions[r['id']] != r:
            raise ValueError('Baseline reaction changed')
        reactions[r['id']] = r
    ledger = {r['id']: r['source_record'] for r in marts_audit['source_ledger']}
    original_evidence = {e['id']: e for e in completion_evidence['rows']}
    additional_evidence = []
    for r in sensitivity['admitted_reactions']:
        if r['id'] not in used:
            continue
        row = original_evidence[r['candidate_evidence_record_id']]
        if row['reaction_id'] != r['id'] or not row['has_candidate_lead']:
            raise ValueError('Completion view evidence mismatch')
        eid = 'completion-sensitivity-evidence:' + row['id']
        additional_evidence.append({'id': eid, 'evidence_class': row['category'],
            'screened_proteins': [{'accession': p} for p in row['screened_cannabis_proteins']],
            'source_report': 'phase1-completion-protein-evidence.json', 'source_record_id': row['id'],
            'validation_blockers': row['validation_blockers'], 'prior_source_reviews': row['prior_source_reviews'],
            'representative_alignments': row['representative_alignments'], 'claim_boundary': row['claim_boundary']})
        sources = [{'source_record_id': sid, 'source_reaction_id': ledger[sid]['rule_id'],
                    'source_urls': [ledger[sid]['source_url']] if ledger[sid]['source_url'] else []}
                   for sid in r['original_source_record_ids']]
        reactions[r['id']] = {**r, 'enzyme_evidence_ids': [eid], 'sources': sources,
                             'is_completion_sensitivity': True}
    targets = [{**t, 'net_status': t['sensitivity_net_status'], 'startup_status': t['sensitivity_startup_status'],
        'certificate_compound_id': t['certificate']['compound_id'] if t['certificate'] else None}
        for t in sensitivity['targets']]
    compounds = {c['id']: {k: v for k, v in c.items() if k != 'labels'} for c in base['compounds']}
    for c in sensitivity['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('View compound identity conflict')
        compounds.setdefault(c['id'], c)
    # Re-resolve all baseline evidence required by newly selected upstream steps.
    baseline_count = base['summary']['target_status_counts']['exact-net-conversion-hypothesis']
    additional_count = sensitivity['summary']['additional_net_target_records']
    startup_changed = any(s['rescued_target_ids'] for s in sensitivity['startup_scenarios'])
    report = {**base, 'targets': targets, 'certificates': certificates,
        'reactions': [r for rid, r in reactions.items() if rid in used], 'compounds': list(compounds.values()),
        'summary': {'target_records': len(targets), 'target_status_counts': dict(Counter(t['net_status'] for t in targets)),
                    'additional_conditional_targets': sensitivity['summary']['additional_net_target_records']},
        'claim_boundary': sensitivity['claim_boundary'], 'view_scenario': 'completion-sensitivity',
        'view_boundary': f'Completion sensitivity: {baseline_count} baseline certificates + {additional_count} conditional additions. Orange edges use original-source homology with unverified completion chemistry. Not confirmed pathways. ' +
            ('See the separate zero-pool startup comparison.' if startup_changed else 'Zero-pool startup does not improve.')}
    return evidence_bundle(report, evidence_sources + [{'enzyme_evidence': additional_evidence}])


def run():
    names = ['phase1-completion-connectivity', 'phase1-candidate-net-flux',
        'phase1-target-hypotheses', 'phase1-screened-enzyme-overlay', 'phase1-route-enzyme-overlay',
        'phase1-completion-protein-evidence', 'phase1-marts-audit']
    paths = [Path('data/reports', n + '.json') for n in names]
    reports = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in reports:
        for p, digest in report.get('source_sha256', {}).items():
            if p in hashes and hashes[p] != digest:
                raise ValueError('View source checksum mismatch')
    bundle = build(reports[0], reports[1], reports[2:5], reports[5], reports[6])
    payload = json.dumps(bundle, separators=(',', ':')) + '\n'
    folder = Path('docs/data/completion-net-view'); folder.mkdir(parents=True, exist_ok=True)
    (folder / 'bundle.json').write_text(payload)
    manifest = {'schema': 'cannabis-carbon.phase1-completion-net-view.v1', 'file': 'bundle.json',
        'sha256': hashlib.sha256(payload.encode()).hexdigest(), 'bytes': len(payload.encode()),
        'source_sha256': hashes, 'summary': bundle['summary']}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps(manifest))


if __name__ == '__main__':
    run()
