"""Current full-model reaction-gap queue; diagnostic fingerprints never add edges."""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from rdkit import RDLogger, rdBase
from .phase1_glycerophospholipid_input_audit import assemble_current
from .phase1_no_producer_audit import diagnostic_keys
from .phase1_scope import orientations

BOUNDARY = ('Current-model diagnostic, not additional reaction or pathway coverage. '
    'All historical identities and direction exclusions are retained. Positive net '
    'production in one equation does not establish supply of all inputs. Stereo removal '
    'and uncharging are search fingerprints only, not identity equivalence, reaction '
    'chemistry or permission to merge. Alternative compounds are separate assertions. '
    'Source conflicts take priority over chemical gap filling; absence from this model '
    'is not biological absence. No solver result is changed or inferred from these joins.')


def producing_steps(reactions, forbidden=()):
    result = defaultdict(set)
    forbidden = set(forbidden)
    for step in orientations(reactions):
        if step['id'] in forbidden:
            continue
        net = defaultdict(Fraction)
        for sign, side in [(-1, 'required_inputs'), (1, 'outputs')]:
            for p in step[side]:
                net[p['compound_id']] += sign * Fraction(p['coefficient'])
        for cid, amount in net.items():
            if amount > 0:
                result[cid].add(step['id'])
    return result


def build(current, compounds, reactions, source_audit):
    allowed = producing_steps(list(reactions.values()), current['forbidden_step_ids'])
    unrestricted = producing_steps(list(reactions.values()))
    conflicts = {r['cannabisdb_id'] for r in source_audit['source_identity_conflicts']}
    source_rows = {r['cannabisdb_id']: r for r in source_audit['targets']}
    indexes = {k: defaultdict(set) for k in ('stereo_removed', 'uncharger', 'uncharger_and_stereo_removed')}
    for cid in sorted(allowed):
        for kind, value in diagnostic_keys(compounds[cid]['smiles']).items():
            indexes[kind][value].add(cid)
    rows = []; used = set()
    for t in current['targets']:
        if t['net_status'] != 'no-net-producing-equation':
            continue
        cid = t['compound_id']
        if allowed[cid]:
            raise ValueError('Stored no-producer status conflicts with current directions')
        alternatives = {kind: sorted(indexes[kind].get(value, set()) - {cid})
                        for kind, value in diagnostic_keys(compounds[cid]['smiles']).items()}
        category = ('source-identity-conflict' if t['cannabisdb_id'] in conflicts else
                    'producer-excluded-by-direction' if unrestricted[cid] else
                    'diagnostic-identity-leads-only' if any(alternatives.values()) else
                    'no-exact-or-diagnostic-producing-match')
        source = source_rows.get(t['cannabisdb_id'], {})
        rows.append({k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'net_status')} | {
            'canonical_smiles': compounds[cid]['smiles'], 'category': category,
            'source_url': source.get('source_url'),
            'source_external_ids': source.get('source_external_ids', {}),
            'excluded_producing_step_ids': sorted(unrestricted[cid]),
            'diagnostic_alternatives': alternatives,
            'next_action': {
                'source-identity-conflict': 'Resolve source structure assertions before adding chemistry.',
                'producer-excluded-by-direction': 'Review excluded direction evidence; do not relax bounds automatically.',
                'diagnostic-identity-leads-only': 'Review exact charge and stereochemistry differences; no diagnostic match is a reaction.',
                'no-exact-or-diagnostic-producing-match': 'Find a source-backed exact producing equation and then audit every required input.'}[category]})
        used.add(cid)
        used.update(a for values in alternatives.values() for a in values)
    target_statuses = defaultdict(set)
    for t in current['targets']:
        target_statuses[t['compound_id']].add(t['net_status'])
    alternatives = sorted({a for row in rows for values in row['diagnostic_alternatives'].values() for a in values})
    used_steps = {s for cid in alternatives for s in allowed[cid]} | {s for row in rows for s in row['excluded_producing_step_ids']}
    used_reactions = {s.rsplit(':', 1)[0] for s in used_steps}
    selected = [r for r in reactions.values() if r['id'] in used_reactions]
    used.update(p['compound_id'] for r in selected for side in ('left', 'right') for p in r[side])
    return {'schema': 'cannabis-carbon.phase1-current-reaction-gaps.v1', 'targets': rows,
        'alternative_producers': [{'compound_id': cid, 'allowed_producing_step_ids': sorted(allowed[cid]),
                                  'historical_target_statuses': sorted(target_statuses[cid])} for cid in alternatives],
        'compounds': [compounds[cid] for cid in sorted(used)], 'reactions': selected,
        'forbidden_step_ids': current['forbidden_step_ids'],
        'historical_inventory_summary_unchanged': current['summary'], 'claim_boundary': BOUNDARY,
        'rdkit_version': rdBase.rdkitVersion,
        'summary': {'inventory_records': len(current['targets']), 'balanced_equations': len(reactions),
            'no_producer_records': len(rows), 'category_counts': dict(Counter(r['category'] for r in rows)),
            'diagnostic_counts': {k: sum(bool(r['diagnostic_alternatives'][k]) for r in rows) for k in indexes},
            'alternative_structures': len(alternatives), 'coverage_gain_claimed': 0}}


def run():
    RDLogger.DisableLog('rdApp.warning')
    root = Path('data/reports')
    current = json.loads((root/'phase1-glycerophospholipid-net.json').read_bytes())
    paths = [root/'phase1-full-balanced-network.json', root/'phase1-marts-completions.json',
             root/'phase1-glycerophospholipid-net.json', root/'phase1-no-producer-audit.json'] + [
                 root/n for n in current['baseline_certificate_reports']]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale gap source')
    network, completions, current, old, *layers = docs
    reactions, compounds, _ = assemble_current(network, completions, current, layers)
    report = build(current, compounds, reactions, old)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (root/'phase1-current-reaction-gaps.json').write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
