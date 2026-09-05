"""Full balanced-network qualitative scope, with explicit carbon-source seeds."""
import copy
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from rdkit import RDLogger
from .balance import _reaction_smiles_balance
from .phase1_balance_reference import concrete_participants
from .phase1_catalog import stable_id
from .phase1_screened_overlay import apply_overlay
from .phase1_target_rhea_coverage import source_records


def extend_network(parent, catalog):
    compounds = {c['id']: copy.deepcopy(c) for c in parent['compounds']}
    reactions = {r['id']: copy.deepcopy(r) for r in parent['reactions']}
    known_sources = {s['coverage_record_id'] for r in reactions.values() for s in r['sources']}
    counts, excluded = Counter(), []
    for source in source_records(catalog):
        smiles = source['reaction_smiles']
        element, charge = _reaction_smiles_balance(smiles)
        status = 'balanced' if element and charge and element['status'] == charge['status'] == 'balanced' else 'imbalanced' if element and charge else 'not-auditable'
        counts[status] += 1
        if status != 'balanced':
            excluded.append({'source_reaction_id': source['source_reaction_id'], 'status': status,
                'source_url': source['source_url'], 'element_balance': element, 'charge_balance': charge})
            continue
        sides = []
        for participants in concrete_participants(smiles):
            members = []
            for p in participants:
                cid = stable_id('structure', p['smiles'])
                compounds.setdefault(cid, {'id': cid, **{k: v for k, v in p.items() if k != 'coefficient'}})
                members.append({'compound_id': cid, 'coefficient': p['coefficient']})
            sides.append(sorted(members, key=lambda m: m['compound_id']))
        flipped = json.dumps(sides[0], sort_keys=True) > json.dumps(sides[1], sort_keys=True)
        if flipped:
            sides.reverse()
        rid = stable_id('balanced-equation', sides)
        r = reactions.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1], 'sources': [],
            'balance_status': 'independently-balanced', 'direction_status': 'unresolved-in-Cannabis; canonical side ordering is not physiology'})
        if source['id'] not in known_sources:
            r['sources'].append({'coverage_record_id': source['id'], 'source_reaction_id': source['source_reaction_id'],
                'source_layer': source['source_layer'], 'source_urls': [source['source_url']],
                'source_left_corresponds_to': 'right' if flipped else 'left', 'evidence_ids': []})
            known_sources.add(source['id'])
    evidence = defaultdict(set)
    for h in parent['hypotheses']:
        evidence[h['reaction_id']].update(h['evidence_ids'])
    for r in reactions.values():
        r['enzyme_evidence_ids'] = sorted(evidence[r['id']])
    return {'schema': 'cannabis-carbon.phase1-full-balanced-network.v1',
        'summary': {'compound_structures': len(compounds), 'balanced_equations': len(reactions),
            'parent_balanced_equations': len(parent['reactions']), 'added_upstream_equations': len(reactions) - len(parent['reactions']),
            'rhea_source_rows': len(catalog), 'rhea_balance_status_counts': dict(counts)},
        'compounds': list(compounds.values()), 'reactions': list(reactions.values()),
        'targets': parent['targets'], 'excluded_rhea_source_records': excluded,
        'enzyme_evidence_sources': ['phase1-target-hypotheses.json', 'phase1-screened-enzyme-overlay.json'],
        'claim_boundary': 'Balanced full-source equations, not all confirmed Cannabis reactions. Exact structures and all coefficients are retained. Source orientations and candidate evidence do not establish physiological direction or enzyme activity.'}


def orientations(reactions):
    rows = []
    for reaction in sorted(reactions, key=lambda r: r['id']):
        for direction, inputs, outputs in [('hypothetical-left-to-right', reaction['left'], reaction['right']),
                                           ('hypothetical-right-to-left', reaction['right'], reaction['left'])]:
            if not inputs or not outputs or any(m['coefficient'] <= 0 for m in inputs + outputs):
                raise ValueError('Scope requires nonempty sides and positive coefficients')
            rows.append({'id': reaction['id'] + ':' + direction, 'reaction_id': reaction['id'],
                'direction_mode': direction, 'required_inputs': inputs, 'outputs': outputs})
    return rows


def expand(reactions, seeds):
    available = set(seeds)
    witnesses = {cid: {'level': 0, 'source': 'explicit-seed'} for cid in available}
    pending = orientations(reactions)
    enabled = []
    level = 0
    while True:
        level += 1
        new, remaining = {}, []
        for step in pending:
            if {m['compound_id'] for m in step['required_inputs']} <= available:
                enabled.append(step['id'])
                for product in step['outputs']:
                    cid = product['compound_id']
                    if cid not in available:
                        new.setdefault(cid, {'level': level, 'source': 'all-reactants-scope', **step})
            else:
                remaining.append(step)
        if not new:
            break
        witnesses.update(new); available.update(new); pending = remaining
    blocked = [{**step, 'missing_inputs': [m for m in step['required_inputs'] if m['compound_id'] not in available]} for step in remaining]
    return {'available': available, 'witnesses': witnesses, 'enabled_orientations': enabled, 'blocked_orientations': blocked}


def assess(network):
    compounds = {c['id']: c for c in network['compounds']}
    co2 = stable_id('structure', 'O=C=O')
    if co2 not in compounds or compounds[co2]['carbon_count'] != 1:
        raise ValueError('Exact CO2 seed missing from network')
    scenarios = []
    for name, seeds, boundary in [
        ('CO2-only', {co2}, 'Only exact CO2 is seeded; no other nutrients, energy carriers or cofactors are supplied.'),
        ('CO2-plus-all-carbon-free-species', {co2} | {cid for cid, c in compounds.items() if c['carbon_count'] == 0},
         'Highly permissive diagnostic reservoir: every cataloged carbon-free species is seeded, including unusual ions and redox species. This is not a plant growth medium or physiological claim. No other carbon-containing compound is seeded.')]:
        result = expand(network['reactions'], seeds)
        producers = defaultdict(list)
        for step in result['blocked_orientations']:
            input_counts = {m['compound_id']: m['coefficient'] for m in step['required_inputs']}
            for output in step['outputs']:
                if output['coefficient'] > input_counts.get(output['compound_id'], 0):
                    producers[output['compound_id']].append({'reaction_id': step['reaction_id'],
                        'direction_mode': step['direction_mode'], 'missing_inputs': step['missing_inputs']})
        targets = [{'cannabisdb_id': t['cannabisdb_id'], 'label': t['label'], 'compound_id': t['compound_id'],
                    'carbon_count': t['carbon_count'],
                    'status': 'explicit-seed' if t['compound_id'] in seeds else 'structural-scope-reachable' if t['compound_id'] in result['available'] else 'blocked',
                    'blocked_producing_steps': producers.get(t['compound_id'], [])} for t in network['targets']]
        scenarios.append({'id': name, 'seed_compound_ids': sorted(seeds), 'seed_boundary': boundary,
            'summary': {'seed_count': len(seeds), 'scope_compounds': len(result['available']),
                'newly_available_compounds': len(result['available'] - seeds),
                'target_status_counts': dict(Counter(t['status'] for t in targets)),
                'carbon_bearing_target_status_counts': dict(Counter(t['status'] for t in targets if t['carbon_count'])),
                'reactant_denominator': len({m['compound_id'] for r in network['reactions'] for side in ('left', 'right') for m in r[side]}),
                'enabled_hypothetical_orientations': len(result['enabled_orientations'])},
            'targets': targets, 'witnesses': result['witnesses'], 'enabled_orientations': result['enabled_orientations'],
            'direction_boundary': 'Every reaction is examined in both directions as a structural upper-bound scenario only. No curated-direction or physiological reversibility claim is made.'})
    return {'schema': 'cannabis-carbon.phase1-all-reactants-scope.v1', 'scenarios': scenarios,
        'claim_boundary': 'Qualitative compound availability, not stoichiometric flux feasibility, thermodynamics, enzyme activity, tissue compatibility or traced carbon provenance. Every reactant is required; finite pool sizes and consumption are not modeled. Full equations and coefficients are preserved. Atom tracing remains deferred.'}


def run():
    RDLogger.DisableLog('rdApp.warning'); RDLogger.DisableLog('rdApp.error')
    parent_path = Path('data/reports/phase1-target-hypotheses.json')
    overlay_path = Path('data/reports/phase1-screened-enzyme-overlay.json')
    metadata_path = Path('data/reports/phase1-balance-reference.json')
    metadata = json.loads(metadata_path.read_text())['catalog']
    raw_path = Path(metadata['snapshot'])
    if hashlib.sha256(raw_path.read_bytes()).hexdigest() != metadata['sha256']:
        raise ValueError('Full Rhea source checksum mismatch')
    overlay = json.loads(overlay_path.read_text())
    if hashlib.sha256(parent_path.read_bytes()).hexdigest() != overlay['source_sha256'][str(parent_path)]:
        raise ValueError('Overlay parent checksum mismatch')
    parent = apply_overlay(json.loads(parent_path.read_text()), overlay)
    network = extend_network(parent, json.loads(raw_path.read_text()))
    network['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [parent_path, overlay_path, raw_path, metadata_path]}
    network['rhea_catalog_provenance'] = metadata
    network_path = Path('data/reports/phase1-full-balanced-network.json')
    network_path.write_text(json.dumps(network, separators=(',', ':')) + '\n')
    report = assess(network)
    report['source_network'] = str(network_path)
    report['source_network_sha256'] = hashlib.sha256(network_path.read_bytes()).hexdigest()
    Path('data/reports/phase1-all-reactants-scope.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({'network': network['summary'], 'scenarios': {s['id']: s['summary'] for s in report['scenarios']}}))


def export_network(report_path, output):
    raw = report_path.read_bytes(); report = json.loads(raw)
    rows = []
    for kind, collection, key in [('compound', 'compounds', 'id'), ('reaction', 'reactions', 'id'),
                                   ('target', 'targets', 'cannabisdb_id'), ('excluded_source', 'excluded_rhea_source_records', 'source_reaction_id')]:
        rows.extend((kind, row[key], row) for row in report[collection])
    rows.append(('metadata', 'network', {k: v for k, v in report.items() if k not in ['compounds', 'reactions', 'targets', 'excluded_rhea_source_records']}))
    return write_rows(rows, hashlib.sha256(raw).hexdigest(), output)


def export_scope(report_path, output):
    raw = report_path.read_bytes(); report = json.loads(raw)
    rows = [('metadata', 'scope', {k: v for k, v in report.items() if k != 'scenarios'})]
    for scenario in report['scenarios']:
        sid = scenario['id']
        rows.append(('scenario', sid, {k: v for k, v in scenario.items() if k not in ['targets', 'witnesses']}))
        rows.extend(('target', sid + ':' + t['cannabisdb_id'], {'scenario_id': sid, **t}) for t in scenario['targets'])
        rows.extend(('witness', sid + ':' + cid, {'scenario_id': sid, 'compound_id': cid, **w}) for cid, w in scenario['witnesses'].items())
    return write_rows(rows, hashlib.sha256(raw).hexdigest(), output)


def write_rows(rows, digest, output):
    with output.open('w') as stream:
        for kind, identifier, row in rows:
            stream.write(json.dumps({'record_kind': kind, 'record_id': identifier,
                'record_json': json.dumps(row, separators=(',', ':')), 'report_sha256': digest}, separators=(',', ':')) + '\n')
    return len(rows)


if __name__ == '__main__':
    run()
