"""Complete carrier-separated candidate network; not validated pathway coverage."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from .phase1_catalog import stable_id
from .phase1_glycerophospholipid_input_audit import assemble_current


def run():
    root = Path('data/reports')
    current_path = root / 'phase1-geranial-reduction-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [root / ('phase1-' + n + '.json') for n in
             ('full-balanced-network', 'marts-completions', 'carrier-reconstruction', 'carrier-scaffolds')]
    paths += [current_path] + [root / n for n in current['baseline_certificate_reports']]
    network, completions, restored, scaffold, current, *layers = [json.loads(p.read_bytes()) for p in paths]
    for doc in (restored, scaffold):
        for p, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale carrier input')
    originals, compounds, old_model = assemble_current(network, completions, current, layers)
    replaced = {j['model_reaction_id'] for r in restored['restored_equations'] for j in r['source_joins']}
    audits = {r['id']: r for r in scaffold['equations']}
    reactions = [{**r, 'carrier_identity_status': 'no-direct-carrier-replacement-not-full-identity-validation'}
                 for rid, r in originals.items() if rid not in replaced]
    old_steps = {}
    for s in old_model.steps:
        old_steps.setdefault(s['reaction_id'], []).append(s)
    steps = [s for s in old_model.steps if s['reaction_id'] not in replaced]
    for r in restored['restored_equations']:
        reactions.append({**r, 'scaffold_audit': audits[r['id']],
                          'carrier_identity_status': 'source-carrier-identities-restored'})
        # Preserve the historical permitted orientation separately from source
        # direction evidence. Never re-enable a forbidden original step.
        by_direction = {}
        for j in r['source_joins']:
            oid = j['model_reaction_id']
            for old in old_steps.get(oid, []):
                if old['required_inputs'] == originals[oid]['left'] and old['outputs'] == originals[oid]['right']:
                    side = 'left'
                elif old['required_inputs'] == originals[oid]['right'] and old['outputs'] == originals[oid]['left']:
                    side = 'right'
                else:
                    raise ValueError('Unknown original step orientation')
                by_direction.setdefault(side, {})[old['id']] = old
        for side, parents in by_direction.items():
            steps.append({'id': stable_id('carrier-separated-step', [r['id'], side]),
                          'reaction_id': r['id'], 'required_inputs': r[side],
                          'outputs': r['right' if side == 'left' else 'left'],
                          'direction_mode': 'hypothetical-' + side + '-to-' + ('right' if side == 'left' else 'left'),
                          'inherited_allowed_step_ids': sorted(parents),
                          'direction_evidence_boundary': 'Inherited sensitivity-model permission, not physiological direction. Source joins remain separate.',
                          'scaffold_status': audits[r['id']]['status']})
    carriers = {c['id']: c for c in restored['carriers']}
    used = {p['compound_id'] for r in reactions for side in ('left', 'right') for p in r[side]}
    for cid in used & carriers.keys():
        c = carriers[cid]
        name = next(p['object']['value'] for p in c['properties'] if p['predicate'].endswith('/name'))
        compounds[cid] = {'id': cid, 'label': name, 'source_carrier': c,
                          'full_formula': None, 'full_charge': None, 'carbon_count': None,
                          'external_uptake_allowed': False,
                          'identity_status': 'source-defined-macromolecular-carrier'}
    exchange = current['external_exchange_compound_ids']
    if set(exchange) & carriers.keys():
        raise ValueError('Carrier uptake forbidden')
    if not used <= compounds.keys():
        raise ValueError('Missing participant identity')
    producers = set()
    for s in steps:
        net = Counter()
        for side, sign in (('required_inputs', -1), ('outputs', 1)):
            for p in s[side]:
                net[p['compound_id']] += sign * p['coefficient']
        producers.update(c for c, n in net.items() if n > 0)
    targets = [{**{k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id')},
                'candidate_participant': t['compound_id'] in used,
                'has_candidate_net_producer': t['compound_id'] in producers,
                'qualified_co2_pathway_status': 'not-recomputed-after-carrier-separation'} for t in current['targets']]
    report = {'schema': 'cannabis-carbon.phase1-carrier-network.v1',
              'reactions': reactions, 'directed_steps': steps, 'compounds': list(compounds.values()),
              'targets': targets, 'replaced_projected_reaction_ids': sorted(replaced),
              'quarantined_original_equations': [originals[rid] for rid in sorted(replaced)],
              'external_exchange_compound_ids': exchange, 'co2_compound_id': current['co2_compound_id'],
              'summary': {'original_projected_equations': len(originals), 'replaced_projected_equations': len(replaced),
                          'carrier_separated_candidate_equations': len(reactions), 'directed_candidate_steps': len(steps),
                          'distinct_internal_source_carriers': len(used & carriers.keys()), 'inventory_records': len(targets),
                          'candidate_participating_records': sum(t['candidate_participant'] for t in targets),
                          'records_without_candidate_net_producer': sum(not t['has_candidate_net_producer'] for t in targets),
                          'qualified_co2_pathway_coverage': None},
              'claim_boundary': 'Candidate network reconstruction, not fully balanced or biologically confirmed pathways. '
                'All directly matched projected carrier equations are replaced; originals remain quarantined for provenance, including any mixed-source alternatives. '
                'Nonmatched catalog equations and inferred templates still need full carrier/polymer auditing. '
                'Carrier species are internal only, with unknown complete formulas and no automatic startup seed. '
                'Scaffold hypotheses do not establish full balance. No historical certificate is reused or promoted. '
                'Inherited external exchanges remain a permissive sensitivity boundary, not a minimum defined medium.',
              'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-carrier-network.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
