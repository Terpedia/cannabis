"""Separate CO2 sensitivity scenario for source-homology completion hypotheses."""
import copy
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
import scipy
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel
from .phase1_scope import expand, write_rows


BOUNDARY = ('Sensitivity scenario only, not a promotion to the baseline network. '
    'Original-source protein homology does not validate product identity, inorganic '
    'stoichiometry, Cannabis activity or direction. All equations are considered in '
    'both hypothetical directions. CO2 is the sole carbon exchange; the unchanged '
    'carbon-free reservoir is not a physiological medium. Net certificates permit '
    'regenerated pre-existing internal pools and do not establish their origin, '
    'energy feasibility, compartment compatibility or zero-pool startup. '
    'Atom tracing remains deferred.')


def assemble(network, startup, completions, overlay):
    """Validate exact joins and build copies; never mutate source chemistry/evidence."""
    compounds = {c['id']: copy.deepcopy(c) for c in network['compounds']}
    for c in completions['compounds']:
        if c['id'] in compounds and compounds[c['id']]['smiles'] != c['smiles']:
            raise ValueError('Compound identity conflict')
        if stable_id('structure', c['smiles']) != c['id']:
            raise ValueError('Compound identity digest mismatch')
        compounds.setdefault(c['id'], copy.deepcopy(c))
    original = {r['id']: r for r in network['reactions']}
    baseline = startup['candidate_reaction_evidence_ids']
    if not baseline.keys() <= original.keys() or any(not ids for ids in baseline.values()):
        raise ValueError('Invalid baseline candidate evidence')
    reactions = {rid: {**copy.deepcopy(original[rid]), 'enzyme_evidence_ids': list(ids)}
                 for rid, ids in baseline.items()}
    evidence = {r['id']: r for r in overlay['rows']}
    if len(evidence) != len(overlay['rows']) or evidence.keys() != {h['id'] for h in completions['completions']}:
        raise ValueError('Completion evidence inventory mismatch')
    variants = {v['id']: v for v in completions['variants']}
    admitted, excluded = [], []
    for h in completions['completions']:
        e = evidence[h['id']]; rid = h['balanced_equation_id']
        if e['reaction_id'] != rid or stable_id('balanced-equation', [h['left'], h['right']]) != rid:
            raise ValueError('Completion equation identity mismatch')
        if not balanced([h['left'], h['right']], compounds):
            raise ValueError('Completion is not element/isotope/charge balanced')
        if rid in original and any(h[side] != original[rid][side] for side in ('left', 'right')):
            raise ValueError('Full equation join mismatch')
        if not e['has_candidate_lead']:
            excluded.append({'id': h['id'], 'reaction_id': rid, 'reason': 'no-screened-candidate-lead'})
            continue
        if e['category'] == 'existing-exact-equation-candidate-evidence':
            if rid not in baseline or set(e['existing_evidence_ids']) != set(baseline[rid]):
                raise ValueError('Existing candidate evidence mismatch')
            continue
        if e['category'] != 'MARTS-source-homology-for-inferred-stoichiometry' or not e['screened_cannabis_proteins']:
            raise ValueError('Invalid completion homology category')
        if rid in reactions:
            raise ValueError('New completion already admitted')
        r = {'id': rid, 'left': copy.deepcopy(h['left']), 'right': copy.deepcopy(h['right']),
            'balance_status': 'independently-balanced; inferred-completion-sensitivity',
            'direction_status': 'both-directions-hypothetical; not physiological',
            'baseline_balanced_equation_exists': rid in original,
            'completion_id': h['id'], 'variant_id': h['variant_id'],
            'original_source_record_ids': variants[h['variant_id']]['source_record_ids'],
            'candidate_evidence_report': 'phase1-completion-protein-evidence.json',
            'candidate_evidence_record_id': e['id'],
            'screened_cannabis_proteins': e['screened_cannabis_proteins'],
            'candidate_evidence_class': e['category'],
            'validation_blockers': e['validation_blockers'],
            'prior_source_reviews': e['prior_source_reviews'],
            'claim_boundary': BOUNDARY}
        reactions[rid] = r; admitted.append(r)
    return list(reactions.values()), compounds, admitted, excluded


def build(network, startup, completions, overlay, baseline_net):
    reactions, compounds, admitted, excluded = assemble(network, startup, completions, overlay)
    scenarios = []
    for old in startup['scenarios']:
        seeds = set(old['seed_compound_ids'])
        result = expand(reactions, seeds)
        old_available = set(old['witnesses'])
        if not old_available <= result['available']:
            raise ValueError('Additive startup scenario lost baseline compounds')
        scenarios.append({'id': old['id'], 'seed_compound_ids': sorted(seeds),
            'seed_boundary': old['seed_boundary'], 'direction_boundary': old['direction_boundary'],
            'newly_available_vs_baseline': sorted(result['available'] - old_available),
            'rescued_target_ids': [t['cannabisdb_id'] for t in old['targets']
                if t['compound_id'] in result['available'] - old_available],
            'witnesses': result['witnesses']})
    exchange = set(baseline_net['external_exchange_compound_ids'])
    scenario = next(s for s in scenarios if s['id'] == 'CO2-plus-all-carbon-free-species')
    if exchange != set(scenario['seed_compound_ids']):
        raise ValueError('External input boundary changed')
    co2 = baseline_net['co2_compound_id']
    if [compounds[c]['smiles'] for c in exchange if compounds[c]['carbon_count']] != ['O=C=O']:
        raise ValueError('CO2 must be the sole carbon input')
    if compounds[co2]['smiles'] != 'O=C=O' or co2 not in exchange:
        raise ValueError('Invalid CO2 identity')
    if [t['cannabisdb_id'] for t in baseline_net['targets']] != [t['cannabisdb_id'] for t in network['targets']]:
        raise ValueError('Baseline target inventory changed')
    model = NetModel(reactions, exchange)
    baseline_certs = {c['compound_id']: c for c in baseline_net['certificates']}
    cache, new_certs, targets = {}, {}, []
    for i, t in enumerate(baseline_net['targets'], 1):
        cid = t['compound_id']
        if cid in baseline_certs:
            # Addition cannot invalidate an existing exact certificate. Preserve
            # that original witness instead of presenting solver churn as progress.
            result = baseline_certs[cid]
            certificate = {'report': 'phase1-candidate-net-flux.json', 'compound_id': cid}
        else:
            if cid not in cache:
                cache[cid] = model.solve(cid)
            result = cache[cid]; certificate = None
            if result['status'] == 'exact-net-conversion-hypothesis':
                carbon_in = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['external_net_consumption'].items())
                carbon_out = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['net_exports'].items())
                if carbon_in <= 0 or carbon_in != carbon_out or Fraction(result['external_net_consumption'].get(co2, '0')) != carbon_in:
                    raise ValueError('Exact CO2 carbon balance failed')
                new_certs[cid] = {'compound_id': cid, **result,
                    'net_carbon_in': str(carbon_in), 'net_carbon_out': str(carbon_out)}
                certificate = {'report': 'phase1-completion-connectivity.json', 'compound_id': cid}
        targets.append({k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'carbon_count')} |
            {'baseline_net_status': t['net_status'], 'sensitivity_net_status': result['status'],
             'baseline_startup_status': t['startup_status'],
             'sensitivity_startup_status': 'explicit-seed' if cid in exchange else
                'hypothesis-scope-reachable' if cid in scenario['witnesses'] else 'blocked',
             'certificate': certificate} |
            {k: result[k] for k in ('solver_status', 'solver_message') if k in result})
        if i % 500 == 0:
            print(f'Completion connectivity: {i}/{len(baseline_net["targets"])} targets, {len(new_certs)} additional certificates', flush=True)
    admitted_ids = {r['id'] for r in admitted}
    used_ids = {s['reaction_id'] for c in new_certs.values() for s in c['steps']}
    used_reactions = [r for r in reactions if r['id'] in used_ids - admitted_ids]
    for r in used_reactions:
        if not balanced([r['left'], r['right']], compounds):
            raise ValueError('Certificate baseline equation failed independent balance')
    for cert in new_certs.values():
        cert['completion_dependencies'] = [r['completion_id'] for r in admitted
            if r['id'] in {s['reaction_id'] for s in cert['steps']}]
        cert['proposed_tests'] = [
            'Resolve the original source product against the exact CannabisDB structure, including absolute stereochemistry; review any retained source-product warnings before assigning specificity.',
            'Express the linked Cannabis candidate proteins and assay the exact original organic substrate. Identify all products against authentic standards and quantify the inferred inorganic coproducts; homology alone is not a positive result.',
            'For a successful enzyme assay, establish compartment compatibility, substrate/cofactor availability, expression and a viable reaction direction for the complete upstream route. This net certificate does not establish startup or pool origin.']
    keep = exchange | {m['compound_id'] for r in admitted + used_reactions for side in ('left', 'right') for m in r[side]}
    keep.update(c for s in scenarios for c in s['witnesses'])
    return {'schema': 'cannabis-carbon.phase1-completion-connectivity.v1',
        'summary': {'target_records': len(targets), 'baseline_candidate_equations': len(startup['candidate_reaction_evidence_ids']),
            'admitted_completion_equations': len(admitted), 'sensitivity_equations': len(reactions),
            'admitted_equations_already_balanced_in_baseline': sum(r['baseline_balanced_equation_exists'] for r in admitted),
            'additional_inferred_equations': sum(not r['baseline_balanced_equation_exists'] for r in admitted),
            'excluded_completions_without_candidate_lead': len(excluded),
            'additional_exact_structure_net_certificates': len(new_certs),
            'additional_net_target_records': sum(t['compound_id'] in new_certs for t in targets),
            'baseline_net_status_counts': dict(Counter(t['baseline_net_status'] for t in targets)),
            'sensitivity_net_status_counts': dict(Counter(t['sensitivity_net_status'] for t in targets)),
            'reactant_denominator_both_hypothetical_directions': len({m['compound_id'] for r in reactions for side in ('left', 'right') for m in r[side]})},
        'targets': targets, 'startup_scenarios': scenarios, 'admitted_reactions': admitted,
        'excluded_completions': excluded, 'additional_net_certificates': list(new_certs.values()),
        'certificate_baseline_reactions': used_reactions,
        'compounds': [compounds[c] for c in sorted(keep)],
        'external_exchange_compound_ids': sorted(exchange), 'co2_compound_id': co2,
        'solver': baseline_net['solver'], 'scipy_version': scipy.__version__,
        'baseline_certificate_policy': 'Retain exact baseline witnesses for previously feasible structures; re-solve all other distinct target structures against the augmented network.',
        'claim_boundary': BOUNDARY}


def run():
    names = ['phase1-full-balanced-network', 'phase1-candidate-scope',
             'phase1-marts-completions', 'phase1-completion-protein-evidence', 'phase1-candidate-net-flux']
    paths = [Path('data/reports', n + '.json') for n in names]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in inputs:
        for path, digest in report.get('source_sha256', {}).items():
            if path in hashes and hashes[path] != digest:
                raise ValueError('Pinned input mismatch: ' + path)
    report = build(*inputs); report['source_sha256'] = hashes
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    for folder in ('data/reports', 'docs/data'):
        Path(folder, 'phase1-completion-connectivity.json').write_text(payload)
    groups = [('target', 'targets', 'cannabisdb_id'), ('startup_scenario', 'startup_scenarios', 'id'),
        ('admitted_reaction', 'admitted_reactions', 'id'), ('excluded_completion', 'excluded_completions', 'id'),
        ('certificate', 'additional_net_certificates', 'compound_id'),
        ('baseline_reaction', 'certificate_baseline_reactions', 'id'), ('compound', 'compounds', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        rows.extend((kind, row[key], row) for row in report[collection])
    digest = hashlib.sha256(payload.encode()).hexdigest()
    count = write_rows(rows, digest, Path('data/derived/phase1-completion-connectivity.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': digest}), flush=True)


if __name__ == '__main__':
    run()
