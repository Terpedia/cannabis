"""Full balanced-catalog net feasibility and selected-certificate enzyme gaps."""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import scipy
from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel
from .phase1_net_view import build as attach_evidence
from .phase1_scope import write_rows

BOUNDARY = ('Chemistry-only diagnostic over the pinned balanced Terpedia catalog, '
    'not a Cannabis pathway-completeness claim. Reactions without candidate enzyme '
    'evidence are explicitly admitted only for this scenario. Every full equation '
    'and coproduct is retained; CO2 is the sole carbon exchange. Both directions '
    'and the unchanged carbon-free reservoir are permissive assumptions. Internal '
    'pools may pre-exist but cannot be depleted overall. Activity, specificity, '
    'physiological direction, compartment compatibility, energy, thermodynamics '
    'and pool origin remain unverified. A selected-certificate gap is not a proven '
    'indispensable reaction or minimum gene set. Atom tracing remains deferred.')


def build(network, candidate_scope, baseline_net, chemistry_scope):
    reactions = {r['id']: r for r in network['reactions']}
    compounds = {c['id']: c for c in network['compounds']}
    evidence = candidate_scope['candidate_reaction_evidence_ids']
    if not evidence.keys() <= reactions.keys():
        raise ValueError('Baseline candidate equation missing from full catalog')
    old_targets = {t['cannabisdb_id']: t for t in baseline_net['targets']}
    old_certs = {c['compound_id']: c for c in baseline_net['certificates']}
    scenario = next(s for s in chemistry_scope['scenarios'] if s['id'] == 'CO2-plus-all-carbon-free-species')
    if [t['cannabisdb_id'] for t in scenario['targets']] != [t['cannabisdb_id'] for t in network['targets']] or set(old_targets) != {t['cannabisdb_id'] for t in network['targets']}:
        raise ValueError('Full target inventory mismatch')
    exchanges = set(baseline_net['external_exchange_compound_ids'])
    if exchanges != set(scenario['seed_compound_ids']):
        raise ValueError('External exchange boundary changed')
    if [compounds[c]['smiles'] for c in exchanges if compounds[c]['carbon_count']] != ['O=C=O']:
        raise ValueError('CO2 must be the sole carbon source')
    co2 = baseline_net['co2_compound_id']
    model = NetModel(list(reactions.values()), exchanges)
    cache, certificates, targets = {}, {}, []
    for i, target in enumerate(scenario['targets'], 1):
        cid = target['compound_id']; old = old_targets[target['cannabisdb_id']]
        if old['compound_id'] != cid:
            raise ValueError('Target chemical identity differs from baseline')
        if cid not in cache:
            cache[cid] = old_certs[cid] if cid in old_certs else model.solve(cid)
        result = cache[cid]
        status = 'no-net-producing-catalog-equation' if result['status'] == 'no-net-producing-candidate-equation' else result['status']
        missing = []
        if status == 'exact-net-conversion-hypothesis':
            missing = sorted({s['reaction_id'] for s in result['steps'] if not evidence.get(s['reaction_id'])})
            carbon_in = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['external_net_consumption'].items())
            carbon_out = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['net_exports'].items())
            if carbon_in <= 0 or carbon_in != carbon_out or Fraction(result['external_net_consumption'].get(co2, '0')) != carbon_in:
                raise ValueError('Exact net CO2 carbon balance failed')
            certificates.setdefault(cid, {'compound_id': cid, **result,
                'net_carbon_in': str(carbon_in), 'net_carbon_out': str(carbon_out),
                'missing_candidate_reaction_ids': missing,
                'evidence_class': 'chemistry-only-with-enzyme-gaps' if missing else 'candidate-linked-net-hypothesis',
                'certificate_origin': 'unchanged-baseline-witness' if cid in old_certs else 'full-catalog-net-solve'})
        targets.append({k: target[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'carbon_count')} |
            {'net_status': status, 'startup_status': target['status'],
             'baseline_candidate_net_status': old['net_status'],
             'certificate_compound_id': cid if cid in certificates else None,
             'missing_candidate_reaction_ids': missing,
             'evidence_class': 'chemistry-only-with-enzyme-gaps' if missing else
                'candidate-linked-net-hypothesis' if cid in certificates else 'unresolved'} |
            {k: result[k] for k in ('solver_status', 'solver_message') if k in result})
        if i % 50 == 0:
            print(f'Full catalog: {i}/{len(scenario["targets"])} targets; {len(certificates)} exact certificates', flush=True)
    used_ids = {s['reaction_id'] for cert in certificates.values() for s in cert['steps']}
    selected = []
    for rid in sorted(used_ids):
        r = reactions[rid]
        if not balanced([r['left'], r['right']], compounds):
            raise ValueError('Selected full equation failed independent isotope/element/charge balance')
        selected.append({**r, 'enzyme_evidence_ids': evidence.get(rid, []),
            'missing_candidate_evidence': not bool(evidence.get(rid))})
    used_compounds = exchanges | {m['compound_id'] for r in selected for side in ('left', 'right') for m in r[side]}
    gap_targets, gap_structures = defaultdict(list), defaultdict(set)
    for t in targets:
        for rid in t['missing_candidate_reaction_ids']:
            gap_targets[rid].append(t['cannabisdb_id']); gap_structures[rid].add(t['compound_id'])
    gaps = [{'id': rid, 'reaction_id': rid, 'selected_certificate_target_ids': ids,
             'selected_certificate_target_count': len(ids), 'selected_certificate_structure_count': len(gap_structures[rid]),
             'status': 'missing-candidate-evidence-in-selected-chemistry-only-certificates',
             'next_test': 'Resolve exact full-reaction references and determine whether the step is enzymatic, spontaneous or only a catalog transformation. If enzymatic, search the Cannabis proteome using correctly linked reference proteins, then assay exact substrates and products. Review direction and every required cofactor; this ranking does not establish necessity.'}
            for rid, ids in gap_targets.items()]
    gaps.sort(key=lambda r: (-r['selected_certificate_structure_count'], r['id']))
    counts = dict(Counter(t['net_status'] for t in targets))
    additional = sum(t['net_status'] == 'exact-net-conversion-hypothesis' and t['baseline_candidate_net_status'] != 'exact-net-conversion-hypothesis' for t in targets)
    return {'schema': 'cannabis-carbon.phase1-catalog-net-gaps.v1',
        'summary': {'target_records': len(targets), 'catalog_equations': len(reactions),
            'baseline_candidate_equations': len(evidence), 'catalog_equations_without_candidate_evidence': len(reactions) - len(evidence),
            'target_status_counts': counts, 'exact_structure_certificates': len(certificates),
            'additional_chemistry_only_target_records': additional,
            'targets_with_missing_candidate_steps': sum(bool(t['missing_candidate_reaction_ids']) for t in targets),
            'selected_missing_candidate_equations': len(gaps),
            'reactant_denominator_both_hypothetical_directions': len({m['compound_id'] for r in reactions.values() for side in ('left', 'right') for m in r[side]})},
        'targets': targets, 'certificates': list(certificates.values()), 'reactions': selected,
        'compounds': [compounds[c] for c in sorted(used_compounds)], 'gap_priorities': gaps,
        'external_exchange_compound_ids': sorted(exchanges), 'co2_compound_id': co2,
        'solver': baseline_net['solver'], 'scipy_version': scipy.__version__,
        'view_scenario': 'full-catalog-chemistry-only',
        'view_boundary': f'Chemistry-only catalog diagnostic: {counts.get("exact-net-conversion-hypothesis", 0)} target records have exact net balances; {additional} are additional to the 101-record candidate baseline. Red edges lack candidate enzyme evidence. Not confirmed Cannabis pathways.',
        'claim_boundary': BOUNDARY,
        'comparison_boundary': 'Compared only with the pinned candidate baseline; inferred completion and proton-transfer hypotheses, including subsequent archive evidence supplements, are not admitted here. Exact existing baseline witnesses are retained. No evidence is inferred for an unsupported step.'}


def run():
    names = ['phase1-full-balanced-network', 'phase1-candidate-scope', 'phase1-candidate-net-flux',
        'phase1-all-reactants-scope', 'phase1-target-hypotheses', 'phase1-screened-enzyme-overlay', 'phase1-route-enzyme-overlay']
    paths = [Path('data/reports', n + '.json') for n in names]
    inputs = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in inputs:
        for p, digest in report.get('source_sha256', {}).items():
            if p in hashes and hashes[p] != digest:
                raise ValueError('Pinned input checksum mismatch')
    if inputs[3]['source_network_sha256'] != hashes[str(paths[0])]:
        raise ValueError('Chemistry scope network mismatch')
    report = attach_evidence(build(*inputs[:4]), inputs[4:]); report['source_sha256'] = hashes
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    output = Path('data/reports/phase1-catalog-net-gaps.json'); output.write_text(payload)
    groups = [('target', 'targets', 'cannabisdb_id'), ('certificate', 'certificates', 'compound_id'),
              ('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id'),
              ('gap_priority', 'gap_priorities', 'id'), ('enzyme_evidence', 'enzyme_evidence', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        rows.extend((kind, row[key], row) for row in report[collection])
    digest = hashlib.sha256(payload.encode()).hexdigest()
    count = write_rows(rows, digest, Path('data/derived/phase1-catalog-net-gaps.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'bytes': len(payload.encode()), 'sha256': digest}), flush=True)


if __name__ == '__main__':
    run()
