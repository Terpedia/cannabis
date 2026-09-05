"""Versioned identity alternatives, tested without rewriting historical targets."""
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

from .phase1_net_flux import NetModel
from .phase1_scope import write_rows
from .phase1_target_coverage import encoded_structure


def build(review, network, candidate):
    compounds = {c['id']: c for c in network['compounds']}
    by_structure = {encoded_structure(c['smiles'])[0]: c['id'] for c in compounds.values()}
    if len(by_structure) != len(compounds):
        raise ValueError('Duplicate exact network structures')
    evidence = candidate['candidate_reaction_evidence_ids']
    reactions = [r for r in network['reactions'] if r['id'] in evidence]
    if len(reactions) != len(evidence):
        raise ValueError('Incomplete candidate equation inventory')
    exchanges = candidate['external_exchange_compound_ids']
    if [c for c in exchanges if compounds[c]['carbon_count']] != [candidate['co2_compound_id']]:
        raise ValueError('Unexpected carbon exchange')
    if compounds[candidate['co2_compound_id']]['smiles'] != 'O=C=O':
        raise ValueError('Carbon source must be CO2')
    models = {s['id']: NetModel(reactions, exchanges, s['forbidden_step_ids']) for s in candidate['scenarios']}
    assertions = {a['id']: a for a in review['assertions']}
    decisions, branches = [], []
    for priority in review['priority_reviews']:
        accession = priority['cannabisdb_id']
        pair = [a for a in review['assertions'] if a['cannabisdb_id'] == accession]
        if len(pair) != 2:
            raise ValueError('Both source alternatives are required')
        choices = priority['exact_structure_supported_assertion_ids']
        if any(i not in {a['id'] for a in pair} for i in choices):
            raise ValueError('Foreign identity choice')
        selected = assertions[choices[0]] if len(choices) == 1 else None
        decisions.append({'cannabisdb_id': accession, 'name': priority['queried_name'],
            'status': 'provisional-named-structure-branch' if selected else 'unresolved',
            'provisional_assertion_id': selected['id'] if selected else None,
            'structure_fields': {k: selected[k] for k in ('canonical_smiles', 'computed_inchikey', 'computed_formula', 'computed_carbon_count')} if selected else None,
            'corroborating_pubchem_cids': selected['name_query_exact_structure_cids'] if selected else [],
            'retained_assertion_ids': [a['id'] for a in pair],
            'inherited_source_xrefs': False,
            'historical_target_replaced': False,
            'occurrence_status': 'not-validated-by-registry-name-match',
            'next_action': 'Verify Cannabis occurrence and primary identification evidence before promoting this provisional identity.' if selected else 'Resolve source stereochemistry against primary identification evidence; neither source is selected.'})
        for assertion in pair:
            cid = by_structure.get(assertion['canonical_smiles'])
            producers = []
            for reaction in network['reactions']:
                amount = sum((sign * Fraction(m['coefficient']) for side, sign in [('left', -1), ('right', 1)] for m in reaction[side] if m['compound_id'] == cid), Fraction())
                if amount:
                    producers.append({'reaction_id': reaction['id'], 'producing_orientation': 'left-to-right' if amount > 0 else 'right-to-left',
                                      'net_amount': str(abs(amount)), 'candidate_evidence_ids': evidence.get(reaction['id'], [])})
            for scenario, model in models.items():
                result = model.solve(cid) if cid else {'status': 'no-exact-structure-in-balanced-network'}
                if result['status'] == 'exact-net-conversion-hypothesis':
                    incoming = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['external_net_consumption'].items())
                    outgoing = sum(Fraction(n) * compounds[c]['carbon_count'] for c, n in result['net_exports'].items())
                    if incoming <= 0 or incoming != outgoing:
                        raise ValueError('Carbon exchange mismatch')
                    result = {**result, 'net_carbon_in': str(incoming), 'net_carbon_out': str(outgoing)}
                branches.append({'id': assertion['id'] + ':' + scenario, 'assertion_id': assertion['id'],
                    'cannabisdb_id': accession, 'scenario_id': scenario, 'compound_id': cid,
                    'canonical_smiles': assertion['canonical_smiles'], 'carbon_count': assertion['computed_carbon_count'],
                    'identity_status': 'provisional-named-structure-branch' if selected and selected['id'] == assertion['id'] else 'unselected-source-alternative',
                    'catalog_producers': producers, 'result': result})
    used = {s['reaction_id'] for b in branches for s in b['result'].get('steps', [])}
    return {'schema': 'cannabis-carbon.phase1-identity-branches.v1', 'decisions': decisions,
        'assertions': [a for a in review['assertions'] if a['cannabisdb_id'] in {d['cannabisdb_id'] for d in decisions}],
        'branches': branches, 'reactions': [{**r, 'enzyme_evidence_ids': evidence[r['id']]} for r in reactions if r['id'] in used],
        'external_exchange_compound_ids': exchanges,
        'scenario_constraints': {s['id']: s['forbidden_step_ids'] for s in candidate['scenarios']},
        'summary': {'priority_accessions': len(decisions), 'provisional_choices': sum(d['provisional_assertion_id'] is not None for d in decisions),
                    'branch_scenario_tests': len(branches), 'status_counts': dict(Counter(b['result']['status'] for b in branches)),
                    'historical_targets_changed': 0},
        'claim_boundary': 'Both source alternatives are tested independently. Provisional named-structure choices are not confirmed Cannabis identities or occurrence. Source labels, xrefs and formulas are never inherited wholesale. Exact net certificates allow pre-existing regenerated pools, hypothetical directions and carbon-free exchanges; they do not establish startup, compartmentation, physiological flux or characterized enzyme activity. Catalog producers are direction-unverified chemical leads, not all-input pathways. Historical counts and atom accounting are unchanged. Atom tracing remains deferred.'}


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-identity-conflict-review', 'phase1-full-balanced-network', 'phase1-replacement-candidate-net')]
    inputs = [json.loads(p.read_text()) for p in paths]
    for report in inputs:
        for source, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(source).read_bytes()).hexdigest() != sha:
                raise ValueError('Source lineage changed')
    report = build(*inputs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-identity-branches.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('decision', 'decisions', 'cannabisdb_id'), ('assertion', 'assertions', 'id'), ('branch', 'branches', 'id'), ('reaction', 'reactions', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, collection, key in groups:
        rows.extend((kind, r[key], r) for r in report[collection])
    count = write_rows(rows, sha, Path('data/derived/phase1-identity-branches.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}))


if __name__ == '__main__':
    run()
