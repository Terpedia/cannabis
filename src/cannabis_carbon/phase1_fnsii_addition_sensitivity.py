"""Factorial CHI/FNS-II stoichiometric sensitivity, never enzyme promotion."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem
from .phase1_chalcone_reference_search import RID as CHI
from .phase1_fnsii_redox_hypothesis import build as redox_build
from .phase1_net_flux import NetModel

FNSII = 'hypothesis:exact-naringenin-apigenin-nadph-fnsii'


def added_reaction(hypothesis, compounds):
    index = {}
    for c in compounds:
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(c['smiles']), isomericSmiles=True)
        index.setdefault(canonical, []).append(c['id'])
    resolved = {}
    for key in hypothesis['net']['stoichiometry']:
        smiles = hypothesis['participants'][key]['smiles']
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
        matches = index.get(canonical, [])
        if len(matches) != 1:
            raise ValueError('Exact unique structure resolution required: ' + key)
        resolved[key] = matches[0]
    stoich = hypothesis['net']['stoichiometry']
    return {'id': FNSII,
        'left': [{'compound_id': resolved[k], 'coefficient': -n} for k, n in stoich.items() if n < 0],
        'right': [{'compound_id': resolved[k], 'coefficient': n} for k, n in stoich.items() if n > 0],
        'balance_status': 'independently-balanced', 'model_eligible': False,
        'direction_status': 'forward-only sensitivity assumption; not established Cannabis physiology',
        'enzyme_evidence_ids': [], 'source_participant_resolution': resolved,
        'sources': [{'source_layer': 'explicit-unverified-redox-hypothesis',
            'source_path': 'data/reports/phase1-fnsii-redox-hypothesis.json',
            'source_reaction_id': FNSII,
            'evidence_boundary': 'Requires unverified exact activity and bound-carrier cycle; not a curated Rhea net reaction'}]}


def run():
    paths = [Path('data/reports/phase1-' + name + '.json') for name in (
        'remaining-candidate-net', 'full-balanced-network', 'fnsii-redox-hypothesis',
        'chalcone-addition-sensitivity')]
    candidate, network, hypothesis, chi_parent = [json.loads(p.read_text()) for p in paths]
    if hypothesis != redox_build():
        raise ValueError('Stale redox derivation')
    base = [r for r in network['reactions'] if r['id'] in candidate['candidate_reaction_evidence_ids']]
    if len(base) != len(candidate['candidate_reaction_evidence_ids']):
        raise ValueError('Incomplete baseline')
    baseline = next(s for s in candidate['scenarios'] if s['id'] == 'eight-reverse-steps-forbidden')
    chi = next(r for r in network['reactions'] if r['id'] == CHI)
    fns = added_reaction(hypothesis, network['compounds'])
    if {CHI, FNSII} & {r['id'] for r in base}:
        raise ValueError('Additions already in baseline')
    exchanges = candidate['external_exchange_compound_ids']
    compounds = {c['id']: c for c in network['compounds']}
    if [compounds[c]['smiles'] for c in exchanges if compounds[c]['carbon_count']] != ['O=C=O']:
        raise ValueError('Unexpected carbon exchange')
    targets = {t['cannabisdb_id']: t for t in network['targets']}
    scenarios = []
    for name, additions in [('baseline', []), ('CHI-only', [chi]),
            ('FNSII-only', [fns]), ('CHI-and-FNSII', [chi, fns])]:
        reactions = base + additions
        forbidden = baseline['forbidden_step_ids'] + [r['id'] + ':hypothetical-right-to-left' for r in additions]
        model = NetModel(reactions, exchanges, forbidden)
        rows = []
        for tid in hypothesis['target_ids']:
            target = targets[tid]
            result = model.solve(target['compound_id'])
            rows.append({k: target[k] for k in ('cannabisdb_id', 'compound_id', 'label')} | result)
        scenario = {'id': name, 'added_reaction_ids': [r['id'] for r in additions],
            'reaction_count': len(reactions), 'forbidden_step_ids': forbidden,
            'rows': rows, 'summary': dict(Counter(r['status'] for r in rows))}
        scenarios.append(scenario)
        print(name, json.dumps(scenario['summary']), flush=True)
    report = {'schema': 'cannabis-fnsii-addition-sensitivity-v1', 'model_eligible': False,
        'scenarios': scenarios, 'reactions': base + [chi, fns], 'compounds': network['compounds'],
        'baseline_candidate_reaction_evidence_ids': candidate['candidate_reaction_evidence_ids'],
        'external_exchange_compound_ids': exchanges,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'hypothetical_fnsii_reaction_id': FNSII, 'hypothetical_chi_reaction_id': CHI,
        'summary': {'target_records_tested': len(hypothesis['target_ids']), 'scenarios': 4,
            'new_exact_enzyme_assignments': 0, 'published_candidate_model_changed': False},
        'claim_boundary': 'Four counterfactuals for two preselected target records, not a metabolome-wide completeness claim. All baseline equations and eight direction exclusions retained. Added CHI and FNS-II equations are forward-only assumptions, not new enzyme assignments. CO2 is the only external carbon source; NADPH is internal and cannot be depleted in net. Conserved internal pools may pre-exist, so certificates do not prove startup, cofactor synthesis, carrier compatibility, compartments, physiological flux or atom provenance. Historical model and confirmed completeness unchanged.'}
    Path('data/reports/phase1-fnsii-addition-sensitivity.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    run()
