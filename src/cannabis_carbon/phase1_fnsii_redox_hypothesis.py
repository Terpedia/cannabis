"""Exact small-molecule balance under an explicit, unverified carrier cycle."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

NARINGENIN = 'structure:7f38b2f3f512da24973fffd733f3f15112403ad6b0a3a5b99675c7ae0e3c968f'
APIGENIN = 'structure:94bde75897edb90cbe6f4eb39155b109ef67a6a3c590d8a965fa652846e12b0a'


def describe(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError('Requires an exact parseable reactive-part structure')
    elements = Counter(a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms())
    return {'smiles': smiles, 'formula': CalcMolFormula(mol),
        'elements': dict(sorted(elements.items())), 'formal_charge': Chem.GetFormalCharge(mol)}


def audit(stoichiometry, participants):
    delta = Counter()
    charge = 0
    for identity, coefficient in stoichiometry.items():
        part = participants[identity]
        for element, count in part['elements'].items():
            delta[element] += coefficient * count
        charge += coefficient * part['formal_charge']
    return {'element_delta': dict(sorted(delta.items())), 'charge_delta': charge,
        'balanced': not any(delta.values()) and charge == 0}


def build():
    paths = [Path(p) for p in (
        'data/reports/phase1-fnsii-alternative-audit.json',
        'data/reports/phase1-full-balanced-network.json',
        'data/raw/phase1-balance-reference-catalog.json',
        'data/curation/cpr-fnsii-carrier-interface-review.json',
        'data/curation/fnsii-rice-primary-assay-review.json')]
    alternative, network, catalog, interface, assay = [json.loads(p.read_text()) for p in paths]
    for parent in (alternative,):
        for path, digest in parent['source_sha256'].items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
                raise ValueError('Changed alternative source')
    fns = next(r['record'] for r in alternative['source_records'] if r['record']['rule_id'] == 'RHEA:57681')
    cpr = next(r for r in catalog if r['rule_id'] == 'RHEA:24041')
    left, right = [side.split('.') for side in fns['reaction_smarts'].split('>>')]
    cpr_left, cpr_right = [side.split('.') for side in cpr['reaction_smarts'].split('>>')]
    compounds = {c['id']: c for c in network['compounds']}
    # Preserve the exact stereochemistry/protonation already used in the gap.
    gap = alternative['parent_fnsi_gap']['reaction']
    assert NARINGENIN in {p['compound_id'] for p in gap['right']}
    assert APIGENIN in {p['compound_id'] for p in gap['left']}
    participants = {
        NARINGENIN: {**describe(compounds[NARINGENIN]['smiles']), 'source_compound': compounds[NARINGENIN]},
        APIGENIN: {**describe(compounds[APIGENIN]['smiles']), 'source_compound': compounds[APIGENIN]},
        'RHEA-COMP:11964': {**describe(left[1]), 'reactive_part': 'CHEBI:57618',
            'scope': 'protein-bound reactive part; invariant protein context assumed, not measured'},
        'RHEA-COMP:11965': {**describe(right[1]), 'reactive_part': 'CHEBI:58210',
            'scope': 'protein-bound reactive part; invariant protein context assumed, not measured'},
        'NADPH': {**describe(cpr_left[2]), 'source_component': 'RHEA:24041:left:2'},
        'NADP': {**describe(cpr_right[2]), 'source_component': 'RHEA:24041:right:2'},
        'oxygen': describe(left[2]), 'water': describe(right[2]), 'proton': describe(right[4]),
    }
    steps = [
        {'id': 'hypothesis:exact-fnsii-bound-carrier',
         'basis': 'Exact substrate/product hypothesis, not automatic generic-Rhea instantiation or experimental stereo assignment',
         'stoichiometry': {NARINGENIN: -1, 'RHEA-COMP:11964': -1, 'oxygen': -1,
             APIGENIN: 1, 'RHEA-COMP:11965': 1, 'water': 2, 'proton': 2}},
        {'id': 'hypothesis:bound-carrier-regeneration',
         'basis': 'Assumed lumped two-electron regeneration; not RHEA:24041 and not a measured elementary mechanism',
         'stoichiometry': {'NADPH': -1, 'RHEA-COMP:11965': -1, 'proton': -2,
             'NADP': 1, 'RHEA-COMP:11964': 1}},
    ]
    total = Counter()
    for step in steps:
        step['model_eligible'] = False
        step['balance'] = audit(step['stoichiometry'], participants)
        if not step['balance']['balanced']:
            raise ValueError('Unbalanced proposed step')
        total.update(step['stoichiometry'])
    net = {k: v for k, v in total.items() if v}
    balance = audit(net, participants)
    if not balance['balanced']:
        raise ValueError('Unbalanced proposed net')
    return {'schema': 'cannabis-fnsii-redox-hypothesis-v1', 'model_eligible': False,
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'source_records': {'generic_fnsii': fns, 'cpr_heme_reduction': cpr},
        'carrier_interface_review': interface, 'primary_assay_review': assay,
        'related_fnsi_gap_id': gap['id'], 'target_ids': ['CDB005071', 'CDB005072'],
        'participants': participants, 'steps': steps,
        'net': {'stoichiometry': net, 'balance': balance,
            'equation': '(2S)-naringenin + NADPH(4-) + O2 -> apigenin(1-) + NADP(3-) + 2 H2O'},
        'canceled_participants': {k: 'same identifier produced/consumed across assumed cycle' for k, v in total.items() if v == 0},
        'protonation_finding': 'Exact apigenin monoanion requires two product protons in the carrier step, versus one in the generic source; both cancel against the assumed regeneration step.',
        'unverified_requirements': ['Exact Cannabis substrate/product activity and stereoselectivity',
            'Compatible protein-bound carrier partnership and compartment',
            'Complete carrier electron/proton cycle and invariant protein context',
            'Endogenous NADPH supply/regeneration and all upstream CO2-only inputs'],
        'compatible_fnsii_partners': [], 'external_organic_seeds': [],
        'summary': {'balanced_hypothetical_steps': 2, 'new_exact_enzyme_assignments': 0,
            'candidate_model_changed': False, 'net_route_rescues_claimed': 0},
        'claim_boundary': 'Stoichiometric hypothesis only. Reactive-part accounting assumes unchanged protein context; it does not audit an entire protein or establish electron-transfer mechanism, carrier pairing, native activity, startup feasibility or atom provenance. No cancellation between distinct heme/flavin carriers, no free-flavin substitution, no network promotion. NADPH carbon remains explicit on both sides; atom tracing is deferred.'}


if __name__ == '__main__':
    report = build()
    Path('data/reports/phase1-fnsii-redox-hypothesis.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({'net': report['net'], 'summary': report['summary']}))
