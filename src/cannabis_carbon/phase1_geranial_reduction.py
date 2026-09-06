"""Substrate-specific external-organism evidence; no Cannabis enzyme assignment."""
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from .phase1_catalog import stable_id
from .phase1_marts_completions import balanced

REFERENCE = 'balanced-equation:87fee30ea8941caec912cf4939585d22d8e11ebcb55b4cdf2a494f0a34fc64fb'
NADPH = 'structure:e25665b655013cd32bb20473bf02ffc0de259c21e04487bcb53f6dfcf320ba4b'
NADP = 'structure:dd271e39eba9fcbd95be488c513e0233a919dd8d6c9b10bd3c6704f65bb39573'
PROTON = 'structure:36d727bf22827a3587ca2cdaeac8f8150bd48f0ca0a6c24c52b7c577df4547e0'
SOURCE = 'https://pubs.acs.org/doi/10.1021/acscatal.1c05334'
BOUNDARY = ('Geranial-to-(S)-citronellal product-forming channel supported by an in vitro '
    'GluER cascade outside Cannabis. Candidate Cannabis chemistry only; no Cannabis enzyme '
    'assignment or in vivo pathway established. Finite enantiomeric excess does not mean '
    'exclusive S production; assay selectivity is not imposed as universal stoichiometry. '
    'Glucose used for assay cofactor regeneration is not an external model carbon input.')


def build(current, network):
    compounds = {c['id']: c for c in current['compounds']}
    ref = next(r for r in network['reactions'] if r['id'] == REFERENCE)
    assert {NADPH, PROTON} <= {p['compound_id'] for p in ref['left']}
    assert NADP in {p['compound_id'] for p in ref['right']}
    smiles = ['CC(C)=CCC/C(C)=C/C=O', 'CC(C)=CCC[C@H](C)CC=O']
    ids = [stable_id('structure', Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=True)) for s in smiles]
    used = [*ids, NADPH, PROTON, NADP]
    if not all(c in compounds for c in used):
        raise ValueError('Exact substrate, product or cofactor absent')
    sides = [[{'compound_id': c, 'coefficient': 1} for c in side]
             for side in ([ids[0], NADPH, PROTON], [ids[1], NADP])]
    if not balanced(sides, compounds):
        raise ValueError('Geranial reduction failed element/isotope/charge balance')
    reaction = {'id': stable_id('geranial-reduction-hypothesis', sides),
        'left': sides[0], 'right': sides[1], 'hypothesis_type': 'geranial-S-citronellal',
        'source_url': SOURCE, 'source_evidence_type': 'substrate-specific-cascade-biochemistry-outside-Cannabis',
        'direction_status': 'proposed-forward-only', 'enzyme_evidence_ids': [],
        'balance_status': 'independently-element-isotope-charge-balanced', 'claim_boundary': BOUNDARY}
    evidence = {'doi': '10.1021/acscatal.1c05334', 'pmcid': 'PMC8787751',
        'source_url': SOURCE, 'source_locations': ['Results and Discussion', 'Figure 3B', 'Abstract'],
        'enzyme_name': 'GluER', 'enzyme_organism': 'Gluconobacter oxydans',
        'substrate_supply': 'CgrAlcOx oxidation of geraniol supplies geranial; sequential cascade separates oxidation and reduction.',
        'reducing_cofactor': 'NADPH', 'assay_regeneration': 'Bacillus subtilis glucose dehydrogenase, glucose and NADP+',
        'reported_cascade_conversion_percent': 95.3, 'reported_S_enantiomeric_excess_percent': 99.2,
        'exclusive_S_product_claimed': False, 'Cannabis_activity_established': False,
        'cofactor_structure_reference_role': 'Exact charged NADPH/NADP identities only; not geranial substrate evidence.'}
    return {'schema': 'cannabis-carbon.phase1-geranial-reduction.v1', 'reactions': [reaction],
        'compounds': [compounds[c] for c in sorted(used)], 'reference_reactions': [ref],
        'biochemical_evidence': [evidence], 'claim_boundary': BOUNDARY,
        'summary': {'balanced_proposed_equations': 1, 'new_compound_structures': 0, 'coverage_gain_claimed': 0}}


def run():
    paths = [Path('data/reports/phase1-' + n + '.json') for n in ('alcohol-acetates-net', 'full-balanced-network')]
    docs = [json.loads(p.read_bytes()) for p in paths]
    for doc in docs:
        for p, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Stale source')
    report = build(*docs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-geranial-reduction.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
