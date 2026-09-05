"""Audit generic FNS-II chemistry separately from exact FNS-I hypotheses."""
import csv
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from .phase1_ureidoglycolate_broad_search import lookup

RAW = Path('data/raw/fnsii-alternative')
QUERY = 'ec:1.14.19.76 AND reviewed:true AND fragment:false'
IDS = {'RHEA:57681', 'RHEA:57682'}


def inspect_components(smiles):
    sides = smiles.split('>>')
    if len(sides) != 2:
        raise ValueError('Expected two reaction sides')
    result = {}
    for side, value in zip(('left', 'right'), sides):
        parts = []
        for index, part in enumerate(value.split('.')):
            mol = Chem.MolFromSmiles(part)
            if mol is None:
                raise ValueError('Unparseable component')
            parts.append({'component_index': index, 'source_smiles': part,
                'dummy_atoms': sum(a.GetAtomicNum() == 0 for a in mol.GetAtoms()),
                'explicit_carbon_atoms': sum(a.GetAtomicNum() == 6 for a in mol.GetAtoms()),
                'formal_charge': Chem.GetFormalCharge(mol)})
        result[side] = parts
    return result


def run():
    RAW.mkdir(parents=True, exist_ok=True)
    item = lookup('reviewed-fnsii', QUERY, RAW)
    catalog_path = Path('data/raw/phase1-balance-reference-catalog.json')
    network_path = Path('data/reports/phase1-full-balanced-network.json')
    gap_path = Path('data/reports/phase1-chalcone-remaining-gaps.json')
    review_path = Path('data/curation/fnsii-carrier-review.json')
    network = json.loads(network_path.read_text())
    catalog = json.loads(catalog_path.read_text())
    selected = [r for r in catalog if r['rule_id'] in IDS]
    if {r['rule_id'] for r in selected} != IDS or len(selected) != 2:
        raise ValueError('Unexpected FNS-II source inventory')
    exclusions = [r for r in network['excluded_rhea_source_records'] if r['source_reaction_id'] in IDS]
    included = [r['id'] for r in network['reactions'] if any(s['source_reaction_id'] in IDS for s in r['sources'])]
    if included or {r['source_reaction_id'] for r in exclusions} != IDS:
        raise ValueError('FNS-II catalog inclusion state changed; review required')
    gap = json.loads(gap_path.read_text())['candidate_gaps'][0]
    records = list(csv.DictReader(Path(item['snapshot']).read_text().splitlines(), delimiter='\t'))
    paths = [catalog_path, network_path, gap_path, review_path, Path(item['snapshot']), RAW / 'reviewed-fnsii-lookup.json']
    report = {'schema': 'cannabis-fnsii-alternative-audit-v1', 'model_eligible': False,
        'parent_fnsi_gap': gap, 'alternative_ec': '1.14.19.76',
        'source_records': [{'record': r, 'component_audit': inspect_components(r['reaction_smarts'])} for r in selected],
        'network_exclusions': exclusions, 'review': json.loads(review_path.read_text()),
        'lookup': item, 'reference_leads': [{'accession': r['Entry'], 'record': r,
            'model_eligible': False, 'match_type': 'generic-FNS-II-EC-class; exact-substrate-review-required'} for r in records],
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'generic_directional_records': len(selected), 'reviewed_reference_leads': len(records),
            'new_exact_reactions': 0, 'new_exact_enzyme_assignments': 0},
        'next_steps': ['Retrieve full reference annotations and screen the full Cannabis proteome for FNS-II leads.',
            'Resolve exact naringenin/apigenin specificity and compatible reductase partners before a separate balanced scenario.',
            'Do not seed organic cofactors or substitute free flavin for a protein-bound carrier.'],
        'claim_boundary': 'Alternative enzyme-class hypothesis, not an exact replacement edge. Generic substituents and protein-carrier context prevent promoting these source SMILES into the exact small-molecule network. Prior FNS-I certificates, exclusions and atom accounting remain unchanged.'}
    Path('data/reports/phase1-fnsii-alternative-audit.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
