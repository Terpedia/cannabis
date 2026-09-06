"""Structure-only inventory sweep beyond the historical TG-labelled review set."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_triglyceride_symmetry import analyze, build as build_symmetry


def build(network, acylation, symmetry, parent):
    covered = {t['cannabisdb_id'] for report in (acylation, symmetry)
               for t in report['targets'] if t['reaction_ids']}
    source = acylation['source_record']
    if source != symmetry['source_record']:
        raise ValueError('Source chemistry mismatch')
    compounds = {c['id']: c for c in parent['compounds']}
    review, scan = [], []
    for target in network['targets']:
        if target['cannabisdb_id'] in covered:
            continue
        mol = Chem.MolFromSmiles(compounds[target['compound_id']]['smiles'])
        result = analyze(mol, source)
        scan.append({'cannabisdb_id': target['cannabisdb_id'],
                     'compound_id': target['compound_id'], 'status': result['status']})
        if result['candidates']:
            review.append({k: target[k] for k in ('cannabisdb_id', 'compound_id', 'label', 'source_smiles')})
    report = build_symmetry({'source_record': source, 'review_targets': review}, parent)
    report['schema'] = 'cannabis-carbon.phase1-triglyceride-inventory-supplement.v1'
    report['inventory_scan'] = scan
    report['summary'].update({'inventory_target_records': len(network['targets']),
        'previously_proposed_target_records': len(covered),
        'remaining_target_records_scanned': len(scan),
        'scan_status_counts': dict(Counter(t['status'] for t in scan))})
    return report


def run():
    RDLogger.DisableLog('rdApp.warning')
    names = ('full-balanced-network', 'triglyceride-acylation',
             'triglyceride-symmetry', 'triglyceride-symmetry-net')
    paths = [Path('data/reports/phase1-' + n + '.json') for n in names]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-triglyceride-inventory-supplement.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
