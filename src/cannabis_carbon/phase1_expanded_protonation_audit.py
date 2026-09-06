"""Review proton-only identity bridges against the complete expanded chemistry model."""
import hashlib
import json
from pathlib import Path
from .phase1_protonation_audit import build
from .phase1_reaction_completion_net import assemble
from .phase1_marts_completions import balanced


def run():
    names = ('full-balanced-network', 'marts-completions', 'lipid-acylation-net',
             'triglyceride-net', 'triglyceride-symmetry-net',
             'phosphatidate-hydrolysis-net', 'glycerolipid-precursors-net', 'cardiolipin-net')
    paths = [Path('data/reports/phase1-' + name + '.json') for name in names]
    inputs = [json.loads(path.read_text()) for path in paths]
    for doc in inputs:
        for source, sha in doc.get('source_sha256', {}).items():
            if hashlib.sha256(Path(source).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    network, completions, *layers = inputs
    reactions, compounds, _ = assemble(network, completions)
    for layer in layers:
        for compound in layer['compounds']:
            prior = compounds.get(compound['id'])
            if prior and prior['smiles'] != compound['smiles']:
                raise ValueError('Conflicting exact identity')
            compounds.setdefault(compound['id'], compound)
        for reaction in layer['added_reactions']:
            if reaction['id'] in reactions and reactions[reaction['id']] != reaction:
                raise ValueError('Conflicting reaction')
            reactions.setdefault(reaction['id'], reaction)
    if len(reactions) != layers[-1]['summary']['balanced_equations']:
        raise ValueError('Incomplete expanded reaction inventory')
    if not all(balanced([r['left'], r['right']], compounds) for r in reactions.values()):
        raise ValueError('Unbalanced expanded reaction')
    report = build({'compounds': list(compounds.values()), 'reactions': list(reactions.values()),
                    'targets': network['targets']})
    status = {t['cannabisdb_id']: t['net_status'] for t in layers[-1]['targets']}
    for target in report['targets']:
        target['current_net_status'] = status[target['cannabisdb_id']]
    report['schema'] = 'cannabis-carbon.phase1-expanded-protonation-audit.v1'
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    report['summary']['expanded_balanced_equations'] = len(reactions)
    report['summary']['no_producer_records_with_review_bridges'] = sum(
        t['current_net_status'] == 'no-net-producing-equation' and bool(t['bridge_ids'])
        for t in report['targets'])
    Path('data/reports/phase1-expanded-protonation-audit.json').write_text(
        json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
