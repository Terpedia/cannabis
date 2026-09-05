"""Explicit, review-only inorganic-participant completion hypotheses for MARTS."""
import hashlib
import json
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_balance_reference import concrete_participants
from .phase1_catalog import stable_id

BOUNDARY = ('Stoichiometric hypothesis only, not a corrected source record, enzyme assignment, '
    'mechanistic proof or CO2 pathway. Exact original organic substrate/product identities are retained. '
    'The reference product may be a different isomer: composition equality is not identity. '
    'All inferred inorganic inputs and outputs require source/publication review. Atom tracing is deferred.')


@lru_cache(maxsize=None)
def composition(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError('Concrete structure required')
    atoms = Counter((a.GetAtomicNum(), a.GetIsotope()) for a in Chem.AddHs(mol).GetAtoms())
    return tuple(sorted(atoms.items())), Chem.GetFormalCharge(mol)


def balanced(sides, compounds):
    totals, charges = [], []
    for side in sides:
        atoms, charge = Counter(), 0
        for p in side:
            n = p['coefficient']
            if not isinstance(n, int) or n <= 0:
                raise ValueError('Positive integral source stoichiometry required')
            counts, q = composition(compounds[p['compound_id']]['smiles'])
            for atom, count in counts:
                atoms[atom] += n * count
            charge += n * q
        totals.append(atoms); charges.append(charge)
    return totals[0] == totals[1] and charges[0] == charges[1]


def organic_pair(sides, compounds):
    organic = [[p for p in side if compounds[p['compound_id']]['carbon_count'] > 0] for side in sides]
    if any(len(side) != 1 or side[0]['coefficient'] != 1 for side in organic):
        return None
    return organic[0][0]['compound_id'], organic[1][0]['compound_id']


def build(audit, network):
    compounds = {c['id']: c for c in network['compounds']}
    existing = {r['id'] for r in network['reactions']}
    template_index = defaultdict(list)
    template_count = 0
    for reaction in network['reactions']:
        if not any(s['source_layer'] == 'terpedia-full-rhea-catalog' for s in reaction['sources']):
            continue
        for direction, sides in [('hypothetical-left-to-right', [reaction['left'], reaction['right']]),
                                  ('hypothetical-right-to-left', [reaction['right'], reaction['left']])]:
            pair = organic_pair(sides, compounds)
            if pair is None:
                continue
            inorganic = [[p for p in side if compounds[p['compound_id']]['carbon_count'] == 0] for side in sides]
            if not any(inorganic):
                continue
            if not balanced(sides, compounds):
                raise ValueError('Reference equation failed independent isotope/element/charge audit')
            key = (pair[0], composition(compounds[pair[1]]['smiles']))
            template_index[key].append({'reference_reaction_id': reaction['id'], 'reference_direction': direction,
                'reference_product_id': pair[1], 'inorganic_participants': inorganic})
            template_count += 1
    variants = {}
    for source in audit['source_ledger']:
        if source['balance_status'] != 'imbalanced':
            continue
        smiles = source['source_record']['reaction_smarts']
        variant = variants.setdefault(smiles, {'id': stable_id('MARTS-incomplete-variant', smiles),
            'original_reaction_smiles': smiles, 'source_record_ids': [], 'completion_ids': []})
        variant['source_record_ids'].append(source['id'])
    completions, used_references = {}, set()
    for smiles, variant in variants.items():
        sides = []
        for side in concrete_participants(smiles):
            members = []
            for p in side:
                cid = stable_id('structure', p['smiles'])
                compounds.setdefault(cid, {'id': cid, **{k: v for k, v in p.items() if k != 'coefficient'}})
                members.append({'compound_id': cid, 'coefficient': p['coefficient']})
            sides.append(members)
        pair = organic_pair(sides, compounds)
        if pair is None or any(len(side) != 1 for side in sides):
            variant['status'] = 'outside-single-organic-pair-with-no-inorganic-participants-scope'
            continue
        key = (pair[0], composition(compounds[pair[1]]['smiles']))
        for template in template_index.get(key, []):
            proposed = [[sides[i][0]] + template['inorganic_participants'][i] for i in range(2)]
            if not balanced(proposed, compounds):
                raise ValueError('Completion failed independent isotope/element/charge audit')
            canonical = [sorted(side, key=lambda p: p['compound_id']) for side in proposed]
            flipped = json.dumps(canonical[0], sort_keys=True) > json.dumps(canonical[1], sort_keys=True)
            if flipped:
                canonical.reverse()
            rid = stable_id('balanced-equation', canonical)
            hid = stable_id('stoichiometric-completion', [variant['id'], rid])
            completion = completions.setdefault(hid, {'id': hid, 'variant_id': variant['id'],
                'balanced_equation_id': rid, 'left': canonical[0], 'right': canonical[1],
                'original_organic_pair': list(pair),
                'marts_forward_direction': 'hypothetical-right-to-left' if flipped else 'hypothetical-left-to-right',
                'inferred_inorganic_participants_in_MARTS_orientation': template['inorganic_participants'],
                'reference_templates': [], 'existing_in_balanced_network': rid in existing,
                'balance_status': 'independently-element-isotope-and-charge-balanced',
                'status': 'structurally-inferred-stoichiometric-hypothesis-review-required',
                'enzyme_evidence_ids': [], 'required_input_status': 'all-input-supply-unestablished',
                'direction_status': 'unestablished-in-Cannabis; source display orientation only',
                'claim_boundary': BOUNDARY,
                'proposed_tests': ['Verify exact source product identity, including absolute/relative stereochemistry and minor-product assignments.',
                    'Confirm complete coproduct stoichiometry in the source publication or by quantitative substrate/product assays.',
                    'Resolve source-protein identity before screening Cannabis candidates; assay specificity, compartment and input supply.']})
            completion['reference_templates'].append({k: v for k, v in template.items() if k != 'inorganic_participants'} |
                {'product_match': 'exact-encoded-product' if template['reference_product_id'] == pair[1] else 'composition-only; not an identity or enzyme join'})
            if hid not in variant['completion_ids']:
                variant['completion_ids'].append(hid)
            used_references.add(template['reference_reaction_id'])
        variant['status'] = 'completion-hypotheses-found' if variant['completion_ids'] else 'no-compatible-reference-template'
    by_source = {sid: v for v in variants.values() for sid in v['source_record_ids']}
    targets = []
    for target in audit['targets']:
        ids = sorted({hid for m in target['marts_unbalanced_matches']
                      for hid in by_source[m['source_record_id']]['completion_ids']})
        targets.append({k: target[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'baseline_exact_balanced_participation')} |
                       {'completion_ids': ids, 'status': 'review-only-completion-hypotheses' if ids else 'no-completion-hypothesis-in-this-scope'})
    references = [r for r in network['reactions'] if r['id'] in used_references]
    needed = {p['compound_id'] for r in list(completions.values()) + references for side in ('left', 'right') for p in r[side]}
    return {'schema': 'cannabis-carbon.phase1-marts-completions.v1', 'claim_boundary': BOUNDARY,
        'summary': {'unbalanced_source_rows': sum(len(v['source_record_ids']) for v in variants.values()),
            'source_variants': len(variants), 'variant_status_counts': dict(Counter(v['status'] for v in variants.values())),
            'eligible_reference_orientations': template_count, 'completion_hypotheses': len(completions),
            'unique_balanced_equations': len({r['balanced_equation_id'] for r in completions.values()}),
            'additional_balanced_equations': len({r['balanced_equation_id'] for r in completions.values() if not r['existing_in_balanced_network']}),
            'targets_with_completions': sum(bool(t['completion_ids']) for t in targets),
            'targets_without_baseline_participation_with_completions': sum(bool(t['completion_ids']) and not t['baseline_exact_balanced_participation'] for t in targets)},
        'variants': list(variants.values()), 'completions': list(completions.values()),
        'targets': targets, 'reference_reactions': references, 'compounds': [compounds[cid] for cid in sorted(needed)],
        'metric_boundary': 'These inferred equations are a separate hypothesis catalog. Neither original source balance statuses nor candidate-enzyme/CO2-pathway metrics are changed.'}


def run():
    RDLogger.DisableLog('rdApp.warning'); RDLogger.DisableLog('rdApp.error')
    paths = [Path('data/reports/phase1-marts-audit.json'), Path('data/reports/phase1-full-balanced-network.json')]
    report = build(*(json.loads(p.read_text()) for p in paths))
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-marts-completions.json').write_text(payload)
    Path('docs/data/phase1-marts-completions.json').write_text(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    keys = [('variant', 'variants'), ('completion', 'completions'), ('target', 'targets'),
            ('reference_reaction', 'reference_reactions'), ('compound', 'compounds')]
    metadata = {k: v for k, v in report.items() if k not in {v for _, v in keys}}
    with Path('data/derived/phase1-marts-completions.ndjson').open('w') as handle:
        for kind, records in [('metadata', [metadata])] + [(k, report[v]) for k, v in keys]:
            for record in records:
                handle.write(json.dumps({'record_kind': kind, 'record_id': record.get('id', record.get('cannabisdb_id', 'metadata')),
                    'record_json': json.dumps(record, separators=(',', ':')), 'report_sha256': digest}) + '\n')
    print(json.dumps({'sha256': digest, 'bytes': len(payload.encode()), **report['summary']}))


if __name__ == '__main__':
    run()
