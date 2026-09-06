"""Restore carrier source IDs without inventing complete macromolecular formulas."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger
from .phase1_catalog import stable_id

RH = 'http://rdf.rhea-db.org/'


def projection(carrier, smiles_index):
    parts = []
    for part in carrier['reactive_parts']:
        chebis = [p['object']['value'] for p in part['properties'] if p['predicate'] == RH + 'chebi']
        if len(chebis) != 1:
            raise ValueError('ambiguous-reactive-part-chebi')
        chebi = chebis[0].rsplit('/', 1)[-1].replace('CHEBI_', 'CHEBI:')
        if chebi not in smiles_index:
            raise ValueError('reactive-part-smiles-unavailable')
        mol = Chem.MolFromSmiles(smiles_index[chebi])
        if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
            raise ValueError('reactive-part-not-fully-defined-small-structure')
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        parts.append({'source_reactive_part_id': part['id'], 'chebi_id': chebi,
                      'projected_compound_id': stable_id('structure', smiles), 'smiles': smiles})
    return parts


def reconstruct(reaction, join, carriers, smiles_index):
    source = join['source_record']
    left = source.get('source_left_corresponds_to')
    if left not in ('left', 'right'):
        raise ValueError('source-side-correspondence-unresolved')
    sides = {s: Counter({p['compound_id']: p['coefficient'] for p in reaction[s]}) for s in ('left', 'right')}
    restored = {s: Counter() for s in sides}; changes = []
    for occurrence in join['carrier_occurrences']:
        role = occurrence['sideRole']
        if role not in (RH + 'substrates', RH + 'products'):
            raise ValueError('unknown-source-side-role')
        side = left if role == RH + 'substrates' else ('right' if left == 'left' else 'left')
        token = occurrence['coefficientPredicate'].removeprefix(RH + 'contains')
        if not token.isdecimal() or int(token) < 1:
            raise ValueError('symbolic-carrier-coefficient-needs-symbolic-balance')
        amount = int(token); cid = occurrence['carrier']
        parts = projection(carriers[cid], smiles_index)
        for part in parts:
            pid = part['projected_compound_id']
            sides[side][pid] -= amount
            if sides[side][pid] < 0:
                raise ValueError('reactive-part-projection-does-not-match-model-side')
        restored[side][cid] += amount
        changes.append({'side': side, 'carrier_id': cid, 'coefficient': amount,
                        'source_occurrence': occurrence, 'reactive_part_projection': parts})
    result = {}
    for side in sides:
        merged = sides[side] + restored[side]
        result[side] = [{'compound_id': cid, 'coefficient': n} for cid, n in sorted(merged.items())]
    return {**result, 'restorations': changes}


def run():
    RDLogger.DisableLog('rdApp.warning')
    paths = [Path('data/reports/phase1-carrier-catalog-audit.json'), Path('data/raw/rhea-chebi-smiles-20260906.tsv')]
    audit = json.loads(paths[0].read_bytes())
    for p, sha in audit['source_sha256'].items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
            raise ValueError('Stale source')
    smiles = dict(line.split('\t', 1) for line in paths[1].read_text().splitlines() if line)
    carriers = {c['id']: c for c in audit['carriers']}
    rows = []; variants = {}
    for equation in audit['equations']:
        for join in equation['carrier_source_joins']:
            row = {'model_reaction_id': equation['id'], 'source_record': join['source_record']}
            try:
                converted = reconstruct(equation['model_reaction'], join, carriers, smiles)
                key = stable_id('carrier-restored-equation', [converted['left'], converted['right']])
                variants.setdefault(key, {'id': key, 'left': converted['left'], 'right': converted['right'],
                    'source_joins': [], 'full_macromolecular_balance_established': False,
                    'status': 'carrier-identities-restored-scaffold-balance-and-direction-review-pending'})['source_joins'].append({
                        **row, 'restorations': converted['restorations']})
                row.update({'status': 'reactive-part-substitution-reversed', 'restored_equation_id': key})
            except ValueError as error:
                row.update({'status': 'unresolved-reconstruction', 'reason': str(error)})
            rows.append(row)
    report = {'schema': 'cannabis-carbon.phase1-carrier-reconstruction.v1', 'source_rows': rows,
        'restored_equations': list(variants.values()), 'carriers': audit['carriers'],
        'summary': {'source_joins_reviewed': len(rows), 'restored_equation_variants': len(variants),
            'source_joins_restored': sum(r['status'] == 'reactive-part-substitution-reversed' for r in rows),
            'unresolved_reason_counts': dict(Counter(r['reason'] for r in rows if 'reason' in r)),
            'corrected_pathway_coverage_established': False},
        'claim_boundary': 'Intermediate identity reconstruction, not a runnable corrected network or a full balance certificate. '
            'Carrier IDs remain distinct from reactive parts and free nutrient species. Every replacement is matched to exact '
            'source-side projected structures; ambiguous, absent or symbolic mappings remain unresolved. Unknown carrier scaffold '
            'formulas are not invented. Separate conservation and scaffold/variable-stoichiometry balance, carrier-specific redox '
            'regeneration, source direction and polymer-only contexts still require review before pathway recomputation.',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    Path('data/reports/phase1-carrier-reconstruction.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
