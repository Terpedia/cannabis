"""Explicit balanced protonation hypotheses, gated by exact source endpoint joins."""
import hashlib
import json
from pathlib import Path
from .phase1_catalog import stable_id
from .phase1_protonation_audit import fingerprint, bridge
from .phase1_marts_completions import balanced

BOUNDARY = ('Source-mapped acid-base sensitivity hypothesis, not a curated biochemical reaction. '
    'Rhea mapping origin is preserved; its pH 7.3 convention does not establish Cannabis tissue pH. '
    'Reversible proton transfer is an explicit model assumption with pKa, concentrations, '
    'compartments and physiological direction unresolved. Exact identities remain distinct. '
    'No enzyme or CO2 route is established by the equation alone.')


def build(audit, joined, parent):
    compounds = {c['id']: c for c in parent['compounds']}
    protons = [c['id'] for c in compounds.values() if c['smiles'] == '[H+]']
    if len(protons) != 1:
        raise ValueError('Expected exactly one explicit proton identity')
    proton = protons[0]
    bridges = {b['id']: b for b in audit['bridges']}
    reactions = []
    for row in joined['rows']:
        if row['status'] != 'exact-endpoints-with-source-mapping':
            continue
        if not row['mapping_evidence']:
            raise ValueError('Missing source mapping evidence')
        b = bridges[row['id']]
        a, z = b['target_compound_id'], b['reaction_participant_compound_id']
        if (a, z) != (row['target_compound_id'], row['reaction_participant_compound_id']):
            raise ValueError('Endpoint mismatch')
        left_fp, _ = fingerprint(compounds[a]['smiles'])
        right_fp, _ = fingerprint(compounds[z]['smiles'])
        replay = bridge(left_fp, right_fp)
        if replay is None or any(b[k] != value for k, value in replay.items()):
            raise ValueError('Exact proton-only bridge replay failed')
        delta = replay['protons_consumed']
        left = [{'compound_id': a, 'coefficient': 1}]
        right = [{'compound_id': z, 'coefficient': 1}]
        (left if delta > 0 else right).append({'compound_id': proton, 'coefficient': abs(delta)})
        if not balanced([left, right], compounds):
            raise ValueError('Protonation equation is not balanced')
        reactions.append({'id': stable_id('source-mapped-protonation', [a, z]),
            'left': left, 'right': right, 'hypothesis_type': 'source-mapped-protonation',
            'direction_status': 'explicit-reversible-acid-base-sensitivity-assumption',
            'balance_status': 'independently-element-isotope-charge-balanced',
            'source_url': 'https://www.rhea-db.org/help/searching-rhea',
            'source_mapping_evidence': row['mapping_evidence'], 'source_bridge_id': row['id'],
            'cannabisdb_ids': row['cannabisdb_ids'], 'enzyme_evidence_ids': [], 'claim_boundary': BOUNDARY})
    used = {m['compound_id'] for r in reactions for side in ('left', 'right') for m in r[side]}
    return {'schema': 'cannabis-carbon.phase1-source-mapped-protonation.v1',
            'reactions': reactions, 'compounds': [compounds[c] for c in sorted(used)],
            'summary': {'balanced_hypotheses': len(reactions),
                        'distinct_target_records': len({t for r in reactions for t in r['cannabisdb_ids']})},
            'claim_boundary': BOUNDARY}


def run():
    paths = [Path('data/reports/phase1-' + n + '.json') for n in
             ('expanded-protonation-audit', 'protonation-verified-join', 'cardiolipin-net')]
    docs = [json.loads(p.read_text()) for p in paths]
    for doc in docs:
        for source, sha in doc['source_sha256'].items():
            if hashlib.sha256(Path(source).read_bytes()).hexdigest() != sha:
                raise ValueError('Changed source snapshot')
    report = build(*docs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    Path('data/reports/phase1-source-mapped-protonation.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
