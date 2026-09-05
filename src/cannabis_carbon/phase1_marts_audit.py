"""Whole-MARTS source audit and exact-identity balanced reaction additions."""
import argparse
import hashlib
import json
import os
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from rdkit import RDLogger, rdBase
from .balance import _reaction_smiles_balance
from .phase1_balance_reference import concrete_participants
from .phase1_catalog import stable_id

TABLE = 'terpedia-489015.terpedia_core.terpene_reaction_smarts_catalog_normalized_current_v2'
QUERY = f"SELECT * FROM `{TABLE}` WHERE STARTS_WITH(rule_id, 'MARTS:') ORDER BY rule_id, reaction_smarts"
BOUNDARY = ('Whole-MARTS catalog audit, not confirmed Cannabis pathways. All required participants and '
    'exact encoded identities are retained. Source protein/EC references are leads, not Cannabis enzyme '
    'assignments. Direction and all-input supply remain unresolved. Atom tracing is deferred.')


def build(network, rows):
    compounds, reactions, ledger = {}, {}, []
    target_ids = {t['compound_id'] for t in network['targets']}
    unbalanced_participation = defaultdict(list)
    existing = {r['id'] for r in network['reactions']}
    baseline_participants = {m['compound_id'] for r in network['reactions'] for side in ('left', 'right') for m in r[side]}
    source_ids = set()
    for source in rows:
        sid = stable_id('marts-source-record', source)
        if sid in source_ids:
            raise ValueError('Duplicate source record')
        source_ids.add(sid)
        smiles = source.get('reaction_smarts')
        element, charge = _reaction_smiles_balance(smiles) if smiles and smiles.count('>>') == 1 else (None, None)
        status = 'balanced' if element and charge and element['status'] == charge['status'] == 'balanced' else 'imbalanced' if element and charge else 'not-auditable'
        record = {'id': sid, 'source_record': source, 'balance_status': status,
                  'element_balance': element, 'charge_balance': charge, 'reaction_id': None}
        if status == 'imbalanced':
            counts = [{stable_id('structure', p['smiles']): p['coefficient']
                       for p in side} for side in concrete_participants(smiles)]
            for cid in (counts[0].keys() | counts[1].keys()) & target_ids:
                unbalanced_participation[cid].append({'source_record_id': sid,
                    'left_coefficient': counts[0].get(cid, 0), 'right_coefficient': counts[1].get(cid, 0),
                    'status': 'excluded-from-balanced-network; full-stoichiometry-review-required'})
        if status == 'balanced':
            sides = []
            for participants in concrete_participants(smiles):
                members = []
                for participant in participants:
                    cid = stable_id('structure', participant['smiles'])
                    compounds[cid] = {'id': cid, **{k: v for k, v in participant.items() if k != 'coefficient'}}
                    members.append({'compound_id': cid, 'coefficient': participant['coefficient']})
                sides.append(sorted(members, key=lambda m: m['compound_id']))
            flipped = json.dumps(sides[0], sort_keys=True) > json.dumps(sides[1], sort_keys=True)
            if flipped:
                sides.reverse()
            rid = stable_id('balanced-equation', sides)
            reaction = reactions.setdefault(rid, {'id': rid, 'left': sides[0], 'right': sides[1],
                'source_links': [], 'existing_in_full_balanced_network': rid in existing,
                'balance_status': 'independently-balanced',
                'direction_status': 'unresolved-in-Cannabis; canonical side ordering is not physiology',
                'required_input_status': 'not-assessed; every full-side input is required',
                'enzyme_evidence_ids': [], 'claim_boundary': BOUNDARY})
            reaction['source_links'].append({'source_record_id': sid,
                'source_left_corresponds_to': 'right' if flipped else 'left'})
            record['reaction_id'] = rid
        ledger.append(record)
    participation = defaultdict(list)
    for reaction in reactions.values():
        left = {m['compound_id']: m['coefficient'] for m in reaction['left']}
        right = {m['compound_id']: m['coefficient'] for m in reaction['right']}
        for cid in sorted(left.keys() | right.keys()):
            participation[cid].append({'reaction_id': reaction['id'], 'left_coefficient': left.get(cid, 0),
                'right_coefficient': right.get(cid, 0), 'has_net_production_in_hypothetical_direction': left.get(cid, 0) != right.get(cid, 0)})
    targets = [{k: t[k] for k in ('cannabisdb_id', 'label', 'compound_id', 'source_url')} |
        {'baseline_exact_balanced_participation': t['compound_id'] in baseline_participants,
         'marts_balanced_matches': participation.get(t['compound_id'], []),
         'marts_unbalanced_matches': unbalanced_participation.get(t['compound_id'], []),
         'status': 'balanced-MARTS-participant' if participation.get(t['compound_id']) else
                   'unbalanced-MARTS-participant' if unbalanced_participation.get(t['compound_id']) else 'no-exact-auditable-MARTS-match'}
        for t in network['targets']]
    added = [r['id'] for r in reactions.values() if not r['existing_in_full_balanced_network']]
    new_targets = [t['cannabisdb_id'] for t in targets if t['marts_balanced_matches'] and not t['baseline_exact_balanced_participation']]
    return {'schema': 'cannabis-carbon.phase1-marts-audit.v1', 'claim_boundary': BOUNDARY,
        'rdkit_version': rdBase.rdkitVersion, 'summary': {'source_rows': len(rows),
            'source_balance_status_counts': dict(Counter(r['balance_status'] for r in ledger)),
            'balanced_equations': len(reactions), 'additional_balanced_equations': len(added),
            'balanced_participant_structures': len(compounds), 'target_records': len(targets),
            'targets_with_balanced_MARTS_matches': sum(bool(t['marts_balanced_matches']) for t in targets),
            'targets_with_unbalanced_MARTS_matches': sum(bool(t['marts_unbalanced_matches']) for t in targets),
            'targets_without_baseline_participation_with_unbalanced_MARTS_matches': sum(bool(t['marts_unbalanced_matches']) and not t['baseline_exact_balanced_participation'] for t in targets),
            'unbalanced_source_equation_strings': len({r['source_record']['reaction_smarts'] for r in ledger if r['balance_status'] == 'imbalanced'}),
            'new_balanced_target_participants': len(new_targets)},
        'new_balanced_target_ids': new_targets, 'additional_reaction_ids': added,
        'targets': targets, 'source_ledger': ledger, 'reactions': list(reactions.values()),
        'compounds': list(compounds.values()),
        'next_tests': ['Resolve full stoichiometry for imbalanced or unauditable source records.',
            'Validate source protein sequences and source reaction annotations before whole-Cannabis-proteome screening.',
            'Assess every input, direction and compartment before proposing a CO2 route.']}


def run(fetch=False):
    RDLogger.DisableLog('rdApp.warning'); RDLogger.DisableLog('rdApp.error')
    raw_path = Path('data/raw/phase1-full-marts-snapshot.json')
    if fetch:
        if raw_path.exists():
            raise ValueError('Snapshot already exists; refusing to overwrite')
        data = subprocess.check_output([os.environ.get('CANNABIS_BQ', 'bq'), '--location=us-central1',
            '--format=json', 'query', '--use_legacy_sql=false', '--maximum_bytes_billed=104857600',
            '--max_rows=100000', QUERY])
        snapshot = {'table': TABLE, 'query': QUERY, 'rows': json.loads(data)}
        raw_path.write_text(json.dumps(snapshot, separators=(',', ':')) + '\n')
    snapshot = json.loads(raw_path.read_text())
    path = Path('data/reports/phase1-full-balanced-network.json')
    report = build(json.loads(path.read_text()), snapshot['rows'])
    report['catalog_provenance'] = {k: snapshot[k] for k in ('table', 'query')}
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (path, raw_path)}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    out = Path('data/reports/phase1-marts-audit.json'); out.write_text(payload)
    Path('docs/data/phase1-marts-audit.json').write_text(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('target', report['targets']), ('source', report['source_ledger']),
              ('reaction', report['reactions']), ('compound', report['compounds'])]
    metadata = {k: v for k, v in report.items() if k not in ('targets', 'source_ledger', 'reactions', 'compounds')}
    with Path('data/derived/phase1-marts-audit.ndjson').open('w') as handle:
        for kind, records in [('metadata', [metadata])] + groups:
            for record in records:
                handle.write(json.dumps({'record_kind': kind, 'record_id': record.get('id', record.get('cannabisdb_id', 'metadata')),
                    'record_json': json.dumps(record, separators=(',', ':')), 'report_sha256': digest}) + '\n')
    print(json.dumps({'sha256': digest, 'bytes': len(payload.encode()), **report['summary']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--fetch', action='store_true')
    run(parser.parse_args().fetch)
