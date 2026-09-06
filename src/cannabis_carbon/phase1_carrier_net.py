"""Full-inventory, checkpointed carrier-regenerating net sensitivity."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from .phase1_net_flux import NetModel, exact_net


def build(network):
    carriers = {c['id'] for c in network['compounds'] if c.get('identity_status') == 'source-defined-macromolecular-carrier'}
    return NetModel([], network['external_exchange_compound_ids'],
                    directed_steps=network['directed_steps'], conserved_ids=carriers)


def validate(result, model):
    if result['status'] != 'exact-net-conversion-hypothesis':
        return
    lookup = {s['id']: s for s in model.steps}
    selected = [lookup[s['step_id']] for s in result['steps']]
    net = exact_net(selected, [s['extent'] for s in result['steps']])
    if any(net.get(c, 0) for c in model.conserved_ids):
        raise ValueError('Carrier pool is depleted or accumulated')
    if any(v < 0 for c, v in net.items() if c not in model.exchange_ids):
        raise ValueError('Internal depletion')
    if net.get(result['compound_id'], 0) < 1:
        raise ValueError('No target yield')
    if {c: str(-v) for c, v in net.items() if v < 0} != result['external_net_consumption']:
        raise ValueError('Input accounting mismatch')
    if {c: str(v) for c, v in net.items() if v > 0} != result['net_exports']:
        raise ValueError('Output accounting mismatch')


def run():
    path = Path('data/reports/phase1-carrier-network.json')
    network = json.loads(path.read_bytes()); sha = hashlib.sha256(path.read_bytes()).hexdigest()
    source_paths = [path, Path(__file__), Path(__file__).with_name('phase1_net_flux.py')]
    source_hashes = {str(p.resolve().relative_to(Path.cwd())): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    checkpoint = Path('data/derived/phase1-carrier-net-checkpoint.ndjson')
    model = build(network); results = {}
    if checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            record = json.loads(line)
            if record['source_sha256'] != source_hashes:
                raise ValueError('Checkpoint model mismatch')
            result = record['result']; validate(result, model)
            if result['compound_id'] in results:
                raise ValueError('Duplicate checkpoint result')
            results[result['compound_id']] = result
    ids = list(dict.fromkeys(t['compound_id'] for t in network['targets']))
    print(json.dumps({'unique_targets': len(ids), 'resumed_results': len(results),
                      'directed_steps': len(model.steps), 'conserved_carriers': len(model.conserved_ids)}), flush=True)
    with checkpoint.open('a') as output:
        for i, cid in enumerate(ids, 1):
            if cid in results:
                continue
            result = {'compound_id': cid, **model.solve(cid)}
            validate(result, model)
            result['carrier_pool_status'] = 'each-source-carrier-net-zero-required'
            result['full_balance_status'] = 'not-established; carrier and polymer source review remains'
            output.write(json.dumps({'source_sha256': source_hashes, 'result': result}, separators=(',', ':')) + '\n')
            output.flush(); results[cid] = result
            if i % 25 == 0 or result['status'] == 'solver-incomplete-or-failed':
                print(json.dumps({'completed': i, 'unique_targets': len(ids),
                                  'status_counts': dict(Counter(r['status'] for r in results.values()))}), flush=True)
    targets = [{**t, 'carrier_regenerating_sensitivity_status': results[t['compound_id']]['status']}
               for t in network['targets']]
    report = {'schema': 'cannabis-carbon.phase1-carrier-net.v1', 'targets': targets,
              'results': list(results.values()), 'conserved_carrier_ids': sorted(model.conserved_ids),
              'source_sha256': source_hashes,
              'summary': {'inventory_records': len(targets), 'unique_targets': len(results),
                          'inventory_status_counts': dict(Counter(t['carrier_regenerating_sensitivity_status'] for t in targets)),
                          'qualified_co2_pathway_coverage': None},
              'claim_boundary': 'Full candidate-network sensitivity with each distinct source carrier constrained to net zero, '
                'not confirmed pathway coverage or complete balance validation. Pre-existing regenerated carrier pools may be required; '
                'their synthesis and startup are not demonstrated. All inherited candidate directions and unresolved scaffold equations '
                'remain in this explicitly permissive scenario; source/polymer/template identity auditing remains incomplete. '
                'Solver infeasibility is model-specific, not biological absence. No historical certificates are reused.'}
    Path('data/reports/phase1-carrier-net.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
