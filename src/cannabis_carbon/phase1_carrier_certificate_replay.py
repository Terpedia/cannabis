"""Exact fixed-flux certificate replay and necessary bounds across carrier choices."""
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from .phase1_net_flux import exact_net


def replay(cert, options, exchange):
    lower = defaultdict(Fraction); upper = defaultdict(Fraction)
    choices = []; ambiguous = False
    for old in cert['steps']:
        alternatives = options.get(old['step_id'], [])
        if not alternatives:
            return {'compound_id': cert['compound_id'], 'status': 'missing-directed-step-reconstruction',
                    'missing_step_id': old['step_id']}
        ambiguous |= len(alternatives) > 1
        amount = Fraction(old['extent'])
        nets = [exact_net([s], [amount]) for s in alternatives]
        for cid in set().union(*(n.keys() for n in nets)):
            lower[cid] += min(n.get(cid, 0) for n in nets)
            upper[cid] += max(n.get(cid, 0) for n in nets)
        choices.append({'original_step_id': old['step_id'], 'extent': old['extent'],
                        'candidate_step_ids': [s['id'] for s in alternatives]})
    forced = {c: str(-v) for c, v in sorted(upper.items()) if c not in exchange and v < 0}
    if forced:
        status = 'fixed-historical-flux-cannot-avoid-internal-depletion'
    elif ambiguous:
        status = 'carrier-choice-allocation-unresolved'
    elif lower.get(cert['compound_id'], 0) < Fraction(cert['target_amount']):
        status = 'fixed-historical-flux-target-yield-not-preserved'
    else:
        status = 'fixed-flux-net-accounting-survives-not-full-pathway-validation'
    return {'compound_id': cert['compound_id'], 'status': status,
            'has_multiple_carrier_choices': ambiguous, 'step_choices': choices,
            'unavoidable_internal_depletion': forced,
            'internal_net_bounds': {c: {'minimum': str(lower[c]), 'maximum': str(upper[c])}
                                   for c in sorted(upper) if c not in exchange and (lower[c] or upper[c])},
            'external_net_bounds': {c: {'minimum': str(lower[c]), 'maximum': str(upper[c])}
                                   for c in sorted(upper) if c in exchange and (lower[c] or upper[c])}}


def run():
    root = Path('data/reports')
    path = root / 'phase1-carrier-network.json'
    network = json.loads(path.read_bytes())
    current_path = root / 'phase1-geranial-reduction-net.json'
    current = json.loads(current_path.read_bytes())
    paths = [path, current_path] + [root / n for n in current['baseline_certificate_reports']]
    for p, sha in network['source_sha256'].items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
            raise ValueError('Stale carrier model')
    layers = [json.loads(p.read_bytes()) for p in paths[2:]] + [current]
    certs = {c['compound_id']: c for layer in layers for c in layer.get('certificates', []) + layer.get('new_certificates', [])}
    options = defaultdict(list)
    for s in network['directed_steps']:
        for parent in s.get('inherited_allowed_step_ids', [s['id']]):
            options[parent].append(s)
    rows = [replay(c, options, set(network['external_exchange_compound_ids'])) for c in certs.values()]
    by_id = {r['compound_id']: r for r in rows}
    targets = [{**t, 'historical_certificate_replay_status': by_id.get(t['compound_id'], {}).get('status', 'no-historical-certificate')}
               for t in network['targets']]
    report = {'schema': 'cannabis-carbon.phase1-carrier-certificate-replay.v1',
              'certificates': rows, 'targets': targets,
              'summary': {'historical_certificates_reviewed': len(rows), 'inventory_records': len(targets),
                          'certificate_status_counts': dict(Counter(r['status'] for r in rows)),
                          'inventory_replay_status_counts': dict(Counter(t['historical_certificate_replay_status'] for t in targets)),
                          'qualified_co2_pathway_coverage': None},
              'claim_boundary': 'Fixed historical reaction extents only; not a full-network feasibility recomputation. '
                'Per-species bounds independently minimize/maximize over source-carrier replacements, including fractional splitting. '
                'A negative internal upper bound proves every such replacement allocation depletes that species at these extents. '
                'Bounds need not be simultaneously attainable; absence of forced depletion does not establish feasible allocation. '
                'A unique surviving replay establishes net species accounting only, not complete formula balance, source identity fidelity of other equations, '
                'startup, physiological direction, or Cannabis activity. Alternate reaction routes or extents may restore supply. '
                'No missing carrier is supplied as an external medium input.',
              'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root / 'phase1-carrier-certificate-replay.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    run()
