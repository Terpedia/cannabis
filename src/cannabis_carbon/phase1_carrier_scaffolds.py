"""Audit symbolic carrier conservation; never infer complete protein formulas."""
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

RH = 'http://rdf.rhea-db.org/'


def scaffold_candidate(carrier):
    names = [p['object']['value'] for p in carrier['properties'] if p['predicate'] == RH + 'name']
    if len(names) != 1:
        return None
    # Exact source context retained: c, cL, c2, generic cytochrome, and distinct
    # Fe-S stoichiometries never share a conservation token.
    match = re.fullmatch(r'(oxidized |reduced |Fe\(II\)-|Fe\(III\)-)(.+)', names[0])
    if not match or '[' not in match[2]:
        return None
    return {'carrier_id': carrier['id'], 'source_name': names[0],
            'symbolic_scaffold': match[2], 'source_redox_prefix': match[1],
            'status': 'source-name-derived-conservation-hypothesis-not-identity-equivalence',
            'source_properties': carrier['properties']}


def audit_equation(equation, candidates):
    residual = Counter(); unknown = set(); used = set()
    for side, sign in (('left', -1), ('right', 1)):
        for participant in equation[side]:
            cid = participant['compound_id']
            if not cid.startswith(RH + 'Compound_'):
                continue
            used.add(cid)
            if cid not in candidates:
                unknown.add(cid)
            else:
                residual[candidates[cid]['symbolic_scaffold']] += sign * participant['coefficient']
    nonzero = {k: v for k, v in sorted(residual.items()) if v}
    status = ('unresolved-carrier-scaffolds' if unknown else
              'candidate-scaffold-nonconservation' if nonzero else
              'conditional-redox-scaffold-conservation')
    return {'id': equation['id'], 'status': status, 'carrier_ids': sorted(used),
            'unresolved_carrier_ids': sorted(unknown), 'symbolic_product_minus_substrate': nonzero,
            'source_reaction_ids': sorted({j['source_record']['source_reaction_id'] for j in equation['source_joins']}),
            'full_macromolecular_balance_established': False}


def run():
    path = Path('data/reports/phase1-carrier-reconstruction.json')
    source = json.loads(path.read_bytes())
    for p, sha in source['source_sha256'].items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
            raise ValueError('Stale reconstruction inputs')
    candidates = {c['id']: item for c in source['carriers'] if (item := scaffold_candidate(c))}
    equations = [audit_equation(e, candidates) for e in source['restored_equations']]
    groups = defaultdict(list)
    for item in candidates.values():
        groups[item['symbolic_scaffold']].append(item['carrier_id'])
    report = {'schema': 'cannabis-carbon.phase1-carrier-scaffolds.v1',
              'source_sha256': {str(path): hashlib.sha256(path.read_bytes()).hexdigest()},
              'carrier_scaffold_candidates': list(candidates.values()),
              'scaffold_groups': [{'symbolic_scaffold': k, 'distinct_carrier_ids': sorted(v)} for k, v in sorted(groups.items())],
              'equations': equations,
              'summary': {'restored_equations_reviewed': len(equations),
                          'equation_status_counts': dict(Counter(e['status'] for e in equations)),
                          'source_carrier_scaffold_candidates': len(candidates),
                          'symbolic_scaffold_groups': len(groups),
                          'corrected_pathway_coverage_established': False},
              'claim_boundary': 'Symbolic conservation hypotheses derived from exact source redox descriptors, not chemical identity merges or full formula validation. '
                'Full carrier species remain distinct, with no external carrier uptake authorized. Equal scaffold counts are necessary under the proposed pairing, '
                'not sufficient evidence of complete balance, physiological direction, regeneration, startup availability, or Cannabis activity. '
                'Other carrier modifications and polymer-only participants remain unresolved; no equation is promoted into a validated network by this audit.'}
    Path('data/reports/phase1-carrier-scaffolds.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps(report['summary']))


if __name__ == '__main__':
    run()
