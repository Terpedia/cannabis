"""Exact-participant supply diagnostics; not an executable plant pathway claim."""
import copy
import hashlib
import json
from collections import Counter
from pathlib import Path

from .phase1_marts_completions import balanced
from .phase1_net_flux import NetModel
from .phase1_scope import write_rows
from .phase1_screened_overlay import build_overlay


def build(network, expanded, references, search, sensitivity):
    compounds = {c['id']: c for c in network['compounds']}
    reactions = {r['id']: r for r in network['reactions']}
    base = set(expanded['candidate_reaction_evidence_ids'])
    if not base <= reactions.keys():
        raise ValueError('Unknown baseline equation')
    source_rows = {r['reaction_id']: r for r in references['rows']}
    if len(source_rows) != len(references['rows']) or len(search['rows']) != len(source_rows):
        raise ValueError('Duplicate or incomplete discovery rows')
    for row in search['rows']:
        if row['reaction_id'] not in source_rows or row['reference_matches'] != source_rows[row['reaction_id']]['reference_matches']:
            raise ValueError('Reference annotation join changed')
    evidence = build_overlay({'reactions': network['reactions'], 'hypotheses': []}, search,
                             'phase1-plant-purine-search.json')['enzyme_evidence']
    searches = {r['reaction_id']: r for r in search['rows']}
    for e in evidence:
        e['evidence_class'] = searches[e['reaction_id']]['evidence_class']
    added = {e['reaction_id'] for e in evidence}
    if len(added) != len(evidence) or base & added:
        raise ValueError('Duplicate or previously admitted supplement')
    for rid in base | added:
        if not balanced([reactions[rid]['left'], reactions[rid]['right']], compounds):
            raise ValueError('Unbalanced candidate equation')
    exchange = set(expanded['external_exchange_compound_ids'])
    if exchange != set(sensitivity['external_exchange_compound_ids']):
        raise ValueError('Exchange boundary changed')
    co2 = expanded['co2_compound_id']
    if {c for c in exchange if compounds[c]['carbon_count']} != {co2} or compounds[co2]['smiles'] != 'O=C=O':
        raise ValueError('Organic carbon exchange forbidden')
    forbidden = {c['id'] for c in sensitivity['constraints']}
    if len(forbidden) != len(sensitivity['constraints']):
        raise ValueError('Duplicate direction restriction')
    pathway_ids = {m['reaction_id'] for a in references['family_audit'] for m in a['equation_matches']}
    if not pathway_ids <= base | added:
        raise ValueError('Plant annotation equation lacks candidate evidence')
    participant_ids = {m['compound_id'] for rid in pathway_ids for side in ('left', 'right') for m in reactions[rid][side]}
    targets = [copy.deepcopy(t) for t in expanded['targets'] if t['new_net_certificate']]
    probe_ids = participant_ids | {t['compound_id'] for t in targets}
    scenarios = []
    for name, ids in [('restricted-candidates', base), ('restricted-plus-plant-hypotheses', base | added),
                      ('restricted-full-catalog-chemistry', set(reactions))]:
        model = NetModel([reactions[rid] for rid in sorted(ids)], exchange, forbidden_step_ids=forbidden)
        results = []
        for cid in sorted(probe_ids):
            results.append({'compound_id': cid, **model.solve(cid)})
        scenarios.append({'id': name, 'equation_count': len(ids),
            'allowed_directed_steps': len(model.steps), 'results': results,
            'participant_status_counts': dict(Counter(r['status'] for r in results if r['compound_id'] in participant_ids)),
            'target_status_counts': dict(Counter(r['status'] for r in results if r['compound_id'] in {t['compound_id'] for t in targets}))})
        print(name, scenarios[-1]['participant_status_counts'], flush=True)
    used_ids = pathway_ids | {s['reaction_id'] for scenario in scenarios for r in scenario['results'] for s in r.get('steps', [])}
    used_compounds = probe_ids | exchange | {m['compound_id'] for rid in used_ids for side in ('left', 'right') for m in reactions[rid][side]}
    aliases = {}
    for t in network['targets']:
        aliases.setdefault(t['compound_id'], []).append({'cannabisdb_id': t['cannabisdb_id'], 'label': t['label']})
    report = {'schema': 'cannabis-carbon.phase1-purine-precursor-audit.v1',
        'summary': {'plant_equations': len(pathway_ids), 'exact_pathway_participants': len(participant_ids),
                    'focused_target_records': len(targets), 'distinct_probe_structures': len(probe_ids), 'added_candidate_equations': len(added)},
        'participant_compound_ids': sorted(participant_ids), 'focused_targets': targets,
        'plant_family_audit': copy.deepcopy(references['family_audit']), 'enzyme_evidence': evidence,
        'baseline_candidate_reaction_ids': sorted(base), 'added_candidate_reaction_ids': sorted(added),
        'scenarios': scenarios, 'constraints': copy.deepcopy(sensitivity['constraints']),
        'external_exchange_compound_ids': sorted(exchange), 'co2_compound_id': co2,
        'reactions': [reactions[rid] for rid in sorted(used_ids)],
        'compounds': [{**compounds[cid], 'exact_cannabisdb_matches': aliases.get(cid, [])} for cid in sorted(used_compounds)],
        'solver': copy.deepcopy(expanded['solver']),
        'claim_boundary': 'Focused diagnostic, not a recomputation of whole-metabolome completeness. Each participant is independently tested for net production; this does not establish joint supply, pathway ordering, zero-pool startup, physiology or enzyme activity. Full-catalog chemistry includes reactions without Cannabis protein evidence. All three scenarios forbid only the same five analyst-selected reverse steps; other directions remain hypothetical. Carbon-free exchanges and regenerated pre-existing internal pools remain permissive. Exact identities and full coefficients are retained. Solver-reported infeasibility is not biological absence. Unreviewed plant homology remains distinct from established activity. Published baseline and expanded scenarios are unchanged. Atom tracing remains deferred.'}
    return add_gap_queue(report)


def add_gap_queue(report):
    """Derive a curation queue from saved witnesses, without rerunning the solver."""
    report = copy.deepcopy(report)
    scenarios = {s['id']: s for s in report['scenarios']}
    candidate = {r['compound_id']: r['status'] for r in scenarios['restricted-plus-plant-hypotheses']['results']}
    evidence = set(report['baseline_candidate_reaction_ids']) | set(report['added_candidate_reaction_ids'])
    reactions = {r['id']: r for r in report['reactions']}
    gaps = {}
    for result in scenarios['restricted-full-catalog-chemistry']['results']:
        if candidate[result['compound_id']] == 'exact-net-conversion-hypothesis':
            continue
        for step in result.get('steps', []):
            rid = step['reaction_id']
            if rid in evidence:
                continue
            gap = gaps.setdefault(rid, {'reaction_id': rid, 'source_joins': reactions[rid]['sources'],
                'left': reactions[rid]['left'], 'right': reactions[rid]['right'], 'selected_uses': [],
                'status': 'no-candidate-link-in-this-model',
                'claim_boundary': 'Selected certificate membership, not necessity or biological absence. Review existing search evidence and exact source direction before new protein searches.'})
            gap['selected_uses'].append({'probe_compound_id': result['compound_id'], **step})
    report['catalog_candidate_gaps'] = sorted(gaps.values(), key=lambda g: (-len(g['selected_uses']), g['reaction_id']))
    report['summary']['selected_catalog_candidate_gaps'] = len(gaps)
    return report


def run():
    names = ('phase1-full-balanced-network', 'phase1-expanded-candidate-net', 'phase1-plant-purine-references',
             'phase1-plant-purine-search', 'phase1-direction-sensitivity')
    paths = [Path('data/reports', n + '.json') for n in names]
    reports = [json.loads(p.read_text()) for p in paths]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for report in reports:
        for path, sha in report.get('source_sha256', {}).items():
            if path in hashes and hashes[path] != sha:
                raise ValueError('Source lineage mismatch')
    if reports[3]['source_discovery_sha256'] != hashes[str(paths[2])]:
        raise ValueError('Protein discovery snapshot mismatch')
    report = build(*reports)
    report['source_sha256'] = hashes
    export_report(report)


def export_report(report):
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-purine-precursor-audit.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('scenario', 'scenarios', 'id'), ('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id'), ('enzyme_evidence', 'enzyme_evidence', 'id'), ('candidate_gap', 'catalog_candidate_gaps', 'reaction_id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, key, id_key in groups:
        rows.extend((kind, r[id_key], r) for r in report[key])
    count = write_rows(rows, sha, Path('data/derived/phase1-purine-precursor-audit.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
