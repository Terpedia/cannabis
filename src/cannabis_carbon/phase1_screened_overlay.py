"""Integrate verified homology evidence without rewriting source chemistry."""
import copy
import hashlib
import json
from pathlib import Path


def build_overlay(parent, search, full_search_report='phase1-new-protein-search.json'):
    reaction_ids = {r['id'] for r in parent['reactions']}
    hypotheses = {h['id']: h for h in parent['hypotheses']}
    alignments = {a['id']: a for a in search['passing_alignments']}
    rows = []
    for row in search['rows']:
        if row['reaction_id'] not in reaction_ids:
            raise ValueError('Search reaction missing from parent chemistry')
        if any(hypotheses[hid]['reaction_id'] != row['reaction_id'] for hid in row['hypothesis_ids']):
            raise ValueError('Search hypothesis/reaction mismatch')
        hits = [alignments[aid] for aid in row['passing_alignment_ids']]
        if {a['cannabis_accession'] for a in hits} != set(row['screened_cannabis_proteins']):
            raise ValueError('Candidate proteins do not match alignments')
        allowed = {r['accession'] for r in row['reference_matches']}
        for hit in hits:
            if hit['reference_accession'] not in allowed or not hit['passes_screen'] or hit['identity_percent'] < 30 or min(hit['query_coverage_percent'], hit['reference_coverage_percent']) < 50 or not 0 <= hit['evalue'] <= 1e-5:
                raise ValueError('Invalid screened alignment or reference join')
        if not hits:
            continue
        representatives = {}
        for hit in sorted(hits, key=lambda a: (-a['bitscore'], -a['identity_percent'], a['id'])):
            representatives.setdefault(hit['cannabis_accession'], hit)
        rows.append({'id': 'proteome-screen:' + row['reaction_id'], 'reaction_id': row['reaction_id'],
            'evidence_class': 'direction-unresolved-reference-homology-candidate',
            'screened_proteins': [{'accession': acc, 'representative_alignment': hit} for acc, hit in sorted(representatives.items())],
            'passing_alignment_ids': row['passing_alignment_ids'], 'reference_matches': row['reference_matches'],
            'screen': search['screen'], 'validation_blockers': row['validation_blockers'],
            'proposed_test': row['proposed_test'],
            'full_search_report': full_search_report,
            'alignment_selection': 'Highest bitscore alignment per protein is displayed as a representative, not a functional rank; every passing alignment ID is retained.',
            'claim_boundary': 'Exact balanced-equation identity joins homology evidence to all its target projections. Protein activity, substrate specificity, direction, and CO2 pathways remain unconfirmed.'})
    return {'schema': 'cannabis-carbon.phase1-screened-enzyme-overlay.v1', 'enzyme_evidence': rows,
        'summary': {'equations_with_screened_evidence': len(rows),
            'distinct_cannabis_proteins': len({p['accession'] for r in rows for p in r['screened_proteins']})},
        'claim_boundary': 'Additional candidate evidence only; no reaction, stoichiometry, direction, target identity, or pathway execution claim is changed.'}


def apply_overlay(parent, overlay):
    report = copy.deepcopy(parent)
    index = {e['reaction_id']: e for e in overlay['enzyme_evidence']}
    if len(index) != len(overlay['enzyme_evidence']):
        raise ValueError('Duplicate overlay reaction')
    if not index.keys() <= {r['id'] for r in report['reactions']}:
        raise ValueError('Overlay reaction missing from parent chemistry')
    existing = {e['id'] for e in report['enzyme_evidence']}
    if existing & {e['id'] for e in index.values()}:
        raise ValueError('Overlay evidence ID collision')
    report['enzyme_evidence'].extend(overlay['enzyme_evidence'])
    baseline = {h['cannabisdb_id'] for h in parent['hypotheses'] if h['has_candidate_enzyme_evidence']}
    added = 0
    for h in report['hypotheses']:
        evidence = index.get(h['reaction_id'])
        if evidence:
            added += not h['has_candidate_enzyme_evidence']
            h['evidence_ids'].append(evidence['id'])
            h['has_candidate_enzyme_evidence'] = True
            h['blockers'] = [b for b in h['blockers'] if b != 'no-candidate-enzyme-evidence-attached']
            h['candidate_evidence_boundary'] = evidence['claim_boundary']
    carbon = {t['cannabisdb_id'] for t in report['targets'] if t['carbon_count']}
    now = {h['cannabisdb_id'] for h in report['hypotheses'] if h['has_candidate_enzyme_evidence']}
    report['summary'].update(
        hypotheses_with_candidate_enzyme_evidence=sum(h['has_candidate_enzyme_evidence'] for h in report['hypotheses']),
        hypotheses_without_candidate_enzyme_evidence=sum(not h['has_candidate_enzyme_evidence'] for h in report['hypotheses']),
        carbon_bearing_targets_with_candidate_enzyme_evidence=len(now & carbon),
        baseline_carbon_bearing_targets_with_candidate_enzyme_evidence=len(baseline & carbon),
        screened_overlay_added_hypotheses=added)
    return report


def run():
    parent_path = Path('data/reports/phase1-target-hypotheses.json')
    search_path = Path('data/reports/phase1-new-protein-search.json')
    search = json.loads(search_path.read_text())
    discovery_path = Path(search['source_discovery'])
    if hashlib.sha256(discovery_path.read_bytes()).hexdigest() != search['source_discovery_sha256']:
        raise ValueError('Discovery source checksum mismatch')
    discovery = json.loads(discovery_path.read_text())
    if hashlib.sha256(parent_path.read_bytes()).hexdigest() != discovery['source_sha256'][str(parent_path)]:
        raise ValueError('Parent hypothesis snapshot mismatch')
    parent = json.loads(parent_path.read_text())
    overlay = build_overlay(parent, search)
    overlay['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [parent_path, search_path, discovery_path]}
    overlay['integrated_summary'] = apply_overlay(parent, overlay)['summary']
    Path('data/reports/phase1-screened-enzyme-overlay.json').write_text(json.dumps(overlay, separators=(',', ':')) + '\n')
    print(json.dumps(overlay['integrated_summary']))


def export_table(report_path, output):
    raw = report_path.read_bytes()
    report = json.loads(raw)
    digest = hashlib.sha256(raw).hexdigest()
    records = [{'record_kind': 'enzyme_evidence', 'record_id': e['id'],
        'record_json': json.dumps(e, separators=(',', ':')), 'report_sha256': digest} for e in report['enzyme_evidence']]
    records.append({'record_kind': 'metadata', 'record_id': 'overlay',
        'record_json': json.dumps({k: v for k, v in report.items() if k != 'enzyme_evidence'}, separators=(',', ':')), 'report_sha256': digest})
    output.write_text(''.join(json.dumps(r, separators=(',', ':')) + '\n' for r in records))
    return len(records)


if __name__ == '__main__':
    run()
