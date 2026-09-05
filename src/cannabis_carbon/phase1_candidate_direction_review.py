"""Explicit source-orientation audit and tests for newly enabled net routes."""
import hashlib
import json
from pathlib import Path
from .phase1_gap_annotations import fetch, assemble as annotations_from_snapshots, RH
from .phase1_reference_discovery import direction_families
from .phase1_scope import write_rows


def build(report, families, annotations):
    annotations = {a['id']: a for a in annotations}
    certificates = {c['compound_id']: c for c in report['certificates']}
    targets = [t for t in report['targets'] if t['new_net_certificate']]
    reviews = []
    for r in report['reactions']:
        if not r['new_catalog_candidate']:
            continue
        masters = {families[s['source_reaction_id']]['RHEA_ID_MASTER'] for s in r['sources']}
        if len(masters) != 1:
            raise ValueError('Ambiguous source master for direction review')
        master = next(iter(masters)); annotation = annotations[master]
        family = families[r['sources'][0]['source_reaction_id']]
        lr = [s for s in r['sources'] if s['source_reaction_id'] == family['RHEA_ID_LR']]
        if len(lr) != 1 or lr[0]['source_left_corresponds_to'] not in ('left', 'right'):
            raise ValueError('Missing unique published left-to-right source mapping')
        source_left = lr[0]['source_left_corresponds_to']
        for s in r['sources']:
            if s['source_reaction_id'] == family['RHEA_ID_RL'] and s['source_left_corresponds_to'] == source_left:
                raise ValueError('Conflicting directional source mappings')
        uses = []
        for target in targets:
            for step in certificates[target['compound_id']]['steps']:
                if step['reaction_id'] != r['id']:
                    continue
                mode = step['direction_mode']
                if mode not in ('hypothetical-left-to-right', 'hypothetical-right-to-left'):
                    raise ValueError('Unknown certificate direction')
                inputs = 'left' if mode == 'hypothetical-left-to-right' else 'right'
                uses.append({'cannabisdb_id': target['cannabisdb_id'], 'label': target['label'],
                    'step_id': step['step_id'], 'direction_mode': mode, 'extent': step['extent'],
                    'source_orientation': 'same-as-source-written' if inputs == source_left else 'opposite-to-source-written',
                    'required_inputs': r[inputs], 'outputs': r['right' if inputs == 'left' else 'left']})
        def values(predicate):
            return sorted({t['object']['value'] for t in annotation['triples'] if t['predicate']['value'] == RH + predicate})
        reviews.append({'id': 'direction-review:' + r['id'], 'reaction_id': r['id'],
            'source_master_id': master, 'source_written_equations': values('equation'),
            'citation_urls': values('citation'), 'published_family': family,
            'source_left_corresponds_to': source_left, 'source_joins': r['sources'],
            'left': r['left'], 'right': r['right'], 'enzyme_evidence_ids': r['enzyme_evidence_ids'],
            'uses': uses, 'status': 'direction-and-specificity-validation-required',
            'warning': ('The selected new routes use this reaction opposite to its source-written equation. Source ordering is not a physiological direction constraint; the reverse catalytic capability remains unestablished.' if any(u['source_orientation'] == 'opposite-to-source-written' for u in uses) else 'Selected uses follow source-written ordering, which does not establish Cannabis activity or physiological direction.'),
            'hypothesis': 'At least one linked Cannabis protein catalyzes the exact full conversion in the selected direction, with the stated cofactors and all coproducts.',
            'discriminating_tests': [
                'Assay each linked protein with the complete selected input set and verify every predicted product against chemical standards; include no-protein and inactive-protein controls.',
                'Run the source-written direction as a separate positive-control assay. Activity only in that direction does not validate the reverse route.',
                'Compare exact cofactors and close substrate analogues rather than treating a family hit as proof of specificity; identify required partners and compartments.',
                'If the selected direction is not supported, search for alternative plant-supported reactions and re-solve; do not interpret failure as absence of the target metabolite.'],
            'source_annotation': annotation})
    target_rows = [{'cannabisdb_id': t['cannabisdb_id'], 'compound_id': t['compound_id'], 'label': t['label'],
        'direction_review_ids': [r['id'] for r in reviews if any(u['cannabisdb_id'] == t['cannabisdb_id'] and u['source_orientation'] == 'opposite-to-source-written' for u in r['uses'])]} for t in targets]
    return {'schema': 'cannabis-carbon.phase1-candidate-direction-review.v1', 'reviews': reviews, 'targets': target_rows,
        'summary': {'reviewed_new_equations': len(reviews), 'new_target_records': len(targets),
            'targets_using_opposite_source_orientation': sum(bool(t['direction_review_ids']) for t in target_rows)},
        'literature_review': {'source_url': 'https://pubmed.ncbi.nlm.nih.gov/19260710/',
            'source_master_id': 'RHEA:27329', 'evidence_scope': 'primary-paper abstract and figure captions',
            'finding': 'The paper characterizes Klebsiella pneumoniae HpxO as an FAD-dependent urate oxidase in purine catabolism. This is not evidence for the reverse Cannabis conversion. Figure assays use NADPH, so exact NADH/NADPH specificity also needs review.',
            'status': 'direction-and-cofactor-review-not-new-Cannabis-activity'},
        'claim_boundary': 'Selected-certificate direction-risk audit, not a thermodynamic proof, gene-essentiality result or physiological direction assignment. Opposite-to-source-written is an exact orientation comparison only. Original candidate certificates, counts and equations are unchanged. No atom tracing.'}


def run():
    source = Path('data/reports/phase1-expanded-candidate-net.json')
    directions = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    report = json.loads(source.read_text()); families = direction_families(directions.read_text())
    masters = sorted({families[s['source_reaction_id']]['RHEA_ID_MASTER'] for r in report['reactions'] if r['new_catalog_candidate'] for s in r['sources']})
    lookups = fetch(masters)
    annotations = annotations_from_snapshots([], lookups)['source_annotations']
    output = build(report, families, annotations); output['lookups'] = lookups
    output['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, directions)}
    payload = json.dumps(output, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-candidate-direction-review.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    rows = [('metadata', 'review', {k: v for k, v in output.items() if k not in ('reviews', 'targets')})]
    rows += [('direction_review', r['id'], r) for r in output['reviews']]
    rows += [('target', r['cannabisdb_id'], r) for r in output['targets']]
    count = write_rows(rows, sha, Path('data/derived/phase1-candidate-direction-review.ndjson'))
    print(json.dumps({'summary': output['summary'], 'sha256': sha, 'rows': count}))


if __name__ == '__main__':
    run()
