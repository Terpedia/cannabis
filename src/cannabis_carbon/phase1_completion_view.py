"""Compact target-only view of the separately labeled completion hypotheses."""
import hashlib
import json
from pathlib import Path


def build(report, audit):
    ids = {hid for t in report['targets'] for hid in t['completion_ids']}
    completions = [h for h in report['completions'] if h['id'] in ids]
    variants = [v for v in report['variants'] if v['id'] in {h['variant_id'] for h in completions}]
    source_ids = {sid for v in variants for sid in v['source_record_ids']}
    references = {x['reference_reaction_id'] for h in completions for x in h['reference_templates']}
    compounds = {p['compound_id'] for h in completions for side in ('left', 'right') for p in h[side]}
    compounds.update(x['reference_product_id'] for h in completions for x in h['reference_templates'])
    return {'schema': 'cannabis-carbon.completion-view.v1', 'summary': report['summary'],
        'claim_boundary': report['claim_boundary'], 'targets': report['targets'], 'completions': completions,
        'variants': variants, 'compounds': [c for c in report['compounds'] if c['id'] in compounds],
        'reference_reactions': [{'id': r['id'], 'sources': r['sources']} for r in report['reference_reactions'] if r['id'] in references],
        'MARTS_sources': [s for s in audit['source_ledger'] if s['id'] in source_ids]}


def run():
    paths = [Path('data/reports/phase1-marts-completions.json'), Path('data/reports/phase1-marts-audit.json')]
    report, audit = [json.loads(p.read_text()) for p in paths]
    if report['source_sha256'][str(paths[1])] != hashlib.sha256(paths[1].read_bytes()).hexdigest():
        raise ValueError('Completion audit input mismatch')
    payload = json.dumps(build(report, audit), separators=(',', ':')) + '\n'
    folder = Path('docs/data/completion-view'); folder.mkdir(parents=True, exist_ok=True)
    (folder / 'bundle.json').write_text(payload)
    manifest = {'file': 'bundle.json', 'bytes': len(payload.encode()),
        'sha256': hashlib.sha256(payload.encode()).hexdigest(),
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (folder / 'index.json').write_text(json.dumps(manifest, separators=(',', ':')) + '\n')
    print(json.dumps(manifest))


if __name__ == '__main__':
    run()
