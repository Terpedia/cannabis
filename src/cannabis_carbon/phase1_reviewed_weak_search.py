"""Apply a checksummed review as an additive layer over immutable search evidence."""
import copy
import hashlib
import json
from pathlib import Path


def build():
    source = Path('data/reports/phase1-weak-nonplant-search.json')
    review_path = Path('data/reports/phase1-weak-hit-domain-review.json')
    original = json.loads(source.read_text())
    review = json.loads(review_path.read_text())
    for path, digest in review['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('Domain review lineage changed')
    report = copy.deepcopy(original)
    matches = [h for h in report['passing_alignments'] if h['id'] == review['alignment']['id']]
    if len(matches) != 1 or matches[0] != review['alignment']:
        raise ValueError('Reviewed alignment differs from search evidence')
    matches[0]['model_eligible'] = review['decision']['model_eligible']
    matches[0]['curation_review'] = {'source_path': str(review_path),
        'source_sha256': hashlib.sha256(review_path.read_bytes()).hexdigest(), 'decision': review['decision']}
    report['schema'] = 'cannabis-reviewed-weak-search-v1'
    report['review_source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, review_path)}
    report['review_summary'] = {'screen_passing_alignments': len(report['passing_alignments']),
        'explicitly_rejected_alignments': sum(h.get('model_eligible') is False for h in report['passing_alignments']),
        'unreviewed_passing_alignments': sum('model_eligible' not in h for h in report['passing_alignments'])}
    report['claim_boundary'] += ' Additive domain review rejects the short accessory-domain hit for model integration. Original search counts and passes_screen flags are preserved, not reinterpreted as eligible enzyme counts.'
    return report


if __name__ == '__main__':
    Path('data/reports/phase1-reviewed-weak-search.json').write_text(json.dumps(build(), separators=(',', ':')) + '\n')
