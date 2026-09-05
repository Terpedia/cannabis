"""Add screened all-target diagnostic reactions without promoting chemical-only routes."""
import hashlib
import json
from pathlib import Path
from .phase1_thiolase_candidate_net import build as extend_net
from .phase1_scope import write_rows


def build(network, parent, search, links):
    result = extend_net(network, parent, search, links, search_filename='phase1-remaining-gap-search.json')
    result['schema'] = 'cannabis-carbon.phase1-remaining-candidate-net.v1'
    for reaction in result['reactions']:
        reaction['new_remaining_gap_candidate'] = reaction.pop('new_thiolase_candidate')
    result['claim_boundary'] = ('All 6,220 historical exact target records and eight distinct synthase participant probes retained. '
        'Adds only passing reviewed-reference homology candidates from the all-target weighted gap search, not characterized Cannabis enzymes. '
        'Previous candidate evidence identifiers and eight reverse exclusions are preserved; other directions remain hypotheses. '
        'CO2 is the only net carbon input; regenerated pre-existing pools are allowed, not zero-pool startup. '
        'No protonation or stereochemical identity merges, no physiological flux or compartment claims, no atom tracing.')
    return result


def run():
    paths = [Path('data/reports', n + '.json') for n in ('phase1-full-balanced-network',
        'phase1-thiolase-candidate-net', 'phase1-remaining-gap-search', 'phase1-synthase-reaction-links')]
    inputs = [json.loads(p.read_text()) for p in paths]
    for report in inputs:
        for p, sha in report.get('source_sha256', {}).items():
            if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
                raise ValueError('Source lineage changed')
        if 'source_discovery' in report and hashlib.sha256(Path(report['source_discovery']).read_bytes()).hexdigest() != report['source_discovery_sha256']:
            raise ValueError('Discovery lineage changed')
    report = build(*inputs)
    report['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-remaining-candidate-net.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('reaction', 'reactions', 'id'), ('compound', 'compounds', 'id'),
        ('enzyme_evidence', 'enzyme_evidence', 'id'), ('synthase_reference_link', 'synthase_reference_links', 'id'),
        ('probe_result', 'probe_results', 'id')]
    rows = [('metadata', 'report', {k: v for k, v in report.items() if k not in {'scenarios', *(g[1] for g in groups)}})]
    for kind, collection, key in groups:
        rows.extend((kind, r[key], r) for r in report[collection])
    for scenario in report['scenarios']:
        rows.append(('scenario', scenario['id'], {k: v for k, v in scenario.items() if k not in ('targets', 'certificates')}))
        for kind, collection, key in [('target', 'targets', 'cannabisdb_id'), ('certificate', 'certificates', 'compound_id')]:
            rows.extend((kind, scenario['id'] + ':' + r[key], r) for r in scenario[collection])
    count = write_rows(rows, sha, Path('data/derived/phase1-remaining-candidate-net.ndjson'))
    print(json.dumps({'summary': report['summary'], 'rows': count, 'sha256': sha}), flush=True)


if __name__ == '__main__':
    run()
