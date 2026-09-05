"""New exact-family reference discovery for previously unsearched weighted gaps."""
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from .phase1_reference_discovery import direction_families, exact_annotations
from .phase1_new_references import attach
from .phase1_new_protein_search import run as screen, export_table


def queue(source, families):
    reactions = {r['id']: r for r in source['reactions']}
    rows = []
    for gap in source['candidate_gaps']:
        if gap['prior_searches']:
            continue
        rid = gap['reaction_id']
        ids = sorted({s['source_reaction_id'] for s in gap['source_joins']})
        rows.append({'reaction_id': rid, 'left': reactions[rid]['left'], 'right': reactions[rid]['right'],
            'source_reaction_ids': ids, 'sources': gap['source_joins'],
            'rhea_families': {sid: families[sid] for sid in ids if sid in families},
            'selected_uses': gap['selected_uses'], 'target_ids': gap.get('target_ids', []),
            'priority_target_ids': gap.get('target_ids', []), 'hypothesis_ids': [],
            'probe_ids': sorted({u['probe_id'] for u in gap['selected_uses'] if 'probe_id' in u}),
            'priority_boundary': 'Selected weighted witnesses, not guaranteed improvement or necessity across all possible pathways.'})
    return rows


def run(source=Path('data/reports/phase1-evidence-weighted-routes.json'), prefix='phase1-weighted-gap'):
    directions = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    original = json.loads(source.read_text())
    for p, sha in original['source_sha256'].items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest() != sha:
            raise ValueError('Source lineage changed')
    rows = queue(original, direction_families(directions.read_text()))
    masters = sorted({f['RHEA_ID_MASTER'] for r in rows for f in r['rhea_families'].values()})
    if not masters:
        raise ValueError('No new exact reference families')
    raw = Path('data/raw', prefix + '-search'); raw.mkdir(parents=True, exist_ok=True)
    expression = '(' + ' OR '.join(f'cc_catalytic_activity:"{rid.lower()}"' for rid in masters) + ') AND reviewed:true AND fragment:false'
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': expression, 'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
    cache = raw / 'lookup.json'
    if not cache.exists():
        with urllib.request.urlopen(url, timeout=45) as response:
            payload = response.read()
        exact_annotations(payload.decode(), set(masters))
        snapshot = raw / 'reference-annotations.tsv'; snapshot.write_bytes(payload)
        lookup = {'url': url, 'requested_master_ids': masters, 'status': 'retrieved',
            'snapshot': str(snapshot), 'sha256': hashlib.sha256(payload).hexdigest(),
            'retrieved_at': datetime.now(timezone.utc).isoformat()}
        cache.write_text(json.dumps(lookup) + '\n')
    lookup = json.loads(cache.read_text())
    if lookup['url'] != url or lookup['requested_master_ids'] != masters:
        raise ValueError('Cached lookup request mismatch')
    proteins = attach(rows, [lookup])
    report = {'schema': 'cannabis-carbon.' + prefix + '-references.v1', 'rows': rows,
        'reference_proteins': list(proteins.values()), 'lookups': [lookup],
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, directions)},
        'claim_boundary': 'Previously unsearched equations only. Exact published Rhea families supply reviewed nonfragment reference annotation leads, not characterized Cannabis activity, physiological direction or net pathway proof.'}
    discovery = Path('data/reports', prefix + '-references.json')
    discovery.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    output = Path('data/reports', prefix + '-search.json')
    print('Reference proteins:', len(proteins), flush=True)
    screen(discovery, raw, output)
    export_table(output, Path('data/derived', prefix + '-search.ndjson'))


if __name__ == '__main__':
    run()
