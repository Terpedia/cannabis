"""New reference discovery for replacement witnesses; retain all prior screens."""
import copy
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from .phase1_new_references import attach
from .phase1_purine_gap_references import PRIOR, queue as new_queue
from .phase1_reference_discovery import direction_families, exact_annotations
from .phase1_scope import write_rows


def queue(source, searches, families):
    reactions = {r['id']: r for r in source['reactions']}
    audit, unseen = [], []
    for gap in source['candidate_gaps']:
        rid = gap['reaction_id']
        prior = [{'report': path, 'row': copy.deepcopy(row)} for path, report in searches.items()
                 for row in report.get('rows', []) if row.get('reaction_id') == rid]
        audit.append({**copy.deepcopy(gap), 'prior_screens': prior,
                      'disposition': 'retain-prior-evidence-without-promotion' if prior else 'new-reference-discovery'})
        if not prior:
            unseen.append({**copy.deepcopy(gap), 'left': reactions[rid]['left'], 'right': reactions[rid]['right']})
    if len({g['reaction_id'] for g in audit}) != len(audit):
        raise ValueError('Duplicate gap')
    rows, _ = new_queue({'catalog_candidate_gaps': unseen, 'focused_targets': source['focused_targets']}, {}, families)
    return rows, audit


def run():
    names = (*PRIOR, 'phase1-purine-gap-search', 'phase1-deferred-search')
    paths = [Path('data/reports', n + '.json') for n in names]
    searches = {str(p): json.loads(p.read_text()) for p in paths}
    source = Path('data/reports/phase1-decay-sensitivity.json')
    directions = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    rows, audit = queue(json.loads(source.read_text()), searches, direction_families(directions.read_text()))
    masters = sorted({f['RHEA_ID_MASTER'] for r in rows for f in r['rhea_families'].values()})
    raw = Path('data/raw/phase1-replacement-references'); raw.mkdir(parents=True, exist_ok=True)
    lookups = []
    for start in range(0, len(masters), 25):
        batch = masters[start:start + 25]
        expression = '(' + ' OR '.join(f'cc_catalytic_activity:"{rid.lower()}"' for rid in batch) + ') AND reviewed:true AND fragment:false'
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': expression, 'format': 'tsv', 'fields': 'accession,rhea,organism_name,protein_name'})
        cache = raw / (hashlib.sha256(url.encode()).hexdigest() + '.json')
        if cache.exists():
            lookup = json.loads(cache.read_text())
            if lookup['url'] != url or lookup['requested_master_ids'] != batch:
                raise ValueError('Request cache mismatch')
        else:
            lookup = {'url': url, 'requested_master_ids': batch, 'retrieved_at': datetime.now(timezone.utc).isoformat()}
            try:
                with urllib.request.urlopen(url, timeout=45) as response:
                    payload = response.read()
                exact_annotations(payload.decode(), set(batch))
                sha = hashlib.sha256(payload).hexdigest()
                snapshot = raw / (sha + '.tsv'); snapshot.write_bytes(payload)
                lookup.update(status='retrieved', snapshot=str(snapshot), sha256=sha)
            except (OSError, ValueError) as error:
                lookup.update(status='retrieval-failed', reason=str(error))
            cache.write_text(json.dumps(lookup, separators=(',', ':')) + '\n')
        lookups.append(lookup)
    proteins = attach(rows, lookups)
    used = {m['accession'] for r in rows for m in r['reference_matches']}
    report = {'schema': 'cannabis-carbon.phase1-replacement-references.v1', 'rows': rows,
              'prior_screen_audit': audit, 'search_report_paths': list(searches), 'lookups': lookups,
              'reference_proteins': [proteins[a] for a in sorted(used)],
              'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [source, directions, *paths]},
              'summary': {'gap_equations': len(audit), 'new_equation_gaps': len(rows),
                          'retained_prior_equations': len(audit)-len(rows), 'requested_master_families': len(masters),
                          'failed_lookups': sum(l['status'] != 'retrieved' for l in lookups),
                          'equations_with_references': sum(bool(r['reference_matches']) for r in rows),
                          'reference_proteins': len(used)},
              'claim_boundary': 'Only equations with no prior search row enter new discovery. All prior results, including withheld partial-reference candidates, remain unchanged and are not promoted. Historical skipped searches remain distinguishable in their source reports. Reviewed references are annotation leads, not Cannabis activity, physiological direction or demonstrated routes. Atom tracing remains deferred.'}
    payload = json.dumps(report, separators=(',', ':')) + '\n'
    Path('data/reports/phase1-replacement-references.json').write_text(payload)
    sha = hashlib.sha256(payload.encode()).hexdigest()
    groups = [('equation_gap', 'rows', 'reaction_id'), ('prior_screen_audit', 'prior_screen_audit', 'reaction_id'),
              ('reference_protein', 'reference_proteins', 'accession'), ('lookup', 'lookups', 'url')]
    records = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, key, identifier in groups:
        records.extend((kind, r[identifier], r) for r in report[key])
    count = write_rows(records, sha, Path('data/derived/phase1-replacement-references.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
