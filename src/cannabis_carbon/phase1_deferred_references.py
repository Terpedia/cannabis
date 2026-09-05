"""Recover skipped reference discovery without rewriting historical search results."""
import copy
import hashlib
import json
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .phase1_new_references import attach
from .phase1_reference_discovery import direction_families, exact_annotations
from .phase1_scope import write_rows


def queue(precursors, candidates, searches, discoveries, families):
    audit, rows = [], []
    for gap in precursors['catalog_candidate_gaps']:
        rid = gap['reaction_id']
        if rid in candidates['candidate_reaction_evidence_ids']:
            continue
        previous = []
        for name, search in searches.items():
            discovery_name = search['source_discovery']
            source_rows = {r['reaction_id']: r for r in discoveries[discovery_name]['rows']}
            for r in search['rows']:
                if r['reaction_id'] != rid:
                    continue
                if r['passing_alignment_ids']:
                    raise ValueError('Prior candidate omitted from model')
                source = source_rows[rid]
                previous.append({'search_report': name, 'discovery_report': discovery_name,
                    'search_status': r['search_status'], 'lookup_status': source.get('lookup_status'),
                    'reference_sequences_present': r['reference_sequences_present'],
                    'reference_match_accessions': [m['accession'] for m in source['reference_matches']]})
        deferred = bool(previous) and all(p['search_status'] == 'no-reference-sequence' and
            p['lookup_status'] == 'not-searched-in-priority-pass' and not p['reference_sequences_present'] and
            not p['reference_match_accessions'] for p in previous)
        audit.append({**copy.deepcopy(gap), 'prior_attempts': previous,
            'disposition': 'recover-skipped-reference-discovery' if deferred else 'retain-other-prior-evidence'})
        if not deferred:
            continue
        source_ids = sorted({s['source_reaction_id'].upper() for s in gap['source_joins']})
        probes = {u['probe_compound_id'] for u in gap['selected_uses']}
        targets = sorted({t['cannabisdb_id'] for t in precursors['focused_targets'] if t['compound_id'] in probes})
        rows.append({'reaction_id': rid, 'left': copy.deepcopy(gap['left']), 'right': copy.deepcopy(gap['right']),
            'sources': copy.deepcopy(gap['source_joins']), 'source_reaction_ids': source_ids,
            'selected_uses': copy.deepcopy(gap['selected_uses']), 'selected_probe_compound_ids': sorted(probes),
            'target_ids': targets, 'priority_target_ids': targets, 'hypothesis_ids': [],
            'rhea_families': {sid: families[sid] for sid in source_ids if sid in families},
            'priority_boundary': 'Previously skipped reference discovery; selected certificate membership is not reaction necessity.'})
    if len({r['reaction_id'] for r in audit}) != len(audit):
        raise ValueError('Duplicate remaining gap')
    return rows, audit


def run():
    # Pin the inventory of prior reports; do not let this module's later output
    # silently become an input and erase the distinction on replay.
    names = ('phase1-new-protein-search', 'phase1-route-protein-search', 'phase1-catalog-protein-search',
             'phase1-backfill-protein-search', 'phase1-completion-protein-search', 'phase1-archived-protein-search',
             'phase1-plant-purine-search', 'phase1-purine-gap-search')
    search_paths = [Path('data/reports', n + '.json') for n in names]
    searches = {str(p): json.loads(p.read_text()) for p in search_paths}
    discoveries = {s['source_discovery']: json.loads(Path(s['source_discovery']).read_text()) for s in searches.values()}
    for s in searches.values():
        if hashlib.sha256(Path(s['source_discovery']).read_bytes()).hexdigest() != s['source_discovery_sha256']:
            raise ValueError('Prior discovery checksum mismatch')
    sources = [Path('data/reports/phase1-purine-precursor-audit.json'), Path('data/reports/phase1-purine-candidate-net.json')]
    directions = Path('data/raw/phase1-reference-discovery/rhea-directions.tsv')
    rows, audit = queue(*[json.loads(p.read_text()) for p in sources], searches, discoveries, direction_families(directions.read_text()))
    masters = sorted({f['RHEA_ID_MASTER'] for r in rows for f in r['rhea_families'].values()})
    raw = Path('data/raw/phase1-deferred-references'); raw.mkdir(parents=True, exist_ok=True)
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
    for row in rows:
        if not row['rhea_families']:
            row['lookup_status'] = 'no-published-Rhea-family-mapping'
    used = {m['accession'] for r in rows for m in r['reference_matches']}
    paths = sources + search_paths + [Path(p) for p in discoveries] + [directions]
    report = {'schema': 'cannabis-carbon.phase1-deferred-references.v1', 'rows': rows, 'prior_attempt_audit': audit,
        'reference_proteins': [proteins[a] for a in sorted(used)], 'lookups': lookups,
        'search_report_paths': list(searches), 'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'summary': {'remaining_gap_equations': len(audit), 'skipped_reference_discovery_equations': len(rows),
            'requested_master_families': len(masters), 'failed_lookups': sum(l['status'] != 'retrieved' for l in lookups),
            'equations_with_reference_leads': sum(bool(r['reference_matches']) for r in rows),
            'distinct_reference_proteins': len(used), 'lookup_status_counts': dict(Counter(r['lookup_status'] for r in rows))},
        'claim_boundary': 'A historical no-reference-sequence row is not necessarily a completed reference search. This audit recovers only equations whose every prior attempt explicitly skipped discovery. Other negative, weak, missing-reference and incomplete results remain unchanged. Reviewed annotations are leads, not Cannabis activity, physiological direction or pathway completion. Atom tracing remains deferred.'}
    payload = json.dumps(report, separators=(',', ':')) + '\n'; sha = hashlib.sha256(payload.encode()).hexdigest()
    Path('data/reports/phase1-deferred-references.json').write_text(payload)
    groups = [('equation_gap', 'rows', 'reaction_id'), ('prior_attempt_audit', 'prior_attempt_audit', 'reaction_id'), ('reference_protein', 'reference_proteins', 'accession'), ('lookup', 'lookups', 'url')]
    records = [('metadata', 'report', {k: v for k, v in report.items() if k not in {g[1] for g in groups}})]
    for kind, key, id_key in groups:
        records.extend((kind, r[id_key], r) for r in report[key])
    count = write_rows(records, sha, Path('data/derived/phase1-deferred-references.ndjson'))
    print(json.dumps({'summary': report['summary'], 'sha256': sha, 'rows': count}), flush=True)


if __name__ == '__main__':
    run()
