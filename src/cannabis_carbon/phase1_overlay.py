"""Compact, exact-variant Phase 1 evidence for the static map."""
import hashlib
import json
from collections import Counter
from pathlib import Path


def key(row):
    return row['reaction_id'], row['reaction_smarts']


def unique_rows(rows):
    result = {key(row): row for row in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate reaction-ID/SMARTS variant')
    return result


def build(queue_path, evidence_path, output):
    queue = json.loads(queue_path.read_text())
    evidence = json.loads(evidence_path.read_text())
    variants = unique_rows(queue['rows'])
    searched = unique_rows(evidence['rows'])
    balanced = {k for k, row in variants.items() if row['balance_status'] == 'balanced'}
    if balanced != searched.keys():
        raise ValueError('Evidence must cover exactly the balanced queue variants')
    rows = []
    for k, variant in variants.items():
        row = searched.get(k, {})
        if row and row['balance_status'] != 'balanced':
            raise ValueError('Search evidence contains a non-balanced variant')
        proteins = sorted({h['cannabis_accession'] for h in row.get('sequence_hits', []) if h['passes_screen']})
        associations = row.get('core_reaction_evidence', [])
        annotated = sorted({p for a in associations for p in a['enzyme_association_ids']})
        candidates = [p for a in associations for p in a.get('candidate_proteins', [])]
        if k not in balanced:
            status = 'balance-unresolved'
        elif proteins:
            status = 'screened-homology'
        elif annotated or candidates:
            status = 'core-association'
        elif row['search_status'] == 'no-reference-sequence':
            status = 'missing-reference'
        elif row['search_status'] == 'no-hits':
            status = 'no-hits'
        else:
            status = 'weak-hits'
        rows.append({
            'reaction_id': k[0], 'reaction_smarts': k[1],
            'balance_status': variant['balance_status'], 'evidence_status': status,
            'screened_proteins': proteins, 'core_enzyme_ids': annotated,
            'core_reaction_ids': sorted({a['core_reaction_id'] for a in associations}),
            'search_status': row.get('search_status', 'not-searched-balance-unresolved'),
            'source_urls': variant['source_urls'],
        })
    with_evidence = sum(r['evidence_status'] in ('screened-homology', 'core-association') for r in rows)
    result = {
        'schema': 'cannabis-carbon.phase1-map-evidence.v1',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (queue_path, evidence_path)},
        'summary': {
            'reaction_variants': len(rows), 'balanced_variants': len(balanced),
            'balance_status_counts': dict(Counter(r['balance_status'] for r in rows)),
            'evidence_status_counts': dict(Counter(r['evidence_status'] for r in rows)),
            'balanced_variants_with_candidate_enzyme_evidence': with_evidence,
            'balanced_variants_without_candidate_enzyme_evidence': len(balanced) - with_evidence,
        },
        'metric_scope': 'Three-hop candidate expansion only; distinct reaction-ID/SMARTS variants, not graph edges or total Cannabis metabolome coverage. Evidence categories are exclusive: screened homology takes priority over core association.',
        'claim_boundary': 'Balance and enzyme annotations/homology do not establish Cannabis reaction specificity, physiological direction, all-substrate availability, or a complete pathway. No confirmed-enzyme coverage is inferred. Atom tracing is deferred.',
        'rows': rows,
    }
    output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    return result['summary']


if __name__ == '__main__':
    print(build(Path('data/reports/phase1-enzyme-discovery-queue.json'),
                Path('data/reports/phase1-core-enzyme-evidence.json'),
                Path('data/reports/phase1-map-evidence.json')))
