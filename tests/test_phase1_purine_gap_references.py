import copy
import hashlib
import json
from pathlib import Path

import pytest

from cannabis_carbon.phase1_new_references import attach
from cannabis_carbon.phase1_purine_gap_references import PRIOR, queue
from cannabis_carbon.phase1_reference_discovery import direction_families

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((ROOT / 'data/reports' / (name + '.json')).read_text())


def test_all_gaps_prior_evidence_and_exact_reference_joins_replay(monkeypatch):
    monkeypatch.chdir(ROOT)
    source = read('phase1-purine-precursor-audit')
    prior = {str(Path('data/reports', n + '.json')): read(n) for n in PRIOR}
    report = read('phase1-purine-gap-references')
    families = direction_families((ROOT / 'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text())
    before = copy.deepcopy((source, prior))
    rows, audit = queue(source, prior, families)
    assert audit == report['prior_screen_audit']
    proteins = attach(rows, report['lookups'])
    for row in rows:
        if not row['rhea_families']:
            row['lookup_status'] = 'no-published-Rhea-family-mapping'
    assert rows == report['rows']
    used = {m['accession'] for row in rows for m in row['reference_matches']}
    assert [proteins[a] for a in sorted(used)] == report['reference_proteins']
    assert (source, prior) == before
    assert len(audit) == 70 and len(rows) == 22
    assert sum(bool(r['prior_screens']) for r in audit) == 48
    assert len(used) == 833
    assert sum(bool(r['reference_matches']) for r in rows) == 16
    for p in report['reference_proteins']:
        assert p['evidence_status'] == 'reviewed-reference-annotation-not-Cannabis-activity'
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    for lookup in report['lookups']:
        assert 'reviewed%3Atrue' in lookup['url']
        assert lookup['status'] == 'retrieved'
        assert hashlib.sha256((ROOT / lookup['snapshot']).read_bytes()).hexdigest() == lookup['sha256']
    source['catalog_candidate_gaps'].append(copy.deepcopy(source['catalog_candidate_gaps'][0]))
    with pytest.raises(ValueError, match='Duplicate candidate gap'):
        queue(source, prior, families)


def test_prior_passing_candidate_cannot_be_silently_ignored():
    source = {'catalog_candidate_gaps': [{'reaction_id': 'R'}]}
    prior = {'prior': {'rows': [{'reaction_id': 'R', 'passing_alignment_ids': ['A']}]}}
    with pytest.raises(ValueError, match='Passing prior candidate'):
        queue(source, prior, {})
