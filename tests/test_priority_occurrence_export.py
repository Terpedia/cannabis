import copy
import hashlib
import json
from pathlib import Path

import pytest
from cannabis_carbon.priority_occurrence_export import validate

ROOT = Path(__file__).resolve().parents[1]


def inputs():
    return [json.loads((ROOT / p).read_text()) for p in (
        'data/curation/priority-occurrence-review.json', 'data/terpedia/cannabisdb-compounds.json',
        'data/reports/phase1-identity-branches.json')]


def test_occurrence_export_preserves_complete_review_and_source_hashes():
    review, table, branches = inputs()
    before = copy.deepcopy(review)
    rows = validate(review, table, branches)
    assert len(rows) == 5
    assert [r[2] for r in rows[1:]] == review['rows']
    assert review == before
    for path, sha in review['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    acetamide = next(r for r in review['rows'] if r['cannabisdb_id'] == 'CDB000546')
    assert acetamide['occurrence_status'] == 'primary-paper-tentative-headspace-identification'
    assert [o['net_match_percent'] for o in acetamide['observations']] == [70, 87]
    assert all(o['pdf_page_one_based'] == 8 for o in acetamide['observations'])


@pytest.mark.parametrize('mutation', ['drop', 'duplicate', 'reference', 'snapshot'])
def test_occurrence_export_rejects_broken_provenance(mutation):
    review, table, branches = inputs()
    if mutation == 'drop':
        review['rows'].pop()
    elif mutation == 'duplicate':
        review['rows'][1] = review['rows'][0]
    elif mutation == 'reference':
        review['rows'][0]['source_reference_pmids'] = ['26657499']
    else:
        review['rows'][-1]['observations'][0]['snapshot'] = 'untracked.pdf'
    with pytest.raises(ValueError):
        validate(review, table, branches)
