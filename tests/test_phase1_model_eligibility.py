import copy
import pytest
from cannabis_carbon.phase1_screened_overlay import build_overlay, apply_overlay


def fixture():
    parent = {'reactions': [{'id': 'r'}], 'hypotheses': []}
    hit = {'id': 'a', 'cannabis_accession': 'p', 'reference_accession': 'ref',
           'passes_screen': True, 'identity_percent': 90, 'query_coverage_percent': 100,
           'reference_coverage_percent': 100, 'evalue': 1e-30, 'bitscore': 100}
    row = {'reaction_id': 'r', 'hypothesis_ids': [], 'passing_alignment_ids': ['a'],
           'screened_cannabis_proteins': ['p'], 'reference_matches': [{'accession': 'ref'}],
           'validation_blockers': [], 'proposed_test': 'Assay exact substrates'}
    return parent, {'passing_alignments': [hit], 'rows': [row], 'screen': {}}


@pytest.mark.parametrize('level', ['report', 'row', 'reference'])
def test_explicit_ineligibility_overrides_passing_homology(level):
    parent, search = fixture()
    assert len(build_overlay(parent, search)['enzyme_evidence']) == 1
    target = {'report': search, 'row': search['rows'][0],
              'reference': search['rows'][0]['reference_matches'][0]}[level]
    target['model_eligible'] = False
    before = copy.deepcopy(search)
    with pytest.raises(ValueError, match='ineligible'):
        build_overlay(parent, search)
    assert search == before


@pytest.mark.parametrize('level', ['report', 'evidence'])
def test_prebuilt_overlay_cannot_bypass_review_rejection(level):
    overlay = {'enzyme_evidence': [{'id': 'ev', 'reaction_id': 'r'}]}
    (overlay if level == 'report' else overlay['enzyme_evidence'][0])['model_eligible'] = False
    with pytest.raises(ValueError, match='ineligible'):
        apply_overlay({}, overlay)
