import pytest
import hashlib
from cannabis_carbon.phase1_reference_discovery import exact_annotations, attach_families


def test_reference_query_results_require_explicit_returned_reaction_ids():
    tsv = 'Entry\tRhea ID\tOrganism\tProtein names\nP1\tRHEA:10020 RHEA:10021\tPlant\tEnzyme\nP2\tRHEA:10022\tOther\tOther enzyme\n'
    result = exact_annotations(tsv, {'RHEA:10020', 'RHEA:10023'})
    assert result['RHEA:10020'][0]['accession'] == 'P1'
    assert result['RHEA:10020'][0]['annotated_rhea_ids'] == ['RHEA:10020', 'RHEA:10021']
    assert result['RHEA:10023'] == []  # Never infer family IDs by arithmetic.
    with pytest.raises(ValueError, match='schema'):
        exact_annotations('<html>error</html>', {'RHEA:10020'})


def test_direction_family_is_candidate_not_exact_direction(tmp_path):
    tsv = b'Entry\tRhea ID\tOrganism\tProtein names\nP1\tRHEA:10020\tPlant\tEnzyme\n'
    snapshot = tmp_path / 'uniprot.tsv'
    snapshot.write_bytes(tsv)
    report = {'uniprot_lookups': [{'status': 'retrieved', 'url': 'https://example.org/query', 'snapshot': str(snapshot), 'sha256': hashlib.sha256(tsv).hexdigest()}],
              'summary': {}, 'rows': [{'reaction_id': rid, 'reference_annotations': []} for rid in ['RHEA:10021', 'RHEA:10024']]}
    mapping = 'RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n10020\t10021\t10022\t10023\n'
    result = attach_families(report, mapping)
    assert result['rows'][0]['reference_annotations'] == []
    candidate = result['rows'][0]['family_reference_annotations'][0]
    assert candidate['accession'] == 'P1'
    assert candidate['direction_status'] == 'not-established-for-requested-direction'
    assert result['rows'][1]['family_reference_annotations'] == []
    snapshot.write_bytes(b'changed')
    with pytest.raises(ValueError, match='checksum'):
        attach_families(report, mapping)
