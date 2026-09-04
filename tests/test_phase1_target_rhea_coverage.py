import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

from cannabis_carbon.balance import _reaction_smiles_balance
from cannabis_carbon.phase1_target_coverage import audit_targets, encoded_structure
from cannabis_carbon.phase1_target_rhea_coverage import source_records


def test_filtered_ledger_counts_all_equations_and_retains_all_targets():
    sources = [{'id': 'test:1', 'source_reaction_id': '1', 'reaction_smiles': 'C=C.O>>CCO'},
               {'id': 'test:2', 'source_reaction_id': '2', 'reaction_smiles': 'C>>C'},
               {'id': 'test:3', 'source_reaction_id': '3', 'reaction_smiles': 'CO.*>>CO'}]
    result = audit_targets([{'id': 'A', 'smiles': 'CCO'}, {'id': 'B', 'smiles': 'CCC'}],
                           {'reactions': []}, {'reactions': []}, sources, True)
    assert result['summary']['source_reaction_records'] == 3
    assert result['summary']['reaction_balance_status_counts'] == {'balanced': 2, 'not-auditable': 1}
    assert [r['id'] for r in result['reaction_ledger']] == ['test:1']
    assert len(result['targets']) == 2
    assert result['targets'][1]['coverage_status'] == 'no-exact-encoded-reaction-match'


def test_source_ids_preserve_distinct_equations_and_evidence_boundary():
    rows = [dict(rule_id='RHEA:1', reaction_smarts=s, source_url='source',
                 source_download_url='download', direction_mode='recorded', source_evidence_type='curated')
            for s in ('C>>C', 'CC>>CC')]
    sources = source_records(rows)
    assert sources[0]['id'] != sources[1]['id']
    assert all('not established Cannabis' in s['orientation_boundary'] for s in sources)


def test_full_catalog_report_provenance_participation_and_balance():
    root = Path(__file__).resolve().parents[1]
    path = root / 'data/reports/phase1-target-rhea-coverage.json'
    report = json.loads(path.read_text())
    assert path.read_bytes() == (root / 'docs/data/phase1-target-rhea-coverage.json').read_bytes()
    targets = json.loads((root / 'docs/data/compounds.json').read_text())['compounds']
    assert [r['id'] for r in targets] == [r['cannabisdb_id'] for r in report['targets']]
    for filename, digest in report['source_sha256'].items():
        source = root / filename
        # Raw GCP snapshot is not shipped to CI; its checksum remains in provenance.
        if source.exists():
            assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
    assert sum(report['summary']['reaction_balance_status_counts'].values()) == report['summary']['source_reaction_records']
    ledger = {r['id']: r for r in report['reaction_ledger']}
    assert len(ledger) == report['summary']['retained_matching_reaction_records']
    snapshot = root / report['catalog_provenance']['snapshot']
    if snapshot.exists():
        sources = source_records(json.loads(snapshot.read_text()))
        assert len(sources) == report['summary']['additional_rhea_source_records']
        source_index = {r['id']: r for r in sources}
        for reaction in ledger.values():
            if reaction['source_layer'] == 'terpedia-full-rhea-catalog':
                assert all(reaction[k] == v for k, v in source_index[reaction['id']].items())
    canonical = lru_cache(None)(lambda s: encoded_structure(s)[0])
    sides = {}
    for rid, reaction in ledger.items():
        element, charge = _reaction_smiles_balance(reaction['reaction_smiles'])
        assert (element, charge) == (reaction['element_balance'], reaction['charge_balance'])
        status = 'balanced' if element and charge and element['status'] == charge['status'] == 'balanced' else 'imbalanced' if element and charge else 'not-auditable'
        assert reaction['computed_balance_status'] == status
        sides[rid] = [Counter(canonical(s) for s in side.split('.'))
                      for side in reaction['reaction_smiles'].split('>>')]
    referenced = set()
    for target in report['targets']:
        assert target['canonical_isomeric_smiles'] == canonical(target['source_smiles'])
        assert target['balanced_reaction_record_count'] == sum(m['computed_balance_status'] == 'balanced' for m in target['reaction_matches'])
        for match in target['reaction_matches']:
            rid = match['reaction_record_id']
            referenced.add(rid)
            assert match['computed_balance_status'] == ledger[rid]['computed_balance_status']
            for role in match['roles']:
                assert sides[rid][0 if role['equation_side'] == 'left' else 1][target['canonical_isomeric_smiles']] == role['coefficient']
    assert referenced == set(ledger)
