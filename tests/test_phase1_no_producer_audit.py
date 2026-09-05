import copy
import gzip
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from rdkit import rdBase
from cannabis_carbon.phase1_no_producer_audit import build, diagnostic_keys
from cannabis_carbon.cannabisdb_xml import _value

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((ROOT / path).read_text())


def test_diagnostics_preserve_isotopes_and_do_not_equate_stereochemical_identity():
    left, right = 'N[C@@H](C)C(=O)O', 'N[C@H](C)C(=O)O'
    assert diagnostic_keys(left)['stereo_removed'] == diagnostic_keys(right)['stereo_removed']
    assert diagnostic_keys(left)['uncharger'] != diagnostic_keys(right)['uncharger']
    assert diagnostic_keys('[13CH3]CO') != diagnostic_keys('CCO')
    assert diagnostic_keys('CC(=O)[O-]')['uncharger'] == diagnostic_keys('CC(=O)O')['uncharger']
    assert diagnostic_keys('CC(=O)[O-]')['stereo_removed'] != diagnostic_keys('CC(=O)O')['stereo_removed']


def test_all_no_producer_targets_and_exact_catalog_classification():
    report = read('data/reports/phase1-no-producer-audit.json')
    model = read('data/reports/phase1-replacement-candidate-net.json')
    network = read('data/reports/phase1-full-balanced-network.json')
    expected = {t['cannabisdb_id'] for t in model['scenarios'][0]['targets'] if t['net_status'] == 'no-net-producing-candidate-equation'}
    assert len(report['targets']) == len(expected) == 5897
    assert {t['cannabisdb_id'] for t in report['targets']} == expected
    participation, producers = defaultdict(set), defaultdict(set)
    for reaction in network['reactions']:
        net = defaultdict(Fraction)
        for side, sign in [('left', -1), ('right', 1)]:
            for member in reaction[side]:
                participation[member['compound_id']].add(reaction['id'])
                net[member['compound_id']] += sign * Fraction(member['coefficient'])
        for cid, amount in net.items():
            if amount:
                producers[cid].add(reaction['id'])
    original = {c['id']: c for c in network['compounds']}
    for row in report['targets']:
        cid = row['compound_id']
        assert row['canonical_smiles'] == original[cid]['smiles']
        assert row['exact_catalog_participation_reaction_ids'] == sorted(participation[cid])
        assert row['exact_catalog_net_producer_reaction_ids'] == sorted(producers[cid])
        assert not producers[cid] & model['candidate_reaction_evidence_ids'].keys()
        for kind, matches in row['diagnostic_alternatives'].items():
            assert cid not in matches and matches == sorted(set(matches))
            for other in matches:
                assert other in participation
                assert diagnostic_keys(original[cid]['smiles'])[kind] == diagnostic_keys(original[other]['smiles'])[kind]
    assert Counter(t['status'] for t in report['targets']) == report['summary']['status_counts']
    assert len(report['source_identity_conflicts']) == 71
    assert report['source_identity_summary']['compared_accessions'] == 6220
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha


def test_report_replay_has_no_source_mutation():
    report = read('data/reports/phase1-no-producer-audit.json')
    inputs = [read(p) for p in report['source_sha256'] if p.endswith('.json')]
    before = copy.deepcopy(inputs)
    result = build(*inputs)
    assert result['rdkit_version'] == rdBase.rdkitVersion
    assert {k: v for k, v in result.items() if k != 'rdkit_version'} == {
        k: v for k, v in report.items() if k not in ('source_sha256', 'rdkit_version')}
    assert inputs == before


def test_xml_assertions_match_archived_source_not_inferred_name_repairs():
    table = {r['accession']: r for r in read('data/terpedia/cannabisdb-compounds.json')['rows']}
    text = gzip.decompress((ROOT / 'data/terpedia/cannabisdb-compounds.xml.gz').read_bytes()).decode()
    seen = set()
    for match in re.finditer(r'<compound>.*?</compound>', text, re.S):
        element = ET.fromstring(match.group(0))
        accession = element.findtext('accession')
        assert accession not in seen
        seen.add(accession)
        for key, tag in [('smiles', 'smiles'), ('formula', 'chemical_formula'), ('inchikey', 'inchikey'), ('name', 'name')]:
            assert table[accession][key] == _value(element.findtext(tag))
    assert seen == set(table)
    conflicts = {r['cannabisdb_id']: r for r in read('data/reports/phase1-no-producer-audit.json')['source_identity_conflicts']}
    glycerol = conflicts['CDB006156']
    assert glycerol['xml_assertion']['smiles'] == 'OCC(O)CO'
    assert glycerol['xml_assertion']['formula'] == 'C3H8O3'
    assert glycerol['sdf_derived_assertion']['formula'] == 'C14H22O'
    assert glycerol['sdf_derived_assertion']['smiles'] != glycerol['xml_assertion']['smiles']
