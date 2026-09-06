import hashlib
import json
from pathlib import Path
import pytest
from rdkit import Chem
from cannabis_carbon.phase1_selenium_source_identity import IDS, parse
from cannabis_carbon.phase1_marts_completions import balanced


def sources():
    folder = Path('data/raw/selenium-kegg')
    return (folder/'R03601.txt').read_text(), {k:(folder/(k+'.mol')).read_text() for k in IDS}


def test_pinned_source_and_balance():
    report = json.loads(Path('data/reports/phase1-selenium-source-identity.json').read_bytes())
    for row in report['retrievals'].values():
        assert hashlib.sha256(Path(row['path']).read_bytes()).hexdigest() == row['sha256']
    for path, sha in report['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    fields, structures, compounds, sides = parse(*sources())
    assert fields == report['source_fields']
    assert balanced(sides, compounds)
    assert sides == [report['reaction']['left'], report['reaction']['right']]
    assert list(compounds.values()) == report['compounds']
    assert structures == {r['kegg_id']:r['smiles'] for r in report['endpoints']}
    assert 'sulfide' in fields['NAME'] and 'selenide' in fields['DEFINITION']
    assert report['summary']['coverage_gain_claimed'] == 0
    assert structures['C05688'] == 'N[C@@H](C[SeH])C(=O)O'


def test_invalid_source_and_sulfur_substitution_rejected():
    reaction, blocks = sources()
    with pytest.raises(ValueError, match='equation changed'):
        parse(reaction.replace('C01528', 'C00001'), blocks)
    with pytest.raises(ValueError, match='Invalid concrete'):
        parse(reaction, {**blocks, 'C05688':'invalid'})
    sulfur = Chem.MolToMolBlock(Chem.MolFromSmiles('N[C@@H](CS)C(=O)O'))
    with pytest.raises(ValueError, match='not element/isotope/charge balanced'):
        parse(reaction, {**blocks, 'C05688':sulfur})
