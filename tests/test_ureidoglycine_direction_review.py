"""Keep literature triage tied to exact source orientation, not enzyme absence."""
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem

from cannabis_carbon.phase1_reference_discovery import direction_families

ROOT = Path(__file__).resolve().parents[1]


def test_review_matches_exact_catalog_direction_and_balanced_structures():
    review = json.loads((ROOT / 'data/curation/ureidoglycine-direction-review.json').read_text())
    for path, sha in review['source_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha
    source = json.loads((ROOT / 'data/reports/phase1-purine-precursor-audit.json').read_text())
    gap = next(g for g in source['catalog_candidate_gaps'] if g['reaction_id'] == review['reaction_id'])
    family = direction_families((ROOT / 'data/raw/phase1-reference-discovery/rhea-directions.tsv').read_text())[review['rhea_master']]
    assert family['RHEA_ID_LR'] == review['rhea_decay_direction']
    decay = next(j for j in gap['source_joins'] if j['source_reaction_id'] == review['rhea_decay_direction'])
    assert decay['source_left_corresponds_to'] == 'right'
    assert review['decay_canonical_direction'] == 'hypothetical-right-to-left'
    assert Counter(s['direction_mode'] for s in gap['selected_uses']) == {
        review['selected_canonical_direction']: review['selected_certificate_uses']}
    assert review['selected_canonical_direction'] == 'hypothetical-left-to-right'
    compounds = {c['id']: c for c in source['compounds']}

    def composition(side):
        totals = Counter()
        for term in gap[side]:
            molecule = Chem.AddHs(Chem.MolFromSmiles(compounds[term['compound_id']]['smiles']))
            for atom in molecule.GetAtoms():
                totals[(atom.GetSymbol(), atom.GetIsotope())] += term['coefficient']
                totals['charge'] += term['coefficient'] * atom.GetFormalCharge()
        return totals

    assert composition('left') == composition('right')
    assert {compounds[t['compound_id']]['smiles'] for t in gap['left']} == {
        'O=CC(=O)[O-]', '[NH4+]', 'NC(N)=O'}
    assert {compounds[t['compound_id']]['smiles'] for t in gap['right']} == {
        'NC(=O)N[C@H]([NH3+])C(=O)[O-]', 'O'}
    snapshot = next(p for p in review['source_sha256'] if p.endswith('.response.json'))
    triples = json.loads((ROOT / snapshot).read_text())['results']['bindings']
    assert any(t['s']['value'] == 'http://rdf.rhea-db.org/33871'
               and t['p']['value'] == 'http://rdf.rhea-db.org/equation'
               and t['o']['value'] == review['source_equation'] for t in triples)
    assert review['model_action'].startswith('review-only;')
