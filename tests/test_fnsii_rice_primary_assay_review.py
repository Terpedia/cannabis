import json
from pathlib import Path


def test_rice_primary_assay_links_reference_but_does_not_assign_partner(monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parents[1])
    review = json.loads(Path('data/curation/fnsii-rice-primary-assay-review.json').read_text())
    annotation = json.loads(Path(review['annotation_snapshot']).read_text())
    assert review['reference_accession'] == annotation['primaryAccession'] == 'Q0JFI2'
    assert review['pmid'] in {r['citation']['id'] for r in annotation['references']}
    xref = next(r for r in annotation['uniProtKBCrossReferences'] if r['database'] == 'EMBL' and r['id'] == review['reported_cdna'])
    assert next(p['value'] for p in xref['properties'] if p['key'] == 'ProteinId') == review['uniprot_linked_protein'] == 'BAG94859.1'
    activity = next(c for c in annotation['comments'] if c['commentType'] == 'CATALYTIC ACTIVITY')
    assert activity['reaction']['ecNumber'] == '1.14.19.76'
    assert {'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': review['pmid']} in activity['reaction']['evidences']
    assert review['assay'] == {'temperature_C': 30, 'duration_min': 30, 'potassium_phosphate_mM': 100, 'pH': 8.0, 'L_glutathione_mM': 2, 'NADPH_mM': 5, 'substrate_uM': 100, 'microsomal_protein_ug': 200}
    assert review['model_eligible'] is False
    assert review['electron_transfer']['native_Cannabis_partner_assigned'] is False
    assert review['electron_transfer']['separately_supplied_reductase_named_in_reviewed_methods'] is False
