from cannabis_carbon.terpedia import cytoscape_elements


def test_terpedia_reaction_becomes_compound_edge():
    network = {"entities": [
        {"id": "chebi:a", "type": "metabolite", "label": "A"},
        {"id": "chebi:b", "type": "metabolite", "label": "B"},
        {"id": "rhea:1", "type": "biochemical_reaction", "label": "A = B", "url": "https://rhea.example/1"},
        {"id": "protein:p", "type": "protein", "label": "enzyme P"},
    ], "statements": [
        {"subjectId": "rhea:1", "predicate": "has_reactant", "objectEntityId": "chebi:a"},
        {"subjectId": "rhea:1", "predicate": "has_product", "objectEntityId": "chebi:b"},
        {"subjectId": "protein:p", "predicate": "catalyzes", "objectEntityId": "rhea:1", "qualifiers": {"directExperimentalEvidence": False}},
    ]}
    graph = cytoscape_elements(network)
    assert graph["stats"] == {"metabolites": 2, "reaction_edges": 1, "reactions": 1}
    assert graph["edges"][0]["data"]["status"] == "candidate"
