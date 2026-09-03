from cannabis_carbon.atom_mapping import map_reaction_smiles


def test_oxidation_preserves_all_product_carbons():
    result = map_reaction_smiles("CCO>>CC=O")
    assert result["status"] == "inferred"
    assert len(result["mappings"]) == 2
    assert result["unresolved_product_carbons"] == []


def test_carboxylation_like_new_carbon_is_unresolved():
    result = map_reaction_smiles("CCO>>CC(=O)O")
    assert result["status"] == "unresolved"
    assert len(result["unresolved_product_carbons"]) >= 1
