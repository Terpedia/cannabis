from cannabis_carbon.atom_mapping import map_reaction_smiles


def test_oxidation_preserves_all_product_carbons():
    result = map_reaction_smiles("CCO>>CC=O")
    assert result["status"] == "inferred"
    assert len(result["mappings"]) == 2
    assert result["unresolved_product_carbons"] == []


def test_carboxylation_like_new_carbon_is_unresolved():
    result = map_reaction_smiles("CCO>>CC(=O)O")
    assert result["status"] == "unresolved"
    assert len(result["mappings"]) == result["product_carbon_atom_count"]
    assert len(result["unresolved_product_carbons"]) >= 1


def test_explicit_co2_carbon_source_is_mapped_to_carboxyl_carbon():
    result = map_reaction_smiles("CC(=O)O.O=C=O>>CC(=O)C(=O)O")
    assert any(mapping["method"].startswith("rdkit-co2-carbon-source") for mapping in result["mappings"])


def test_unique_mcs_resolves_a_carbon_rearrangement():
    result = map_reaction_smiles(
        "O=C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-]>>"
        "O=C(CO)[C@H](O)[C@H](O)COP(=O)([O-])[O-]"
    )
    assert result["status"] == "inferred"
    assert all(mapping["status"] == "inferred" for mapping in result["mappings"])
    assert any(mapping["method"] == "rdkit-mcs-carbon-conservation" for mapping in result["mappings"])
