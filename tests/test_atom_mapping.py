from cannabis_carbon.atom_mapping import apply_reaction_specific_candidate_mapping, map_reaction_smiles, map_identity_pair_smiles


def test_oxidation_preserves_all_product_carbons():
    result = map_reaction_smiles("CCO>>CC=O")
    assert result["status"] == "inferred"
    assert len(result["mappings"]) == 2
    assert result["unresolved_product_carbons"] == []


def test_heteroatom_addition_with_conserved_carbon_skeleton_is_inferred():
    result = map_reaction_smiles("CCO>>CC(=O)O")
    assert result["status"] == "inferred"
    assert {row["method"] for row in result["mappings"]} == {"rdkit-carbon-skeleton-mcs"}


def test_unique_carbon_skeleton_maps_retained_product_when_carbon_is_lost():
    result = map_reaction_smiles("CCOC>>CC(=O)O")
    assert result["status"] == "inferred"
    assert len(result["mappings"]) == 2
    assert all(row["method"] == "rdkit-carbon-skeleton-mcs" for row in result["mappings"])


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
    assert any(mapping["method"] in ("rdkit-mcs-carbon-conservation-relaxed-bond", "rdkit-full-carbon-mcs-conservation") for mapping in result["mappings"])


def test_unique_product_substructure_maps_retained_carbons():
    result = map_reaction_smiles("CC(=O)O>>CC")
    assert result["status"] == "inferred"
    assert len(result["mappings"]) == 2
    assert all(mapping["method"] == "rdkit-unique-product-substructure" for mapping in result["mappings"])


def test_decarboxylation_maps_released_carbon_to_co2():
    result = map_reaction_smiles("CC(=O)O>>C=O.O=C=O")
    assert result["status"] == "inferred"
    assert len(result["mappings"]) == 2
    assert any(mapping["method"] == "rdkit-decarboxylation-released-carbon" for mapping in result["mappings"])


def test_unresolved_product_co2_gets_only_a_candidate_carboxyl_source():
    reaction = "CC(=O)O>>C.O=C=O"
    result = apply_reaction_specific_candidate_mapping("rhea:test", reaction, map_reaction_smiles(reaction))
    co2 = [row for row in result["mappings"] if row["product_index"] == 1][0]
    assert co2["status"] == "candidate"
    assert co2["method"] == "rdkit-co2-product-decarboxylation-candidate"
    assert co2["alternatives"] == [{"reactant_index": 0, "reactant_atom": 1}]
    assert any(row["status"] == "unresolved" for row in result["mappings"] if row["product_index"] == 0)


def test_full_carbon_mcs_maps_oxidation_with_oxygen_only_cofactor():
    result = map_reaction_smiles("CCCCCC1=CC(O)=C(CC=C(C)CCC=C(C)C)C(O)=C1C(O)=O.O=O>>CCCCCC1=C(C(O)=O)C(O)=C([C@@H]2C=C(C)CC[C@H]2C(C)=C)C(O)=C1.OO")
    assert result["status"] == "inferred"
    assert len(result["mappings"]) == 22
    assert all(row["method"] == "rdkit-full-carbon-mcs-conservation" for row in result["mappings"])


def test_equivalent_redox_carbon_skeleton_is_candidate_not_unresolved():
    result = map_reaction_smiles("CC(=O)O.CC>>CC(O)O.CC")
    assert result["status"] == "candidate"
    assert result["unresolved_product_carbons"] == []
    assert any(item["method"] == "rdkit-equivalent-carbon-skeleton-candidate" for item in result["mappings"])


def test_identity_pair_mapper_handles_heteroatom_change_without_cofactor_noise():
    result = map_identity_pair_smiles("CCO", "CC=O")
    assert result["status"] == "inferred"
    assert len(result["mappings"]) == 2
    assert result["unresolved_product_carbons"] == []


def test_identity_pair_mapper_keeps_carbon_gain_unresolved():
    result = map_identity_pair_smiles("CC", "CCC")
    assert result["status"] == "unresolved"
    assert result["reason"] == "identity-pair-carbon-count-delta"
    assert len(result["mappings"]) == 3
    assert len(result["unresolved_product_carbons"]) == 1
