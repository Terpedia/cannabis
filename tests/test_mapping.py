from rdkit import Chem

from cannabis_carbon.mapping import map_conserved_carbons, unresolved_product_carbons


def test_conserved_carbon_mapping():
    reactant = Chem.MolFromSmiles("CCO")
    product = Chem.MolFromSmiles("CC=O")
    mapping = map_conserved_carbons(reactant, product)
    assert len(mapping) == 2
    assert unresolved_product_carbons(product, mapping) == []


def test_new_carbon_is_unresolved():
    reactant = Chem.MolFromSmiles("CCO")
    product = Chem.MolFromSmiles("CCCO")
    mapping = map_conserved_carbons(reactant, product)
    assert len(unresolved_product_carbons(product, mapping)) >= 1
