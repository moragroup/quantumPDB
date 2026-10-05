from qp.analyze.cofactor_dataset import classify_contacts, ligand_roles


def test_reaction_mapping_preserves_stereochemistry_and_excludes_products():
    compounds = [{"chebi_id": "1", "type": "reactant"},
                 {"chebi_id": "2", "type": "product"}]
    chebi = {"1": {"default_structure": {"standard_inchi_key": "AAA-BBB-M"}},
             "2": {"default_structure": {"standard_inchi_key": "CCC-DDD-N"}}}
    components = {"SUB": {"inchi_key": "AAA-BBB-N", "weight": 200},
                  "ISO": {"inchi_key": "AAA-EEE-N", "weight": 200},
                  "PRO": {"inchi_key": "CCC-DDD-N", "weight": 200},
                  "FAD": {"cross_links": [{"resource": "ChEBI", "resource_id": "CHEBI:1"}], "weight": 700},
                  "UNK": {"cross_links": None, "weight": None}}
    assert ligand_roles(compounds, components, {"FAD"}, chebi) == ({"SUB"}, {"PRO"})


def test_cofactor_only_rejects_inhibitors_and_retains_cofactor_sets():
    args = ({"SUB"}, {"FAD", "ZN"}, {"SO4"})
    assert classify_contacts({"FAD", "ZN", "SO4"}, *args)[:3] == (
        "substrate_free", {"FAD", "ZN"}, set())
    assert classify_contacts({"FAD", "SUB"}, *args)[:3] == (
        "substrate_bound", {"FAD"}, {"SUB"})
    assert classify_contacts({"FAD", "INH"}, *args)[0] == "excluded"
    assert classify_contacts({"SUB"}, *args)[0] == "excluded"
