"""Numerical and mapping checks for the offline apo/holo analysis."""

import gzip

import numpy as np

from qp.analyze.catalytic_movement import atom_rmsd, fit_coordinates, load_structure, sequence_mapping


def test_rigid_transform_is_removed_but_local_movement_remains():
    fixed = np.array([[0., 0., 0.], [3., 0., 0.], [0., 4., 0.], [0., 0., 5.]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    offset = np.array([10., -20., 3.])
    mobile = fixed @ rotation + offset
    r, t, rmsd = fit_coordinates(mobile, fixed)
    assert rmsd < 1e-10
    assert np.linalg.det(r) > 0
    apo = {"CA": np.array([1., 2., 3.])}
    holo = {"CA": (apo["CA"] + np.array([2., 0., 0.])) @ rotation + offset}
    assert np.isclose(atom_rmsd(apo, holo, {"CA"}, r, t, "GLY"), 2.)


def test_equivalent_carboxylate_labels_do_not_count_as_movement():
    apo = {"CB": np.array([0., 0., 0.]), "CG": np.array([1., 0., 0.]),
           "OD1": np.array([2., 1., 0.]), "OD2": np.array([2., -1., 0.])}
    holo = dict(apo, OD1=apo["OD2"], OD2=apo["OD1"])
    assert atom_rmsd(apo, holo, set(apo), np.eye(3), np.zeros(3), "ASP") == 0


def test_aromatic_swaps_are_coupled():
    apo = {"CD1": np.array([0., 1., 0.]), "CD2": np.array([0., -1., 0.]),
           "CE1": np.array([1., 1., 0.]), "CE2": np.array([1., -1., 0.])}
    holo = dict(apo, CD1=apo["CD2"], CD2=apo["CD1"])
    assert atom_rmsd(apo, holo, set(apo), np.eye(3), np.zeros(3), "PHE") > 0


def test_reference_mapping_handles_truncation_and_tags():
    reference = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
    sequence = reference[8:32] + "HHHHHH"
    mapping, identity, coverage = sequence_mapping(reference, sequence)
    assert mapping[9] == 1
    assert mapping[32] == 24
    assert identity > 0.95
    assert coverage >= 24 / len(reference)


def test_mmcif_uses_label_positions_entity_first_model_and_occupancy(tmp_path):
    directory = tmp_path / "ab"
    directory.mkdir()
    fields = ["label_asym_id", "label_seq_id", "label_comp_id", "label_atom_id", "type_symbol",
              "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "label_alt_id", "pdbx_PDB_model_num",
              "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code"]
    text = "data_test\nloop_\n_entity_poly.entity_id\n_entity_poly.pdbx_seq_one_letter_code_can\n1 AC\n2 AC\n"
    text += "loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2\nL 3\nloop_\n"
    text += "\n".join("_atom_site." + f for f in fields) + "\n"
    text += "A 2 CYS CA C 1 2 3 0.3 A 1 X 502 A\n"
    text += "A 2 CYS CA C 4 5 6 0.7 B 1 X 502 A\n"
    text += "A 2 CYS CA C 9 9 9 1 . 2 X 502 A\n"
    text += "B 2 CYS CA C 8 8 8 1 . 1 Y 99 ?\n"
    text += "L . ZN ZN ZN 5 5 5 1 . 1 Z 1 ?\n"
    with gzip.open(directory / "1abc.cif.gz", "wt") as handle:
        handle.write(text)
    chains, ligands, _ = load_structure(tmp_path, "1ABC", {"1"})
    assert set(chains) == {"A"}
    residue = chains["A"]["residues"][2]
    assert residue["author_pos"] == "502A"
    assert np.allclose(residue["atoms"]["CA"], [4, 5, 6])
    assert [code for code, xyz in ligands] == ["ZN"]
