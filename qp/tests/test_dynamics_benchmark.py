"""Check coordinate mapping and the benchmark's fixed-atom conclusion."""

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, PDBIO, Residue, Structure

from qp.analyze.dynamics_benchmark import measure_output, read_xyz


def test_relaxation_mapping_and_rigid_alignment(tmp_path):
    structure = Structure.Structure("test")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    structure.add(model)
    model.add(chain)
    coordinates = []
    for number, xyz in enumerate([(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2)], 1):
        residue = Residue.Residue((" ", number, " "), "ALA", " ")
        chain.add(residue)
        for name, coordinate in [("CA", xyz)] + ([("CB", (1, 1, 2))] if number == 4 else []):
            residue.add(Atom.Atom(name, np.array(coordinate, dtype=float), 0, 1, " ", name, len(coordinates) + 1, element="C"))
            coordinates.append(np.array(coordinate, dtype=float))
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(tmp_path / "0.pdb"))
    final = np.array(coordinates) + [10, -3, 2]
    final[-1] += [1, 0, 0]
    xyz = tmp_path / "final.xyz"
    xyz.write_text(str(len(final)) + "\n\n" + "\n".join("C {} {} {}".format(*p) for p in final) + "\n")
    result = measure_output(tmp_path, xyz, "A", "4", "ALA")
    assert result["status"] == "measured"
    assert result["simulated_sidechain_relaxation_A"] == pytest.approx(1)
    assert result["cluster_fit_rmsd_A"] < 1e-12
    assert measure_output(tmp_path, xyz, "A", "5", "ALA")["status"] == "residue_absent"
    xyz.write_text(xyz.read_text().replace("C ", "N ", 1))
    with pytest.raises(ValueError, match="order"):
        measure_output(tmp_path, xyz, "A", "4", "ALA")


def test_xyz_rejects_multiple_frames(tmp_path):
    path = tmp_path / "trajectory.xyz"
    path.write_text("1\n\nC 0 0 0\n1\n\nC 1 0 0\n")
    with pytest.raises(ValueError, match="one final"):
        read_xyz(path)
