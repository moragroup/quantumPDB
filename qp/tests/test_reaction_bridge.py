"""Chemical/indexing/convergence regressions for the reaction bridge."""

import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from qp.reaction.jobs import MODES, orca_input
from qp.reaction.mapping import build_cluster_model, graph_mappings, validate_mapping, write_xyz
from qp.reaction.results import (endpoint_matches, mode_overlap, normal_modes,
                                 parse_orca, validate_step)
from qp.reaction.schema import HypothesisError, reaction_edits, validate, write_json
from qp.reaction.workflow import _edit_scans, collect, compare, prepare, stage_jobs


def exchange():
    return {"schema_version": "1.0", "reaction_id": "hydrogen_exchange",
            "hypothesis_id": "fixture", "hydrogen_policy": "explicit",
            "interpretation": {"kind": "elementary", "step_id": "1"},
            "provenance": {"source": "synthetic test, not a model prediction"},
            "participants": [{"id": "hydrogen", "role": "substrate"}],
            "atoms": [{"map_id": i, "element": "H", "participant": "hydrogen"} for i in (11, 23, 47)],
            "endpoints": {s: {"atom_states": [{"map_id": i, "formal_charge": 0,
                                               "radical_electrons": int(i == radical)} for i in (11, 23, 47)],
                                "bonds": [{"atoms": pair, "order": 1.0}]}
                          for s, pair, radical in (("reactant", [11, 23], 47), ("product", [23, 47], 11))},
            "electronic_state": {"charge": 0, "multiplicity": 2}}


def model():
    return {"elements": ["H", "H", "H"], "charge": 0, "multiplicity": 2,
            "reactant_xyz": "r.xyz", "product_xyz": "p.xyz", "ts_xyz": "ts.xyz",
            "frozen_indices": [], "atom_mapping": [
                {"map_id": a, "qm_index": i, "structure_key": {"participant": "hydrogen", "atom": a}}
                for i, a in enumerate((11, 23, 47))]}


def test_schema_conservation_and_net_rejection():
    record = exchange()
    assert validate(record) is record
    assert reaction_edits(record)[0]["atoms"] == [11, 23]
    altered = copy.deepcopy(record)
    altered["interpretation"]["kind"] = "net"
    with pytest.raises(HypothesisError, match="decomposition"):
        validate(altered)
    assert validate(altered, require_elementary=False)
    altered = copy.deepcopy(record)
    altered["endpoints"]["product"]["atom_states"].pop()
    with pytest.raises(HypothesisError, match="same unique"):
        validate(altered)
    altered = copy.deepcopy(record)
    altered["electronic_state"]["multiplicity"] = 1
    with pytest.raises(HypothesisError, match="parity"):
        validate(altered)


def test_mapping_symmetry_and_stereo():
    record = exchange()
    graph = {"atoms": [{"id": str(i), "element": "H", "formal_charge": 0,
                        "radical_electrons": int(i == 3)} for i in (1, 2, 3)],
             "bonds": [{"atoms": ["1", "2"], "order": 1}]}
    matches = graph_mappings(record, graph)
    assert len(matches) == 2
    assert all(m[47] == "3" for m in matches)
    graph["atoms"][0]["stereo"] = "R"
    with pytest.raises(HypothesisError, match="stereochemistry"):
        graph_mappings(record, graph)
    graph["atoms"][0].pop("stereo")
    graph["bonds"].append({"atoms": ["2", "1"], "order": 1})
    with pytest.raises(HypothesisError, match="prepared graph bond"):
        graph_mappings(record, graph)


def test_reactive_boundaries_and_element_order():
    record, prepared = exchange(), model()
    validate_mapping(record, prepared)
    prepared["frozen_indices"] = [1]
    with pytest.raises(HypothesisError, match="cannot be frozen"):
        validate_mapping(record, prepared)
    prepared["frozen_indices"] = []
    prepared["atom_mapping"][0]["is_cap"] = True
    with pytest.raises(HypothesisError, match="caps"):
        validate_mapping(record, prepared)
    prepared["atom_mapping"][0]["is_cap"] = False
    prepared["elements"][0] = "C"
    with pytest.raises(HypothesisError, match="element"):
        validate_mapping(record, prepared)
    prepared = model()
    prepared["charge"] = 0.0
    with pytest.raises(HypothesisError, match="integer cluster/environment"):
        validate_mapping(record, prepared)


def test_all_job_modes_move_heavy_atoms_and_index_scans():
    prepared = model()
    for mode in MODES:
        text = orca_input(mode, prepared, {}, scans=[{"indices": [0, 1], "start": .75, "end": 2., "steps": 8}])
        assert "optimizehydrogens" not in text.lower()
        assert "* xyzfile 0 2" in text
    assert "B 0 1 =" in orca_input("constrained_scan", prepared, {}, scans=[{"indices": [0, 1], "start": .75, "end": 2., "steps": 8}])
    prepared["frozen_indices"] = [2]
    assert "{ C 2 C }" in orca_input("endpoint_opt", prepared, {})
    assert "Partial_Hess" in orca_input("frequency", prepared, {})
    assert "Partial_Hess { 2 } end" in orca_input("frequency", prepared, {})
    assert "NumHessTransInvar false" in orca_input("frequency", prepared, {})
    with pytest.raises(HypothesisError, match="Frozen-boundary"):
        orca_input("irc", prepared, {})


def test_parser_and_normal_modes(tmp_path):
    output = tmp_path / "job.out"
    output.write_text("Program Version 6.1.1\nFINAL SINGLE POINT ENERGY -1.60\n"
                      "THE OPTIMIZATION HAS CONVERGED\n  0: -900.0 cm**-1\n 1: 0.0 cm**-1\n"
                      "TOTAL RUN TIME: 0 days 0 hours 0 minutes 2 seconds 100 msec\nORCA TERMINATED NORMALLY\n")
    parsed = parse_orca(output, "frequency")
    assert parsed["converged"]
    assert parsed["frequencies_cm_1"][0]["frequency"] == -900
    assert parsed["engine_wall_seconds"] == pytest.approx(2.1)
    output.write_text(output.read_text() + "Number of displaced atoms ... 2\nList of displaced atoms ... 1 3\n\nNumber of displacements ... 12\n")
    assert parse_orca(output, "frequency")["displaced_atoms_0based"] == [0, 2]
    output.write_text(output.read_text() + "SCF NOT CONVERGED\n")
    assert not parse_orca(output, "endpoint_opt")["converged"]
    hessian = tmp_path / "job.hess"
    hessian.write_text("$normal_modes\n3 3\n 0 1\n0 1.0 0.0\n1 0.0 1.0\n2 0.0 0.0\n 2\n0 0.0\n1 0.0\n2 1.0\n# ORCA atom section\n$atoms\n")
    assert np.allclose(normal_modes(hessian), np.eye(3))


def test_chemical_validation_requires_irc_and_reaction_mode():
    record, prepared = exchange(), model()
    xyz = np.asarray([[0., 0, 0], [.74, 0, 0], [3.7, 0, 0]])
    assert endpoint_matches(record, prepared, xyz, "reactant")
    assert not endpoint_matches(record, prepared, xyz, "product")
    assert mode_overlap(record, prepared, xyz, [[-1., 0, 0], [1., 0, 0], [0., 0, 0]]) > .99
    evidence = {s: {"converged": True, "energy_hartree": e} for s, e in
                (("reactant", -1.7), ("product", -1.7), ("ts", -1.6))}
    evidence["reactant"]["matches"], evidence["product"]["matches"] = ["reactant"], ["product"]
    for side in ("reactant", "product", "ts"):
        evidence[side + "_frequency"] = {"converged": True, "frequencies_cm_1": [
            {"index": 0, "frequency": -900. if side == "ts" else 100.}], "reaction_mode_overlap": .9}
    assert validate_step(record, prepared, evidence)["status"] == "not_converged"
    evidence.update(irc={"converged": True}, connected_a={"converged": True, "matches": ["reactant"]},
                    connected_b={"converged": True, "matches": ["product"]})
    result = validate_step(record, prepared, evidence)
    assert result["status"] == "validated_step"
    assert result["forward_barrier"] == pytest.approx(62.7509474)
    evidence["connected_b"]["matches"] = ["reactant"]
    assert validate_step(record, prepared, evidence)["status"] == "wrong_connectivity"


def test_workflow_dependencies_and_incompatible_restart(tmp_path):
    record, prepared = exchange(), model()
    write_json(tmp_path / "hypothesis.json", record)
    write_json(tmp_path / "model.json", prepared)
    for name, xyz in (("r", [[0, 0, 0], [.74, 0, 0], [3.74, 0, 0]]),
                      ("p", [[0, 0, 0], [3., 0, 0], [3.74, 0, 0]]),
                      ("ts", [[0, 0, 0], [.93, 0, 0], [1.86, 0, 0]])):
        write_xyz(tmp_path / (name + ".xyz"), prepared["elements"], xyz)
    root = prepare(tmp_path / "hypothesis.json", tmp_path / "model.json", tmp_path / "bundle", strategy="ts_guess")
    state = json.loads((root / "state.json").read_text())
    assert {j["job"] for j in state["jobs"]} == {"reactant", "product", "ts"}
    assert collect(root)["status"] == "not_converged"
    settings = root / "settings.json"
    settings.write_text('{"method":"HF"}')
    with pytest.raises(HypothesisError, match="changed"):
        stage_jobs(root)


def test_fluxion_snapshot_keeps_maps_and_rejects_electron_loss(tmp_path):
    path = Path(__file__).resolve().parents[3] / "Fluxion/scripts/export_qm_hypotheses.py"
    if not path.exists():
        pytest.skip("Sibling Fluxion repository unavailable")
    spec = importlib.util.spec_from_file_location("exporter", path)
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    snapshot = {"elements": ["H"] * 3, "map_ids": [11, 23, 47],
                "reactant_be": [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
                "candidates": [{"product_be": [[1, 0, 0], [0, 0, 1], [0, 1, 0]], "count": 5},
                               {"product_be": [[0, 0, 0], [0, 0, 1], [0, 1, 0]], "count": 1}]}
    metadata = {"reaction_id": "fixture", "interpretation": {"kind": "elementary", "step_id": "1"},
                "electronic_state": {"charge": 0, "multiplicity": 2}}
    report = exporter.export_snapshot(snapshot, metadata, tmp_path)
    assert len(report["exported"]) == 1 and len(report["rejected"]) == 1
    record = json.loads(Path(report["exported"][0]).read_text())
    validate(record)
    assert [a["map_id"] for a in record["atoms"]] == [11, 23, 47]
    assert any(t["hydrogen"] == 23 and t["donors"] == [11] and t["acceptors"] == [47]
               for t in record["proton_transfers"])


def test_cluster_mapping_uses_topology_anchors_and_sphere_order(tmp_path):
    from Bio.PDB import Atom, Chain, Model, PDBIO, Residue, Structure
    structure = Structure.Structure("cluster")
    parent = Model.Model(0)
    chain = Chain.Chain("B")
    structure.add(parent)
    parent.add(chain)
    residue = Residue.Residue(("H_LIG", 9, " "), "LIG", " ")
    chain.add(residue)
    for i, name in enumerate(("HA", "HB", "HC")):
        residue.add(Atom.Atom(name, np.array([float(i), 0., 0.]), 0, 1, " ", name, 17 + i, "H"))
    sphere = tmp_path / "0.pdb"
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(sphere))
    specification = {"charge": 0, "multiplicity": 2, "frozen_indices": [], "graphs": [{
        "participant": "hydrogen", "atoms": [{"id": name, "qm_index": i, "element": "H",
            "formal_charge": 0, "radical_electrons": int(i == 2)} for i, name in enumerate(("HA", "HB", "HC"))],
        "bonds": [{"atoms": ["HA", "HB"], "order": 1}]}]}
    with pytest.raises(HypothesisError, match="Ambiguous"):
        build_cluster_model(exchange(), [sphere], specification, tmp_path / "ambiguous")
    specification["graphs"][0]["anchors"] = {"11": "HB"}
    prepared = build_cluster_model(exchange(), [sphere], specification, tmp_path / "selected")
    assert {a["map_id"]: a["qm_index"] for a in prepared["atom_mapping"]} == {11: 1, 23: 0, 47: 2}
    assert prepared["atom_mapping"][0]["structure_key"]["chain"] == "B"
    assert "qm_index_1based" in (tmp_path / "selected/atom_mapping.tsv").read_text()


def test_guided_scan_indices_follow_persistent_maps(tmp_path):
    prepared = model()
    write_xyz(tmp_path / "r.xyz", prepared["elements"], [[0, 0, 0], [.75, 0, 0], [3.75, 0, 0]])
    write_xyz(tmp_path / "p.xyz", prepared["elements"], [[0, 0, 0], [3, 0, 0], [3.75, 0, 0]])
    scans = _edit_scans(exchange(), prepared, tmp_path / "r.xyz", tmp_path / "p.xyz", {})
    assert [s["indices"] for s in scans] == [[0, 1], [1, 2]]
    assert scans[0]["start"] == .75 and scans[1]["end"] == .75
    assert "Simul_Scan true" in orca_input("constrained_scan", prepared, {}, scans=scans)


def test_parser_rejects_partial_neb_refinement_and_irc_summary(tmp_path):
    output = tmp_path / "job.out"
    output.write_text("THE NEB OPTIMIZATION HAS CONVERGED\nORCA TERMINATED NORMALLY\n")
    assert not parse_orca(output, "path_search")["converged"]
    output.write_text("IRC PATH SUMMARY\nORCA TERMINATED NORMALLY\n")
    assert not parse_orca(output, "irc")["converged"]


def test_geometry_generation_preserves_environment_and_product_graph():
    pytest.importorskip("rdkit")
    from qp.reaction.geometry import generate_product
    prepared = model()
    prepared["elements"].append("He")
    prepared["frozen_indices"] = [3]
    xyz = np.array([[0., 0., 0.], [.75, 0., 0.], [3.75, .05, 0.], [10., 3., 2.]])
    product = generate_product(exchange(), prepared, xyz)
    assert np.array_equal(product[3], xyz[3])
    assert endpoint_matches(exchange(), prepared, product, "product")
    record = exchange()
    record["provenance"]["stereochemistry_audit_required"] = True
    with pytest.raises(HypothesisError, match="audited mapped product"):
        generate_product(record, prepared, xyz)


def test_restart_detects_job_tampering_and_comparison_checks_geometries(tmp_path):
    record, prepared = exchange(), model()
    write_json(tmp_path / "hypothesis.json", record)
    write_json(tmp_path / "model.json", prepared)
    for name, xyz in (("r", [[0, 0, 0], [.74, 0, 0], [3.74, 0, 0]]),
                      ("p", [[0, 0, 0], [3., 0, 0], [3.74, 0, 0]]),
                      ("ts", [[0, 0, 0], [.93, 0, 0], [1.86, 0, 0]])):
        write_xyz(tmp_path / (name + ".xyz"), prepared["elements"], xyz)
    first = prepare(tmp_path / "hypothesis.json", tmp_path / "model.json", tmp_path / "first")
    write_xyz(tmp_path / "r.xyz", prepared["elements"], [[0, 0, 0], [.76, 0, 0], [3.74, 0, 0]])
    second = prepare(tmp_path / "hypothesis.json", tmp_path / "model.json", tmp_path / "second")
    assert not compare([first, second], tmp_path / "comparison.json")["matched_settings"]
    script = first / "candidates/001/reactant/run.sh"
    script.write_text(script.read_text() + "# edited launch protocol\n")
    with pytest.raises(HypothesisError, match="Staged job inputs changed"):
        stage_jobs(first)


def test_heavy_capping_receives_chain(monkeypatch):
    from Bio.PDB import Atom, Chain, Model, Residue, Structure
    from qp.cluster import spheres
    chain = Chain.Chain("A")
    residues = []
    for i in range(1, 4):
        residue = Residue.Residue((" ", i, " "), "ALA", " ")
        for name, element in (("N", "N"), ("CA", "C"), ("C", "C")):
            residue.add(Atom.Atom(name, np.asarray([float(i), 0., 0.]), 0, 1, " ", name, i, element))
        chain.add(residue)
        residues.append(residue)
    model_object = Model.Model(0)
    model_object.add(chain)
    structure = Structure.Structure("test")
    structure.add(model_object)
    calls = []
    monkeypatch.setattr(spheres, "build_heavy", lambda *args: calls.append(args) or args[2])
    spheres.cap_chains(model_object, {residues[1]}, 2)
    assert len(calls) == 2
    assert all(args[0] is chain and args[1] is residues[1] for args in calls)
