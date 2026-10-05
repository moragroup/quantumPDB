import copy

import numpy as np
import pytest

from qp.reaction.cyp_adapter import adapt
from qp.reaction.cyp_preparation import oxo_position
from qp.reaction.schema import HypothesisError


def test_oxo_is_distal_to_axial_sulfur():
    fe = np.array([2., 3., 4.])
    sulfur = fe + [0., 0., -2.2]
    oxo = oxo_position(fe, sulfur)
    assert np.linalg.norm(oxo - fe) == pytest.approx(1.65)
    assert np.dot(oxo - fe, sulfur - fe) < 0
    with pytest.raises(ValueError):
        oxo_position(fe, fe)


def test_metal_hat_adapter_preserves_audited_graph_and_rejects_other_edits():
    # Minimal bookkeeping fixture, not a physically prepared heme model.
    template = {
        "schema_version": "1.0", "reaction_id": "fixture", "hypothesis_id": "hat",
        "hydrogen_policy": "explicit", "interpretation": {"kind": "elementary", "step_id": "HAT"},
        "provenance": {"source": "test fixture"}, "participants": [{"id": "test"}],
        "atoms": [{"map_id": i, "element": e, "participant": "test"}
                  for i, e in enumerate(("C", "H", "O", "Fe"), 1)],
        "electronic_state": {"charge": 0, "multiplicity": 2}, "endpoints": {}}
    for side in ("reactant", "product"):
        template["endpoints"][side] = {
            "atom_states": [{"map_id": i, "formal_charge": 0,
                             "radical_electrons": int(i == (3 if side == "reactant" else 1))}
                            for i in range(1, 5)],
            "bonds": [{"atoms": [1, 2] if side == "reactant" else [2, 3], "order": 1},
                      {"atoms": [3, 4], "order": 2 if side == "reactant" else 1}]}
    left = np.zeros((4, 4))
    left[0, 1] = left[1, 0] = 1
    left[2, 2] = 1
    right = np.zeros((4, 4))
    right[1, 2] = right[2, 1] = 1
    right[0, 0] = 1
    snapshot = {"elements": ["C", "H", "O", "Fe"], "map_ids": [1, 2, 3, 4],
                "reactant_be": left.tolist(), "candidates": [{"product_be": right.tolist()}]}
    spec = {"oxo_map": 3, "substrate_carbon_maps": [1], "metal_state_audit": "Fixture only",
            "allowed_scaffold_be_edits": [{"atoms": [1, 1], "delta": 1},
                                          {"atoms": [3, 3], "delta": -1}]}
    result = adapt(template, snapshot, spec)
    assert result["endpoints"] == template["endpoints"]
    bad = copy.deepcopy(snapshot)
    bad["candidates"][0]["product_be"][3][3] += 1
    bad["candidates"][0]["product_be"][0][0] -= 1
    with pytest.raises(HypothesisError, match="Unaudited BE edit"):
        adapt(template, bad, spec)
    bad = copy.deepcopy(spec)
    bad["substrate_carbon_maps"] = []
    with pytest.raises(HypothesisError, match="outside"):
        adapt(template, snapshot, bad)
