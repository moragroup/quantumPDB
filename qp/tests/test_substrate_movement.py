"""Substrate-state checks independent of the global apo/holo labels."""

import numpy as np

from qp.analyze.substrate_movement import classify_chain


def test_substrate_state_is_chain_specific_and_requires_matching_cofactor():
    reference = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
    chain = {"sequence": reference, "residues": {1: {"atoms": {"CA": np.zeros(3)}}}}
    ligands = [("FAD", np.array([1., 0., 0.])), ("TRP", np.array([20., 0., 0.]))]
    classify = lambda c, l: classify_chain(c, l, reference, "TRP", "FAD", {"CTE"}, 4.)
    assert classify(chain, ligands)[0] == "substrate_free"
    assert classify(chain, ligands + [("TRP", np.array([2., 0., 0.]))])[0] == "substrate_bound"
    assert classify(chain, [("FDA", np.ones(3))])[0] == "excluded"
    assert classify(chain, ligands + [("CTE", np.ones(3))])[0] == "excluded"
    mutant = dict(chain, sequence=reference[:20] + "A" + reference[21:])
    assert classify(mutant, ligands)[0] == "excluded"
