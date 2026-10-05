"""Reconcile a CYP BE HAT prediction against an audited physical template.

Metal electron bookkeeping is never translated with organic valence rules.
The caller supplies an audited physical HAT endpoint template and explicitly
declares every permitted scaffold BE edit. Other prediction edits are rejected.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from qp.reaction.schema import (HypothesisError, reaction_edits, sha256,
                                validate, write_json)


def adapt(template, snapshot, specification, candidate_index=0):
    validate(template)
    maps = snapshot["map_ids"]
    elements = snapshot["elements"]
    if len(maps) != len(set(maps)) or len(maps) != len(elements):
        raise HypothesisError("Invalid prediction inventory")
    inventory = {a["map_id"]: a["element"] for a in template["atoms"]}
    if any(inventory.get(i) != e for i, e in zip(maps, elements)):
        raise HypothesisError("Prediction maps must explicitly match audited physical template")
    x0 = np.asarray(snapshot["reactant_be"], float)
    x1 = np.asarray(snapshot["candidates"][candidate_index]["product_be"], float)
    for x in (x0, x1):
        if (x.shape != (len(maps), len(maps)) or not np.isfinite(x).all()
                or not np.allclose(x, x.T, atol=1e-8, rtol=0)):
            raise HypothesisError("Invalid BE matrix")
    if not np.isclose(x0.sum(), x1.sum(), atol=1e-6, rtol=0):
        raise HypothesisError("Prediction does not conserve BE electrons")
    delta = x1 - x0
    oxo = specification["oxo_map"]
    if inventory.get(oxo) != "O":
        raise HypothesisError("Oxo map must identify oxygen")
    broken, formed = [], []
    observed = []
    for i in range(len(maps)):
        for j in range(i, len(maps)):
            if abs(delta[i, j]) < 1e-6:
                continue
            edit = {"atoms": [maps[i], maps[j]], "delta": float(delta[i, j])}
            observed.append(edit)
            pair = {elements[i], elements[j]}
            if i != j and pair == {"C", "H"} and delta[i, j] < 0:
                broken.append((maps[i] if elements[i] == "C" else maps[j],
                               maps[j] if elements[j] == "H" else maps[i]))
            if i != j and oxo in (maps[i], maps[j]) and pair == {"O", "H"} and delta[i, j] > 0:
                formed.append(maps[j] if elements[j] == "H" else maps[i])
    if len(broken) != 1 or formed != [broken[0][1]]:
        raise HypothesisError("Prediction must specify exactly one matched C-H / oxo-H transfer")
    carbon, hydrogen = broken[0]
    if carbon not in specification["substrate_carbon_maps"]:
        raise HypothesisError("Transfer donor is outside the declared substrate")
    allowed = {tuple(sorted(e["atoms"])): float(e["delta"])
               for e in specification.get("allowed_scaffold_be_edits", [])}
    for edit in observed:
        a, b = edit["atoms"]
        if {a, b} in ({carbon, hydrogen}, {oxo, hydrogen}):
            continue
        # Atomic bookkeeping at the radical donor and transferring H must also
        # be specified; the adapter does not silently reinterpret these edits.
        if not np.isclose(allowed.get(tuple(sorted((a, b))), float("inf")), edit["delta"], atol=1e-6):
            raise HypothesisError("Unaudited BE edit: " + str(edit))
    record = copy.deepcopy(template)
    template_edits = reaction_edits(record)
    template_broken = [e for e in template_edits if e["product_order"] == 0
                       and inventory[e["atoms"][0]] != "Fe" and inventory[e["atoms"][1]] != "Fe"]
    template_formed = [e for e in template_edits if e["reactant_order"] == 0]
    if (len(template_broken) != 1 or len(template_formed) != 1
            or set(template_broken[0]["atoms"]) != {carbon, hydrogen}
            or set(template_formed[0]["atoms"]) != {oxo, hydrogen}):
        raise HypothesisError("An independently audited site/H-specific physical template is required")
    record["provenance"].update({
        "cyp_conversion_policy": "Audited physical endpoint template; BE edits reconciled without Fe valence conversion",
        "raw_prediction_provenance": snapshot.get("provenance", {}),
        "candidate_index": candidate_index, "observed_be_edits": observed,
        "metal_state_audit": specification["metal_state_audit"],
        "physical_template_is_independently_required": True})
    return validate(record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True)
    parser.add_argument("--snapshot", required=True)
    parser.add_argument("--specification", required=True)
    parser.add_argument("--candidate-index", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    record = adapt(*(json.loads(Path(p).read_text()) for p in
                     (args.template, args.snapshot, args.specification)), args.candidate_index)
    record["provenance"]["adapter_inputs_sha256"] = {
        str(Path(p).resolve()): sha256(p) for p in (args.template, args.snapshot, args.specification)}
    write_json(args.output, record)


if __name__ == "__main__":
    main()
