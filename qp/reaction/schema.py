"""Validate the JSON contract shared with reaction-prediction programs.

No prediction/training packages are imported here. Atom maps are persistent
chemical identities; QM indices are assigned only during structure mapping.
"""

import hashlib
import json
import math
from pathlib import Path

from qp.reaction import SCHEMA_VERSION


class HypothesisError(ValueError):
    """A proposal is unsuitable for a physical reaction search."""


ATOMIC_NUMBERS = dict(zip(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Mo I W".split(),
    list(range(1, 37)) + [42, 53, 74],
))
VALENCE_ELECTRONS = {"H": 1, "B": 3, "C": 4, "N": 5, "O": 6,
                     "F": 7, "Si": 4, "P": 5, "S": 6, "Cl": 7,
                     "Br": 7, "I": 7}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def content_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                     separators=(",", ":")).encode()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def bonds(endpoint):
    result = {}
    for bond in endpoint["bonds"]:
        a, b = bond["atoms"]
        key = tuple(sorted((a, b)))
        if a == b or key in result:
            raise HypothesisError("Self/duplicate bond: {}".format(key))
        order = bond["order"]
        if not isinstance(order, (int, float)) or not math.isfinite(order) or order <= 0:
            raise HypothesisError("Invalid bond order")
        result[key] = float(order)
    return result


def reaction_edits(record):
    left, right = (bonds(record["endpoints"][s]) for s in ("reactant", "product"))
    return [{"atoms": list(key), "reactant_order": left.get(key, 0.0),
             "product_order": right.get(key, 0.0)}
            for key in sorted(set(left) | set(right))
            if left.get(key, 0.0) != right.get(key, 0.0)]


def reactive_maps(record):
    maps = {i for edit in reaction_edits(record) for i in edit["atoms"]}
    left = {a["map_id"]: a for a in record["endpoints"]["reactant"]["atom_states"]}
    right = {a["map_id"]: a for a in record["endpoints"]["product"]["atom_states"]}
    for i in left:
        if any(left[i].get(k) != right[i].get(k)
               for k in ("formal_charge", "radical_electrons", "stereo")):
            maps.add(i)
    return maps


def proton_transfers(record):
    elements = {a["map_id"]: a["element"] for a in record["atoms"]}
    edits = reaction_edits(record)
    result = []
    for h, element in elements.items():
        if element != "H":
            continue
        broken = [next(i for i in e["atoms"] if i != h) for e in edits
                  if h in e["atoms"] and e["product_order"] == 0]
        formed = [next(i for i in e["atoms"] if i != h) for e in edits
                  if h in e["atoms"] and e["reactant_order"] == 0]
        if broken or formed:
            result.append({"hydrogen": h, "donors": broken, "acceptors": formed})
    return result


def validate(record, require_elementary=True):
    """Return the record, or fail before generating jobs (no silent repair)."""
    if record.get("schema_version") != SCHEMA_VERSION:
        raise HypothesisError("Unsupported schema_version")
    for name in ("reaction_id", "hypothesis_id", "provenance", "participants",
                 "atoms", "endpoints", "electronic_state", "interpretation"):
        if name not in record:
            raise HypothesisError("Missing {}".format(name))
    interpretation = record["interpretation"]
    if interpretation.get("kind") not in ("elementary", "net"):
        raise HypothesisError("interpretation.kind must be elementary or net")
    if require_elementary and interpretation["kind"] != "elementary":
        raise HypothesisError("Net transformations require explicit elementary-step decomposition")
    if interpretation["kind"] == "elementary" and not interpretation.get("step_id"):
        raise HypothesisError("An elementary step requires step_id")
    if record.get("hydrogen_policy") != "explicit":
        raise HypothesisError("All hydrogens, including transferring protons, must be explicit")
    participant_ids = [p["id"] for p in record["participants"]]
    if len(set(participant_ids)) != len(participant_ids):
        raise HypothesisError("Duplicate participant ID")
    atoms = record["atoms"]
    ids = [a["map_id"] for a in atoms]
    if not ids or any(type(i) is not int or i < 1 for i in ids) or len(set(ids)) != len(ids):
        raise HypothesisError("Atom maps must be unique positive integers")
    for atom in atoms:
        if atom["element"] not in ATOMIC_NUMBERS or atom["participant"] not in participant_ids:
            raise HypothesisError("Unknown element or participant")
    state = record["electronic_state"]
    if type(state.get("charge")) is not int or type(state.get("multiplicity")) is not int:
        raise HypothesisError("electronic_state_unresolved: specify integer charge and multiplicity")
    electrons = sum(ATOMIC_NUMBERS[a["element"]] for a in atoms) - state["charge"]
    if electrons < 0 or not 1 <= state["multiplicity"] <= electrons + 1 or (electrons - (state["multiplicity"] - 1)) % 2:
        raise HypothesisError("Charge/multiplicity electron parity mismatch")
    for side in ("reactant", "product"):
        endpoint = record["endpoints"][side]
        states = endpoint["atom_states"]
        if len(states) != len(ids) or {a["map_id"] for a in states} != set(ids):
            raise HypothesisError("Endpoints must contain the same unique atom maps")
        for atom in states:
            if type(atom.get("formal_charge")) is not int or type(atom.get("radical_electrons")) is not int or atom["radical_electrons"] < 0:
                raise HypothesisError("Explicit integer charge/radical state required per atom")
        if sum(a["formal_charge"] for a in states) != state["charge"]:
            raise HypothesisError("Endpoint charge does not match electronic state")
        edges = bonds(endpoint)
        if any(a not in ids or b not in ids for a, b in edges):
            raise HypothesisError("Bond refers to absent participant atom")
        # BE exporters retain nonbonding electrons for independently checkable accounting.
        if all("nonbonding_electrons" in a for a in states):
            elements = {a["map_id"]: a["element"] for a in atoms}
            for atom in states:
                nonbonding = atom["nonbonding_electrons"]
                if type(nonbonding) is not int or nonbonding < 0:
                    raise HypothesisError("Nonbonding electron counts must be nonnegative integers")
                valence = sum(v for k, v in edges.items() if atom["map_id"] in k)
                element = elements[atom["map_id"]]
                if element not in VALENCE_ELECTRONS:
                    raise HypothesisError("BE electron accounting for this element needs an explicit conversion policy")
                if abs(VALENCE_ELECTRONS[element] - valence - nonbonding - atom["formal_charge"]) > 1e-6:
                    raise HypothesisError("Atomic BE charge accounting is inconsistent")
                if atom["radical_electrons"] != nonbonding % 2:
                    raise HypothesisError("BE radical conversion disagrees with nonbonding parity")
            total = sum(a["nonbonding_electrons"] for a in states) + 2 * sum(edges.values())
            expected = sum(VALENCE_ELECTRONS[a["element"]] for a in atoms
                           if a["element"] in VALENCE_ELECTRONS) - state["charge"]
            if any(a["element"] not in VALENCE_ELECTRONS for a in atoms):
                raise HypothesisError("BE electron accounting for this element needs an explicit conversion policy")
            if abs(total - expected) > 1e-6:
                raise HypothesisError("BE valence-electron total is inconsistent")
    if "edits" in record and record["edits"] != reaction_edits(record):
        raise HypothesisError("Declared edits disagree with endpoint graphs")
    if "proton_transfers" in record and record["proton_transfers"] != proton_transfers(record):
        raise HypothesisError("Declared proton transfers disagree with endpoint graphs")
    if not reactive_maps(record):
        raise HypothesisError("No chemical change in proposal")
    return record


def load(path, require_elementary=True):
    try:
        return validate(json.loads(Path(path).read_text()), require_elementary)
    except (KeyError, TypeError) as exc:
        raise HypothesisError("Malformed hypothesis: {}".format(exc)) from exc
