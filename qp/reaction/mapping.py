"""Persistent graph/structure/QM atom mapping and boundary checks."""

import csv
import json
import math
import shutil
from pathlib import Path

import numpy as np

from qp.reaction.schema import ATOMIC_NUMBERS, HypothesisError, bonds, reactive_maps, sha256, write_json


def read_xyz(path, last=False):
    """Strict XYZ reader; optionally select the final trajectory frame."""
    lines = Path(path).read_text().splitlines()
    frames = []
    offset = 0
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1
            continue
        n = int(lines[offset])
        block = lines[offset + 2:offset + 2 + n]
        if n < 1 or len(block) != n:
            raise HypothesisError("Incomplete XYZ frame")
        elements, xyz = [], []
        for line in block:
            fields = line.split()
            if len(fields) != 4 or fields[0] not in ATOMIC_NUMBERS:
                raise HypothesisError("Invalid XYZ atom")
            elements.append(fields[0])
            xyz.append([float(v) for v in fields[1:]])
        xyz = np.asarray(xyz, dtype=float)
        if not np.isfinite(xyz).all():
            raise HypothesisError("Non-finite XYZ coordinates")
        frames.append((elements, xyz))
        offset += n + 2
    if not frames or (len(frames) != 1 and not last):
        raise HypothesisError("Expected one XYZ frame")
    return frames[-1]


def write_xyz(path, elements, xyz, comment="persistent QM atom order"):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("{}\n{}\n".format(len(elements), comment) + "".join(
        "{} {:.10f} {:.10f} {:.10f}\n".format(e, *r) for e, r in zip(elements, xyz)))


def cluster_from_spheres(paths):
    """Match QuantumPDB's sphere concatenation order, not PDB serial numbers."""
    from Bio.PDB import PDBParser
    atoms, elements, coords = [], [], []
    for path in paths:
        model = PDBParser(QUIET=True).get_structure("cluster", str(path))[0]
        for atom in model.get_atoms():
            residue = atom.get_parent()
            key = {"sphere": Path(path).name, "chain": residue.get_parent().id,
                   "hetero": residue.id[0], "residue_number": residue.id[1],
                   "insertion_code": residue.id[2], "residue_name": residue.resname,
                   "atom_name": atom.name, "altloc": atom.altloc}
            elements.append(atom.element.title())
            coords.append(atom.coord)
            atoms.append({"qm_index": len(atoms), "structure_key": key,
                          "is_cap": residue.resname in ("ACE", "NME")})
    return elements, np.asarray(coords, dtype=float), atoms


def graph_mappings(record, structure_graph, participant=None, max_matches=128, anchors=None):
    """Enumerate topology/stereo-preserving mappings; never guess atom names.

    structure_graph uses atoms {id, element, formal_charge, radical_electrons,
    stereo?} and bonds {atoms:[id,id], order, stereo?}. Prepared ligand graphs
    can come from a CCD/SDF record; coordinates alone do not establish topology.
    IDs are caller-owned structure keys. Symmetric alternatives are returned.
    """
    source_atoms = [a for a in record["atoms"]
                    if participant is None or a["participant"] == participant]
    states = {a["map_id"]: a for a in record["endpoints"]["reactant"]["atom_states"]}
    source = {a["map_id"]: dict(a, **{k: v for k, v in states[a["map_id"]].items()
                                    if k != "map_id"}) for a in source_atoms}
    target = {a["id"]: a for a in structure_graph["atoms"]}
    if len(target) != len(structure_graph["atoms"]):
        raise HypothesisError("Duplicate prepared graph atom ID")
    if any(type(a.get("formal_charge")) is not int or
           type(a.get("radical_electrons")) is not int or a["radical_electrons"] < 0
           for a in target.values()):
        raise HypothesisError("Prepared graph requires explicit integer charge/radical states")
    left = {k: v for k, v in bonds(record["endpoints"]["reactant"]).items()
             if all(i in source for i in k)}
    right = {}
    for bond in structure_graph["bonds"]:
        pair, value = bond["atoms"], bond["order"]
        key = frozenset(pair)
        if (len(pair) != 2 or len(key) != 2 or key in right or not key <= set(target)
                or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0):
            raise HypothesisError("Invalid/duplicate prepared graph bond")
        right[key] = float(value)
    left_stereo = {tuple(sorted(b["atoms"])): b.get("stereo")
                   for b in record["endpoints"]["reactant"]["bonds"]}
    right_stereo = {frozenset(b["atoms"]): b.get("stereo") for b in structure_graph["bonds"]}
    if len(source) != len(target):
        raise HypothesisError("Participant atom inventory differs after preparation")
    candidates = {}
    anchors = anchors or {}
    if any(i not in source or j not in target for i, j in anchors.items()):
        raise HypothesisError("Mapping anchor refers to an absent chemical/prepared graph atom")
    for i, a in source.items():
        candidates[i] = [j for j, b in target.items() if all(
            a.get(k, 0 if k in ("formal_charge", "radical_electrons") else None) ==
            b.get(k, 0 if k in ("formal_charge", "radical_electrons") else None)
            for k in ("element", "formal_charge", "radical_electrons", "stereo"))]
        if i in anchors:
            candidates[i] = [j for j in candidates[i] if j == anchors[i]]
    order = sorted(source, key=lambda i: (len(candidates[i]), -sum(i in k for k in left)))
    matches = []

    def search(mapping):
        if len(matches) > max_matches:
            raise HypothesisError("Mapping symmetry exceeds limit; select explicit chemical anchors")
        if len(mapping) == len(order):
            matches.append(dict(mapping))
            return
        i = order[len(mapping)]
        for j in candidates[i]:
            if j in mapping.values():
                continue
            valid = True
            for k, other in mapping.items():
                skey, tkey = tuple(sorted((i, k))), frozenset((j, other))
                if left.get(skey, 0) != right.get(tkey, 0) or left_stereo.get(skey) != right_stereo.get(tkey):
                    valid = False
                    break
            if valid:
                mapping[i] = j
                search(mapping)
                del mapping[i]
    search({})
    if not matches:
        raise HypothesisError("No topology/stereochemistry-preserving structure mapping")
    return matches


def validate_mapping(record, model):
    """Check a selected explicit mapping against the prepared model inventory."""
    elements = model["elements"]
    if not elements or any(e not in ATOMIC_NUMBERS for e in elements):
        raise HypothesisError("Unknown prepared cluster element")
    rows = model["atom_mapping"]
    chemical = {a["map_id"]: a for a in record["atoms"]}
    if len(rows) != len(chemical) or {r["map_id"] for r in rows} != set(chemical):
        raise HypothesisError("Every chemical participant atom must map exactly once")
    indices = [r["qm_index"] for r in rows]
    if any(type(i) is not int or not 0 <= i < len(elements) for i in indices) or len(set(indices)) != len(indices):
        raise HypothesisError("QM indices must be unique, zero-based, and in range")
    for row in rows:
        if chemical[row["map_id"]]["element"] != elements[row["qm_index"]]:
            raise HypothesisError("Mapped element differs from prepared atom")
        if not row.get("structure_key"):
            raise HypothesisError("Persistent prepared structure key required")
    frozen = model.get("frozen_indices", [])
    if len(set(frozen)) != len(frozen) or any(type(i) is not int or not 0 <= i < len(elements) for i in frozen):
        raise HypothesisError("Invalid frozen boundary indices")
    mapped = {r["map_id"]: r for r in rows}
    for i in reactive_maps(record):
        if mapped[i]["qm_index"] in frozen or mapped[i].get("is_cap"):
            raise HypothesisError("Reactive atoms cannot be frozen or replaced by caps")
    if type(model.get("charge")) is not int or type(model.get("environment_charge", 0)) is not int:
        raise HypothesisError("electronic_state_unresolved: integer cluster/environment charges required")
    if model["charge"] != record["electronic_state"]["charge"] + model.get("environment_charge", 0):
        raise HypothesisError("Cluster charge does not equal chemical plus environment charge")
    if type(model.get("multiplicity")) is not int or model["multiplicity"] < 1:
        raise HypothesisError("electronic_state_unresolved: cluster multiplicity required")
    ne = sum(ATOMIC_NUMBERS[e] for e in elements) - model["charge"]
    if ne < 0 or model["multiplicity"] > ne + 1 or (ne - model["multiplicity"] + 1) % 2:
        raise HypothesisError("Cluster charge/spin parity inconsistent")
    return model


def build_cluster_model(record, sphere_paths, specification, output, source_dir=None):
    """Map audited participant graphs onto QuantumPDB's actual capped atom order.

    Each graph atom carries its actual ``qm_index``. Topology, atom/bond stereo,
    charge and radicals are matched before assigning chemical map IDs. Symmetric
    matches require explicit chemical-map-to-graph-ID anchors in the recipe.
    """
    output = Path(output).resolve()
    if (output / "model.json").exists():
        raise HypothesisError("Mapped model already exists; choose a new output directory")
    output.mkdir(parents=True, exist_ok=True)
    source_dir = Path(source_dir or ".").resolve()
    elements, xyz, structure_atoms = cluster_from_spheres(sphere_paths)
    rows, alternatives = [], {}
    graph_participants = [g["participant"] for g in specification["graphs"]]
    if len(set(graph_participants)) != len(graph_participants):
        raise HypothesisError("One prepared graph per participant required")
    for graph in specification["graphs"]:
        for atom in graph["atoms"]:
            i = atom["qm_index"]
            if type(i) is not int or not 0 <= i < len(elements) or elements[i] != atom["element"]:
                raise HypothesisError("Prepared graph does not match actual sphere atom order/elements")
        anchors = {int(i): key for i, key in graph.get("anchors", {}).items()}
        matches = graph_mappings(record, graph, graph["participant"], anchors=anchors)
        alternatives[graph["participant"]] = matches
        write_json(output / "mapping_alternatives.json", alternatives)
        if len(matches) != 1:
            raise HypothesisError("Ambiguous/absent mapping for {}; add chemical anchors (see mapping_alternatives.json)".format(graph["participant"]))
        atoms = {a["id"]: a for a in graph["atoms"]}
        for chemical_map, key in matches[0].items():
            i = atoms[key]["qm_index"]
            rows.append(dict(structure_atoms[i], map_id=chemical_map, graph_atom_id=key))
    prepared = {key: specification[key] for key in
                ("charge", "multiplicity", "environment_charge", "frozen_indices", "covalent_environment_links")
                if key in specification}
    prepared.update(elements=elements, atom_mapping=rows, reactant_xyz="reactant.xyz",
                    mapping_provenance={"method": "exact participant topology/stereo with explicit symmetry anchors",
                                        "sphere_sha256": {str(Path(p).resolve()): sha256(p) for p in sphere_paths},
                                        "prepared_graphs": specification["graphs"]})
    validate_mapping(record, prepared)
    write_xyz(output / "reactant.xyz", elements, xyz)
    for key in ("product_xyz", "ts_xyz"):
        if specification.get(key):
            target = key.replace("_xyz", ".xyz")
            shutil.copyfile(source_dir / specification[key], output / target)
            prepared[key] = target
    write_json(output / "model.json", prepared)
    write_mapping(output / "atom_mapping.tsv", record, prepared)
    return prepared


def write_mapping(path, record, model):
    chemical = {a["map_id"]: a for a in record["atoms"]}
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["map_id", "element", "participant",
                                "qm_index_0based", "qm_index_1based", "movable", "structure_key"], delimiter="\t")
        writer.writeheader()
        for row in model["atom_mapping"]:
            i = row["qm_index"]
            atom = chemical[row["map_id"]]
            writer.writerow({"map_id": row["map_id"], "element": atom["element"],
                             "participant": atom["participant"], "qm_index_0based": i,
                             "qm_index_1based": i + 1, "movable": i not in model.get("frozen_indices", []),
                             "structure_key": json.dumps(row["structure_key"], sort_keys=True)})
