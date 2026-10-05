"""Engine evidence, reaction connectivity, and conservative cost accounting."""

import csv
import re
from pathlib import Path

import numpy as np

from qp.reaction.mapping import read_xyz
from qp.reaction.schema import HypothesisError, bonds, reaction_edits


HARTREE_TO_KCAL_MOL = 627.509474
# Covalent radii, Angstrom. Metal bond perception is intentionally not automatic.
RADII = {"H": .31, "B": .84, "C": .76, "N": .71, "O": .66, "F": .57,
         "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39}


def parse_orca(path, mode):
    text = Path(path).read_text(errors="replace")
    energies = re.findall(r"FINAL SINGLE POINT ENERGY\s+([-+\d.Ee]+)", text)
    frequencies = [(int(i), float(f)) for i, f in re.findall(
        r"^\s*(\d+):\s*([-+\d.]+)\s*cm(?:\*\*-1|\*\*1|\^-1|-1)", text, re.M)]
    normal = "ORCA TERMINATED NORMALLY" in text
    optimizer = "THE OPTIMIZATION HAS CONVERGED" in text
    scf_failed = any(s in text.upper() for s in ("SCF NOT CONVERGED", "SCF CONVERGENCE FAILURE", "SCF DID NOT CONVERGE"))
    runtime = re.search(r"TOTAL RUN TIME:\s*(\d+) days\s*(\d+) hours\s*(\d+) minutes\s*(\d+) seconds\s*(\d+) msec", text)
    seconds = None if runtime is None else sum(int(v) * w for v, w in zip(runtime.groups(), (86400, 3600, 60, 1, .001)))
    converged = normal and not scf_failed
    if mode in ("endpoint_opt", "ts_opt", "path_search", "constrained_scan"):
        converged = converged and optimizer
    if mode == "path_search":
        converged = converged and "THE NEB OPTIMIZATION HAS CONVERGED" in text
    if mode == "constrained_scan":
        # Every scan point needs a converged constrained minimum.
        points = len(re.findall(r"(?im)^\s*(?:RELAXED SURFACE SCAN STEP|RELAXED SURFACE SCAN:)", text))
        converged = converged and (not points or text.count("THE OPTIMIZATION HAS CONVERGED") >= points)
    if mode == "frequency":
        converged = converged and bool(frequencies)
    if mode == "irc":
        # A summary is also printed after exhausting the iteration limit.
        converged = converged and text.upper().count("THE IRC HAS CONVERGED") >= 2
    displaced = re.search(r"List of displaced atoms\s*\.\.\.\s*([\d\s]+?)(?=\n\s*\n|\n\s*[A-Za-z])", text)
    return {"normal_termination": normal, "converged": converged,
            "optimizer_converged": optimizer, "scf_failure": scf_failed,
            "energy_hartree": float(energies[-1]) if energies else None,
            "frequencies_cm_1": [{"index": i, "frequency": f} for i, f in frequencies],
            "engine_wall_seconds": seconds,
            "displaced_atoms_0based": [int(i) - 1 for i in displaced.group(1).split()] if displaced else None,
            "cost": {"energy_evaluations_observed": len(energies),
                     "gradient_evaluations_observed": len(re.findall(r"CARTESIAN GRADIENT", text)),
                     "hessian_evaluations_observed": len(re.findall(r"(?m)^\s*(?:CARTESIAN HESSIAN|HESSIAN MATRIX)\s*$", text)),
                     "counter_semantics": "output markers, lower bounds; engine-internal calls may not be printed"},
            "engine_version": (re.search(r"Program Version\s+(\S+)", text).group(1)
                               if re.search(r"Program Version\s+(\S+)", text) else None)}


def normal_modes(path):
    """Read ORCA's blocked $normal_modes matrix (columns are mode indices)."""
    lines = Path(path).read_text().splitlines()
    start = lines.index("$normal_modes") + 1
    while not lines[start].strip():
        start += 1
    nr, nc = map(int, lines[start].split())
    result = np.zeros((nr, nc))
    row = start + 1
    while row < len(lines) and not lines[row].startswith("$"):
        fields = lines[row].split()
        row += 1
        if not fields or fields[0].startswith("#"):
            continue
        columns = list(map(int, fields))
        for _ in range(nr):
            values = lines[row].split()
            row += 1
            i = int(values[0])
            for j, value in zip(columns, values[1:]):
                result[i, j] = float(value.replace("D", "E"))
    return result


def mode_overlap(record, model, xyz, displacement):
    """Projection onto independent changed-bond stretch directions."""
    index = {r["map_id"]: r["qm_index"] for r in model["atom_mapping"]}
    vectors = []
    for edit in reaction_edits(record):
        a, b = (index[i] for i in edit["atoms"])
        delta = xyz[b] - xyz[a]
        length = np.linalg.norm(delta)
        if length < 1e-8:
            raise HypothesisError("Overlapping reaction atoms")
        vector = np.zeros_like(xyz)
        vector[a], vector[b] = -delta / length, delta / length
        vectors.append(vector.ravel())
    if not vectors:
        return 0.0
    u, s, _ = np.linalg.svd(np.asarray(vectors).T, full_matrices=False)
    basis = u[:, s > 1e-8]
    v = np.asarray(displacement).ravel()
    return float(np.linalg.norm(basis.T @ v) / max(np.linalg.norm(v), 1e-12))


def perceived_edges(record, model, xyz, scale=1.25):
    """Distance-based connectivity only: does not establish bond order/stereo."""
    rows = {r["map_id"]: r["qm_index"] for r in model["atom_mapping"]}
    atoms = record["atoms"]
    if any(a["element"] not in RADII for a in atoms):
        raise HypothesisError("Metal connectivity requires audited endpoint graph evidence")
    edges = set()
    for offset, a in enumerate(atoms):
        for b in atoms[offset + 1:]:
            distance = np.linalg.norm(xyz[rows[a["map_id"]]] - xyz[rows[b["map_id"]]])
            if distance < .35:
                raise HypothesisError("Unphysical overlapping endpoint atoms")
            if distance <= scale * (RADII[a["element"]] + RADII[b["element"]]):
                edges.add(tuple(sorted((a["map_id"], b["map_id"]))))
    return edges


def endpoint_matches(record, model, xyz, side, audited_graph=None):
    expected = bonds(record["endpoints"][side])
    if audited_graph is not None:
        states = lambda graph: {a["map_id"]: a for a in graph["atom_states"]}
        stereo = lambda graph: {tuple(sorted(b["atoms"])): b.get("stereo") for b in graph["bonds"]}
        return (bonds(audited_graph) == expected and states(audited_graph) == states(record["endpoints"][side])
                and stereo(audited_graph) == stereo(record["endpoints"][side]))
    # For simple single-bond non-stereochemical reactions this is sufficient
    # connectivity evidence. Other chemistry needs explicit endpoint graph audit.
    states = record["endpoints"][side]["atom_states"]
    if (record["provenance"].get("stereochemistry_audit_required") or
            any(v != 1 for v in expected.values()) or any(a.get("stereo") for a in states)
            or any(b.get("stereo") for b in record["endpoints"][side]["bonds"])):
        raise HypothesisError("Bond-order/stereo validation requires audited connected endpoint graphs")
    return perceived_edges(record, model, xyz) == set(expected)


def validate_step(record, model, evidence, imaginary_cutoff=-50., minimum_overlap=.2):
    result = {"status": "not_converged", "reasons": [], "evidence": evidence,
              "validation_scope": "movable_subspace" if model.get("frozen_indices") else "full_cluster",
              "barrier_type": "electronic_cluster_energy", "barrier_units": "kcal/mol"}
    required = ("reactant", "product", "ts", "reactant_frequency", "product_frequency", "ts_frequency", "irc", "connected_a", "connected_b")
    if any(name not in evidence or not evidence[name].get("converged") for name in required):
        result["reasons"].append("Missing or unconverged stationary-point/IRC evidence")
        return result
    if any(side not in evidence[side].get("matches", []) for side in ("reactant", "product")):
        result["status"] = "wrong_connectivity"
        result["reasons"].append("Optimized initial endpoints differ from proposed graphs or need a graph audit")
        return result
    for name in ("reactant_frequency", "product_frequency"):
        if any(f["frequency"] < imaginary_cutoff for f in evidence[name]["frequencies_cm_1"]):
            result["reasons"].append("Endpoint is not a minimum in evaluated subspace")
            return result
    imaginary = [f for f in evidence["ts_frequency"]["frequencies_cm_1"] if f["frequency"] < imaginary_cutoff]
    if len(imaginary) != 1 or evidence["ts_frequency"].get("reaction_mode_overlap", 0) < minimum_overlap:
        result["reasons"].append("TS needs exactly one significant reaction-aligned imaginary mode")
        return result
    matches = [evidence[k].get("matches", []) for k in ("connected_a", "connected_b")]
    if not (("reactant" in matches[0] and "product" in matches[1]) or
            ("product" in matches[0] and "reactant" in matches[1])):
        result["status"] = "wrong_connectivity"
        result["reasons"].append("Optimized IRC endpoints do not match both proposed graphs")
        return result
    energies = [evidence[s].get("energy_hartree") for s in ("reactant", "product", "ts")]
    if any(e is None for e in energies):
        result["reasons"].append("Missing comparable energies")
        return result
    er, ep, et = energies
    result.update(status="validated_step", forward_barrier=(et - er) * HARTREE_TO_KCAL_MOL,
                  reverse_barrier=(et - ep) * HARTREE_TO_KCAL_MOL,
                  reaction_energy=(ep - er) * HARTREE_TO_KCAL_MOL,
                  imaginary_frequencies_cm_1=imaginary)
    return result


def write_cost_report(path, records):
    fields = ["job", "mode", "status", "cached", "wall_seconds", "queue_seconds",
              "allocated_cpu_hours", "allocated_gpu_hours", "energy_evaluations_observed",
              "gradient_evaluations_observed", "hessian_evaluations_observed"]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in records:
            writer.writerow(dict(row, **row.get("cost", {})))
