"""Prepare and evaluate the six cofactor-retained QuantumPDB examples.

python -m qp.analyze.dynamics_benchmark --help
"""

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from scipy.spatial import cKDTree

from qp.analyze.catalytic_movement import (
    SIDECHAIN, atom_rmsd, fit_coordinates, read_table, write_table,
)

EXAMPLES = ("P00183", "P09788", "Q46AN5", "Q40577", "P0A2K1", "P95480")


class FirstModel(Select):
    def accept_model(self, model):
        return model.id == 0

    def accept_atom(self, atom):
        if not atom.is_disordered():
            return True
        alternatives = atom.parent[atom.name].disordered_get_list()
        chosen = min(alternatives, key=lambda a: (
            -(a.occupancy or 0), 0 if a.altloc == " " else 1 if a.altloc == "A" else 2, a.altloc))
        return atom is chosen


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atoms_by_residue(path):
    model = PDBParser(QUIET=True).get_structure("input", str(path))[0]
    return {(r.parent.id, str(r.id[1]) + r.id[2].strip()): r for r in model.get_residues()}


def prepare(args):
    root = Path(args.output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    pairs = read_table(Path(args.results_dir) / "pairs.tsv", "\t")
    residues = read_table(Path(args.results_dir) / "residues.tsv", "\t")
    cases, targets, hashes = [], [], {}
    for accession in EXAMPLES:
        eligible = [p for p in pairs if p["uniprot_id"] == accession and p["status"] == "complete"]
        if not eligible:
            raise ValueError("No complete pair for " + accession)
        pair = max(eligible, key=lambda p: (float(p["max_sidechain_rmsd_A"]), p["seq_group"]))
        case = dict(pair, case_id=accession)
        selected = [r for r in residues if r["seq_group"] == pair["seq_group"]]
        for state in ("substrate_free", "substrate_bound"):
            pdb = pair[state + "_pdb"].lower()
            work = root / accession / state
            out = work / "clusters"
            source = out / pdb / (pdb + ".pdb")
            source.parent.mkdir(parents=True, exist_ok=True)
            cif = Path(pair[state + "_file"])
            with gzip.open(cif, "rt") as handle:
                structure = MMCIFParser(QUIET=True).get_structure(pdb, handle)
            writer = PDBIO()
            writer.set_structure(structure)
            writer.save(str(source), FirstModel())
            hashes[str(cif)] = digest(cif)
            hashes[str(source)] = digest(source)
            chain = pair[state + "_auth_chain"]
            by_residue = atoms_by_residue(source)
            centers = []
            for row in selected:
                pos = row[state + "_author_position"]
                residue = by_residue[(chain, pos)]
                if residue.id[2].strip():
                    raise ValueError("Center syntax cannot encode insertion code: " + pos)
                centers.append("{}_{}{}".format(residue.resname, chain, residue.id[1]))
                targets.append(dict(row, case_id=accession, state=state,
                                    auth_chain=chain, author_pos=pos))
            protein = [a.coord for r in structure[0][chain] if r.id[0] == " " for a in r]
            tree = cKDTree(protein)
            for residue in structure[0].get_residues():
                if residue.id[0] == " " or residue.resname not in pair["cofactors"].split(";"):
                    continue
                if min(tree.query([a.coord for a in residue])[0]) <= 4.0:
                    if residue.id[2].strip():
                        raise ValueError("Cofactor center has insertion code")
                    centers.append("{}_{}{}".format(residue.resname, residue.parent.id, residue.id[1]))
            # A dash forces exact-match mode in CenterResidue (even for one token).
            center = "-".join(centers) if len(centers) > 1 else centers[0] + "-" + centers[0]
            electronic_review = accession in {"P00183", "P09788"}
            input_csv = work / "input.csv"
            with input_csv.open("w", newline="") as handle:
                writer_csv = csv.DictWriter(handle, fieldnames=["pdb_id", "center", "oxidation", "multiplicity"])
                writer_csv.writeheader()
                writer_csv.writerow(dict(pdb_id=pdb, center=center, oxidation=0,
                                         multiplicity="" if electronic_review else 1))
            config = dict(input=str(input_csv), output_dir=str(out), modeller=True,
                          optimize_select_residues=1, protoss=True, coordination=True,
                          skip="all", number_of_spheres=2, radius_of_first_sphere=4.0,
                          additional_ligands=sorted(set((pair["cofactors"] + ";" + pair["substrates"]).split(";"))),
                          include_ligands=2, capping_method=2, compute_charges=True,
                          count_residues=True, write_xyz=True, smoothing_method=2,
                          qm_program="orca", method="b3lyp", basis="def2-SVP",
                          optimization=True, dispersion="D3BJ", aux_basis="def2/J",
                          dielectric=10, use_implicit_solvent=True, charge_embedding=False,
                          scheduler="slurm", partition="cpu", nprocs=16, memory="64G",
                          time_limit="2-00:00:00", create_jobs=True, submit_jobs=False)
            (work / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
            case[state + "_config"] = str(work / "config.yaml")
            case[state + "_input_pdb"] = str(source)
        case["electronic_review"] = "required" if accession in {"P00183", "P09788"} else "closed_shell_starting_assumption"
        cases.append(case)
    write_table(root / "cases.tsv", cases)
    write_table(root / "targets.tsv", targets)
    manifest = dict(script_sha256=digest(__file__), inputs=hashes,
                    selection="maximum complete-pair catalytic sidechain RMSD per example",
                    interpretation="retrospective endpoint benchmark, not a dynamics simulation",
                    python=sys.version)
    for name in ("pairs.tsv", "residues.tsv"):
        manifest["inputs"][str(Path(args.results_dir).resolve() / name)] = digest(Path(args.results_dir) / name)
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    evaluate(argparse.Namespace(output_dir=str(root), outputs=None))


def run(args):
    root = Path(args.output_dir).resolve()
    if args.stage == "run" and importlib.util.find_spec("modeller") is None:
        raise RuntimeError("Modeller is missing. Install it in this Python environment before preparation.")
    for case in read_table(root / "cases.tsv", "\t"):
        for state in ("substrate_free", "substrate_bound"):
            config = Path(case[state + "_config"])
            if args.stage == "submit":
                data = yaml.safe_load(config.read_text())
                electronic = read_table(data["input"])
                if any(not r["multiplicity"] or not r["oxidation"] for r in electronic):
                    raise ValueError("Specify reviewed oxidation/multiplicity in " + data["input"])
            # CLI has no python -m entrypoint; invoke Click explicitly.
            command = [sys.executable, "-c", "from qp.cli import cli; cli()", args.stage, "-c", str(config)]
            log = config.parent / (args.stage + ".log")
            with log.open("w") as handle:
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
            text = log.read_text()
            if result.returncode or "CRITICAL FAILURE" in text or " errors:" in text:
                raise RuntimeError("Pipeline failed; inspect " + str(log))
            if args.stage == "run" and not list(Path(yaml.safe_load(config.read_text())["output_dir"]).glob("*/*/*.xyz")):
                raise RuntimeError("No clusters produced; inspect " + str(log))


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    count = int(lines[0])
    rows = [line.split() for line in lines[2:] if line.strip()]
    if len(rows) != count or any(len(row) != 4 for row in rows):
        raise ValueError("Expected one final XYZ frame in " + str(path))
    coordinates = np.array([[float(v) for v in row[1:]] for row in rows])
    if not np.isfinite(coordinates).all():
        raise ValueError("Non-finite XYZ coordinates")
    return [r[0].upper() for r in rows], coordinates


def cluster_atoms(cluster):
    spheres = sorted(Path(cluster).glob("[0-9]*.pdb"), key=lambda p: int(p.stem))
    if not spheres:
        raise ValueError("No numbered sphere PDB files in " + str(cluster))
    parser = PDBParser(QUIET=True)
    return [a for sphere in spheres for a in parser.get_structure(sphere.stem, str(sphere))[0].get_atoms()]


def measure_output(cluster, final_xyz, chain, position, code):
    """Measure relaxation, preserving the pipeline's sphere/XYZ atom order."""
    atoms = cluster_atoms(cluster)
    elements, final = read_xyz(final_xyz)
    if elements != [a.element.upper() for a in atoms]:
        raise ValueError("Final XYZ atom count/order does not match numbered spheres")
    initial = np.array([a.coord for a in atoms], dtype=float)
    target_indices = [i for i, a in enumerate(atoms) if a.parent.parent.id == chain
                      and str(a.parent.id[1]) + a.parent.id[2].strip() == position
                      and a.parent.resname == code]
    if not target_indices:
        return dict(status="residue_absent")
    anchors = [i for i, a in enumerate(atoms) if a.name == "CA" and i not in target_indices]
    if len(anchors) < 3 or np.linalg.matrix_rank(initial[anchors] - initial[anchors].mean(axis=0)) < 2:
        return dict(status="insufficient_fit_anchors")
    rotation, translation, fit = fit_coordinates(final[anchors], initial[anchors])
    before = {atoms[i].name: initial[i] for i in target_indices}
    after = {atoms[i].name: final[i] for i in target_indices}
    names = set(SIDECHAIN.get(code, "").split()) or {"CA"}
    if not names <= before.keys():
        return dict(status="incomplete_sidechain")
    return dict(status="measured", simulated_sidechain_relaxation_A=atom_rmsd(
        before, after, names, rotation, translation, code), cluster_fit_rmsd_A=fit)


def evaluate(args):
    root = Path(args.output_dir).resolve()
    cases = {c["case_id"]: c for c in read_table(root / "cases.tsv", "\t")}
    targets = read_table(root / "targets.tsv", "\t")
    outputs = read_table(args.outputs) if args.outputs else []
    rows = []
    for target in targets:
        case = cases[target["case_id"]]
        config = yaml.safe_load(Path(case[target["state"] + "_config"]).read_text())
        row = {k: target[k] for k in ("case_id", "state", "uniprot_position", "residue_code", "sidechain_rmsd_A")}
        row.update(protocol="hydrogen_only_optimization" if config["optimization"] else "single_point",
                   capability="cannot_sample_heavy_atom_motion",
                   status="not_run", predicted_heavy_atom_motion_A=0.0,
                   interpretation="Heavy atoms fixed; separate endpoint calculations do not demonstrate recovery")
        pdb = case[target["state"] + "_pdb"].lower()
        clusters = list((Path(config["output_dir"]) / pdb).glob("*/*.xyz"))
        covered = []
        for xyz in clusters:
            atoms = cluster_atoms(xyz.parent)
            present = {a.name for a in atoms if a.parent.parent.id == target["auth_chain"]
                       and str(a.parent.id[1]) + a.parent.id[2].strip() == target["author_pos"]
                       and a.parent.resname == target["residue_code"]}
            expected = set(SIDECHAIN.get(target["residue_code"], "").split()) or {"CA"}
            if expected <= present:
                covered.append(str(xyz.parent))
        row["clusters_with_complete_sidechain"] = ";".join(covered)
        if clusters:
            row["status"] = "cluster_prepared" if covered else "sidechain_not_in_clusters"
        matches = [o for o in outputs if o["case_id"] == target["case_id"] and o["state"] == target["state"]]
        if len(matches) > 1:
            raise ValueError("Provide only one selected cluster output per case/state")
        if matches:
            output = matches[0]
            if output["converged"].lower() != "true":
                row["status"] = "not_converged"
            else:
                row.update(measure_output(output["cluster_dir"], output["final_xyz"],
                                          target["auth_chain"], target["author_pos"], target["residue_code"]))
                if row.get("simulated_sidechain_relaxation_A", 0) > 0.01:
                    row["interpretation"] = "Motion conflicts with fixed-heavy-atom protocol; check actual input and atom mapping"
        rows.append(row)
    write_table(root / "evaluation.tsv", rows)
    summary = dict(examples=len(cases), endpoint_inputs=len(cases) * 2,
                   outputs_supplied=len(outputs), heavy_atom_dynamics_supported=False,
                   conclusion="Current QuantumPDB optimizations freeze catalytic heavy atoms and cannot recover their rearrangements.",
                   modeller_available=importlib.util.find_spec("modeller") is not None,
                   executable_paths={p: shutil.which(p) for p in ("sbatch", "terachem", "orca")})
    (root / "evaluation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    prep.add_argument("--results-dir", default="substrate_movement_results/all_mcsa")
    prep.add_argument("--output-dir", default="dynamics_benchmark_results")
    prep.set_defaults(function=prepare)
    execution = commands.add_parser("run")
    execution.add_argument("--output-dir", default="dynamics_benchmark_results")
    execution.add_argument("--stage", choices=("run", "submit"), default="run")
    execution.set_defaults(function=run)
    evaluation = commands.add_parser("evaluate")
    evaluation.add_argument("--output-dir", default="dynamics_benchmark_results")
    evaluation.add_argument("--outputs", help="CSV: case_id,state,cluster_dir,final_xyz,converged")
    evaluation.set_defaults(function=evaluate)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
