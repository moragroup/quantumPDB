"""Stage a boundary-capped CYP21A2 Compound I model and spin diagnostics.

The charge is a declared fragment-state starting hypothesis, not an inferred
SCF state. Single points must be reviewed before optimization or HAT searches.
"""
import argparse
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO

from qp.cluster.spheres import cap_chains, write_pdbs
from qp.reaction.jobs import write_job
from qp.reaction.mapping import cluster_from_spheres, write_xyz
from qp.reaction.schema import sha256, write_json


def oxo_position(fe, sulfur, distance=1.65):
    direction = np.asarray(fe, float) - np.asarray(sulfur, float)
    if np.linalg.norm(direction) < 1:
        raise ValueError("Invalid axial Fe-S geometry")
    return np.asarray(fe, float) + distance * direction / np.linalg.norm(direction)


def stage(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError("Compound I model already exists; do not overwrite")
    structure = PDBParser(QUIET=True).get_structure("4y8w", str(source))
    chain = structure[0]["A"]
    selected = {r for r in chain if r.resname in ("HEM", "STR", "FE")
                or (r.id[0] == " " and r.id[1] in (429, 92, 427, 288))}
    expected = {("CYS", 429), ("ARG", 92), ("ARG", 427), ("ASP", 288)}
    if not expected <= {(r.resname, r.id[1]) for r in selected}:
        raise ValueError("Missing prescribed catalytic-environment residues")
    cysteine = next(r for r in selected if r.resname == "CYS")
    if "HG" in cysteine:
        raise ValueError("Axial cysteine must be a thiolate")
    protonation_changes = []
    for res in selected:
        if res.resname == "ASP":
            for name in ("HD2", "HOD1", "HOD2"):
                if name in res:
                    protonation_changes.append({"residue": "ASP_A288", "removed_atom": name,
                                                "reason": "Declared deprotonated Asp starting hypothesis; sensitivity required"})
                    res.detach_child(name)
    for res in selected:
        if res.resname == "ARG" and not {"HE", "HH11", "HH12", "HH21", "HH22"} <= {a.name for a in res}:
            raise ValueError("Starting charge requires protonated ARG")
    caps = cap_chains(structure[0], selected, 2)
    output.mkdir(parents=True)
    writer = PDBIO()
    writer.set_structure(structure)
    write_pdbs(writer, selected | caps, str(output / "cluster.pdb"))
    elements, xyz, atoms = cluster_from_spheres([output / "cluster.pdb"])
    fe_index = elements.index("Fe")
    sulfur = next(i for i, a in enumerate(atoms) if a["structure_key"]["residue_name"] == "CYS"
                  and a["structure_key"]["atom_name"] == "SG")
    oxo = oxo_position(xyz[fe_index], xyz[sulfur])
    distances = np.linalg.norm(xyz - oxo, axis=1)
    collisions = [i for i, d in enumerate(distances) if i != fe_index and d < 1.25]
    if collisions:
        raise ValueError("Inserted oxo clashes with existing atoms: " + str(collisions))
    elements = elements + ["O"]
    coords = np.vstack((xyz, oxo))
    frozen = [i for i, a in enumerate(atoms) if a["is_cap"] or (
        a["structure_key"]["residue_name"] in ("CYS", "ARG", "ASP") and
        a["structure_key"]["atom_name"] in ("N", "CA", "C", "O"))]
    write_xyz(output / "compound_I.xyz", elements, coords)
    # HEM propionates -2, porphyrin dianion -2 + radical-cation +1;
    # Fe(IV) +4, oxo -2, thiolate -1, two Arg +2, Asp -1 => -1.
    fragments = {"heme_propionates": -2, "porphyrin_dianion": -2,
                 "porphyrin_radical_cation": 1, "Fe_IV": 4, "oxo": -2,
                 "axial_thiolate": -1, "two_arginines": 2, "aspartate": -1}
    from qp.reaction.schema import ATOMIC_NUMBERS
    electrons = sum(ATOMIC_NUMBERS[e] for e in elements) - sum(fragments.values())
    if electrons % 2 != 1:
        raise ValueError("Compound I candidate must have odd electron count")
    settings = {"nprocs": 16, "maxcore_mb": 3900, "walltime": "24:00:00",
                "job_timeout_seconds": 86400, "method": "UKS B3LYP"}
    jobs = []
    for multiplicity in (2, 4):
        directory = output / ("spin_{}".format(multiplicity))
        model = {"elements": elements, "charge": -1, "multiplicity": multiplicity,
                 "frozen_indices": frozen}
        # Generate launch scripts with the shared writer, then use an explicit
        # single point: optimization waits for SCF/spin/coordination diagnostics.
        engine = dict(settings, method="B3LYP", basis="def2-SVP", keywords=["UKS", "D3BJ", "TightSCF"], dielectric=10)
        write_job(directory, "endpoint_opt", model, engine)
        batch = (directory / "submit.sh").read_text().replace(
            "#SBATCH --cpus-per-task=16", "#SBATCH --ntasks=16\n#SBATCH --cpus-per-task=1")
        (directory / "submit.sh").write_text(batch)
        inp = (directory / "job.inp").read_text().split("%geom")[0]
        inp = inp.replace("TightSCF Opt", "TightSCF SP")
        inp += "%scf MaxIter 300 end\n* xyzfile -1 {} input.xyz\n".format(multiplicity)
        (directory / "job.inp").write_text(inp)
        write_xyz(directory / "input.xyz", elements, coords)
        jobs.append({"multiplicity": multiplicity, "directory": str(directory),
                     "input_sha256": sha256(directory / "job.inp")})
    write_json(output / "model_audit.json", {
        "status": "provisional_Compound_I_spin_screen", "atom_count": len(elements),
        "electron_count": electrons, "charge": -1, "fragment_charge_ledger": fragments,
        "source": str(source), "source_sha256": sha256(source),
        "protonation_changes": protonation_changes,
        "selected_residues": sorted([r.resname + "_A" + str(r.id[1]) for r in selected]),
        "atoms": atoms, "oxo_qm_index": len(atoms), "oxo_position_A": oxo.tolist(),
        "frozen_indices": frozen, "jobs": jobs,
        "required_review": ["SCF convergence", "Fe/oxo/porphyrin spin populations", "S_squared",
                            "Alternative broken-symmetry initial guesses", "Protonation and cluster sensitivity"],
        "not_yet_validated_for_HAT": True})
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    print(stage(args.source, args.output))


if __name__ == "__main__":
    main()
