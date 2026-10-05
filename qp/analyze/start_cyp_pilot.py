"""Freeze the CYP21A2 pilot and stage its QuantumPDB preparation inputs."""
import argparse
import gzip
from pathlib import Path

import yaml
from Bio.PDB import MMCIFParser, PDBIO

from qp.analyze.dynamics_benchmark import FirstModel
from qp.reaction.schema import sha256, write_json


class ChainA(FirstModel):
    def accept_chain(self, chain):
        return chain.id == "A"


def prepare(root, fluxion, mmcif):
    root = Path(root).resolve()
    fluxion = Path(fluxion).resolve()
    source = Path(mmcif) / "y8/4y8w.cif.gz"
    work = root / "models/P08686/chain_A"
    pdb = work / "clusters/4y8w/4y8w.pdb"
    if pdb.exists():
        raise ValueError("Preparation already staged; resume the existing config")
    pdb.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rt") as handle:
        structure = MMCIFParser(QUIET=True).get_structure("4y8w", handle)
    chain = structure[0]["A"]
    assert chain[("H_HEM", 603, " ")].resname == "HEM"
    assert chain[("H_STR", 604, " ")].resname == "STR"
    assert chain[(" ", 429, " ")].resname == "CYS"
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(pdb), ChainA())
    csv = work / "input.csv"
    csv.write_text("pdb_id,center,oxidation,multiplicity\n4y8w,HEM_A603-STR_A604-CYS_A429,,\n")
    config = dict(input=str(csv), output_dir=str(work / "clusters"), modeller=True,
                  optimize_select_residues=1, protoss=True, coordination=True,
                  skip="all", number_of_spheres=1, radius_of_first_sphere=4.0,
                  additional_ligands=["HEM", "STR"], include_ligands=2,
                  capping_method=2, compute_charges=True, count_residues=True,
                  write_xyz=True, smoothing_method=2, create_jobs=False,
                  submit_jobs=False)
    (work / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    checkpoint = fluxion / "experiments/flow/checkpoints/cyp_reps/srs/s00/best_acc.pt"
    split = fluxion / "data/publication/reproducibility/Fig3/input/data/splits/cyp_subholdout/split.json"
    write_json(root / "protocol.json", {
        "status": "preparation_started", "enzyme": "P08686", "pdb": "4Y8W", "chain": "A",
        "reaction_id": "1000337", "sites": ["C21", "C17"], "spin_multiplicities": [2, 4],
        "checkpoint_selection": "archived SRS s00 best_acc; no substitution",
        "checkpoint": str(checkpoint), "checkpoint_available": checkpoint.is_file(),
        "split": str(split), "split_sha256": sha256(split),
        "source_structure": str(source.resolve()), "source_sha256": sha256(source),
        "preparation_config": str(work / "config.yaml"),
        "method": "B3LYP", "basis": "def2-SVP", "dispersion": "D3BJ", "dielectric": 10,
        "max_initial_concurrent_jobs": 2, "nprocs_per_job": 16,
        "memory_per_job_GB": 64, "walltime_hours": 24,
        "gates": ["Complete protonated cluster", "Audited Compound I charge and spin",
                  "Mapped HAT hypotheses", "Exact checkpoint provenance"],
        "preparation_is_not_reaction_validation": True})
    return work


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="regioselectivity_evaluation_results")
    parser.add_argument("--fluxion-root", default="../Fluxion")
    parser.add_argument("--mmcif-root", default="/mnt/labs/shared/enzymes/mmCIF")
    args = parser.parse_args()
    print(prepare(args.output, args.fluxion_root, args.mmcif_root))


if __name__ == "__main__":
    main()
