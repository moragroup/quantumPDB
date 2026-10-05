"""Offline structure/benchmark audit for the progesterone CYP pilot.

Run with an interpreter containing RDKit and Biopython. This does not run QM.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem, rdBase


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def structure_audit(path, substrate):
    with gzip.open(path, "rt") as handle:
        data = MMCIF2Dict(handle)
    atom_rows = list(zip(*(data["_atom_site." + field] for field in (
        "label_comp_id", "auth_asym_id", "auth_seq_id", "label_atom_id",
        "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num"))))
    chains = sorted({r[1] for r in atom_rows if r[0] == "STR" and r[7] == "1"})
    names, elements = [], []
    for comp, name, element in zip(data["_chem_comp_atom.comp_id"],
                                   data["_chem_comp_atom.atom_id"],
                                   data["_chem_comp_atom.type_symbol"]):
        if comp == "STR" and element not in ("H", "D"):
            names.append(name)
            elements.append(element)
    mol = Chem.RWMol()
    for element in elements:
        mol.AddAtom(Chem.Atom(element))
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for comp, a, b, order in zip(*(data["_chem_comp_bond." + field] for field in
                                  ("comp_id", "atom_id_1", "atom_id_2", "value_order"))):
        if comp == "STR" and a in names and b in names:
            mol.AddBond(names.index(a), names.index(b), orders[order.upper()])
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    xyz = {r[3]: np.array(r[4:7], float) for r in atom_rows
           if r[0] == "STR" and r[1] == chains[0] and r[7] == "1"}
    conf = Chem.Conformer(len(names))
    for i, name in enumerate(names):
        conf.SetAtomPosition(i, xyz[name])
    mol.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(mol)
    matches = mol.GetSubstructMatches(substrate, useChirality=True)
    if len(matches) != 1:
        raise ValueError("STR must map uniquely to benchmark progesterone including stereo")
    ranks = list(Chem.CanonicalRankAtoms(substrate))
    mapping = [{"benchmark_atom_index_0based": i, "canonical_site_rank": int(ranks[i]),
                "ccd_atom_name": names[j]} for i, j in enumerate(matches[0])]
    distances = []
    for chain in chains:
        rows = [r for r in atom_rows if r[1] == chain and r[7] == "1"]
        fe = next(r for r in rows if r[0] == "HEM" and r[3].upper() == "FE")
        fexyz = np.array(fe[4:7], float)
        for name in ("C21", "C17", "C16", "C6", "C7"):
            carbon = next(r for r in rows if r[0] == "STR" and r[3] == name)
            distances.append({"chain": chain, "carbon": name,
                              "fe_carbon_distance_A": float(np.linalg.norm(
                                  fexyz - np.array(carbon[4:7], float)))})
    return {"pdb": path.name.split(".")[0].upper(), "atom_mapping": mapping,
            "fe_carbon_distances": distances,
            "mutation_annotation": data.get("_entity.pdbx_mutation", [])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fluxion-root", type=Path, required=True)
    parser.add_argument("--mmcif-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    benchmark = args.fluxion_root / "CYPs/CYP_benchmark_regioselectivity_train.tsv"
    split_path = args.fluxion_root / "data/publication/reproducibility/Fig3/input/data/splits/cyp_subholdout/split.json"
    archive = args.fluxion_root / "cyp_subholdout_2048_final_repeats.zip"
    with open(benchmark) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    split = json.loads(split_path.read_text())
    cases = []
    for index, accession, site in ((337, "P08686", "C21"), (333, "P05093", "C17")):
        row = rows[index]
        if row["protein"] != accession or index not in split["test"]:
            raise ValueError("Benchmark row or held-out membership changed")
        cases.append({"source_row_0based": index, "reaction_id": str(1000000 + index),
                      "source_id": row["id"], "uniprot": accession,
                      "observed_site": site, "split": "test", "substrate": row["subs"],
                      "product": row["prods"]})
    substrate = Chem.MolFromSmiles(cases[0]["substrate"])
    # The benchmark evaluator ranks the unmapped substrate. Atom-map numbers
    # otherwise influence RDKit canonical ranks and silently change site IDs.
    for atom in substrate.GetAtoms():
        atom.SetAtomMapNum(0)
    structures = []
    inputs = [benchmark, split_path, archive]
    for pdb in ("4y8w", "4nkx"):
        path = args.mmcif_root / pdb[1:3] / (pdb + ".cif.gz")
        inputs.append(path)
        structures.append(structure_audit(path, substrate))
    site_names = {str(r["canonical_site_rank"]): r["ccd_atom_name"]
                  for r in structures[0]["atom_mapping"]}
    predictions = []
    with zipfile.ZipFile(archive) as z:
        for name in sorted(z.namelist()):
            if "/srs/" not in name or not name.endswith("_sites.csv") or name.startswith("__MACOSX"):
                continue
            for row in csv.DictReader(io.StringIO(z.read(name).decode())):
                if row["reaction_id"] in {c["reaction_id"] for c in cases}:
                    predictions.append(dict(row, carbon=site_names[row["site_rank"]], archive_member=name))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {"scope": "offline preflight; no enzyme QM calculations", "rdkit_version": rdBase.rdkitVersion,
              "cases": cases, "split_counts": {k: len(v) for k, v in split.items()},
              "structures": structures, "archived_predictions": predictions,
              "input_sha256": {str(p.resolve()): digest(p) for p in inputs},
              "script_sha256": digest(__file__)}
    (args.output_dir / "preflight.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"cases": cases, "mapping": structures[0]["atom_mapping"],
                      "prediction_rows": len(predictions)}, indent=2))


if __name__ == "__main__":
    main()
