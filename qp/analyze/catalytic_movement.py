"""Reproducible, offline M-CSA catalytic-residue apo/holo comparison.

Run ``python -m qp.analyze.catalytic_movement --help`` for inputs/options.
Coordinates are in Angstroms. Results describe structural differences, not
proof that ligand binding caused them. One best-resolution pair per input
group is used; this is not an exhaustive search across all deposited pairs.
"""

import argparse
import csv
import gzip
import hashlib
import json
import platform
import re
from collections import defaultdict
from pathlib import Path

import Bio
import numpy as np
import scipy
from Bio import Align, SeqIO
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


SIDECHAIN = {
    "ALA": "CB", "ARG": "CB CG CD NE CZ NH1 NH2", "ASN": "CB CG OD1 ND2",
    "ASP": "CB CG OD1 OD2", "CYS": "CB SG", "GLN": "CB CG CD OE1 NE2",
    "GLU": "CB CG CD OE1 OE2", "GLY": "", "HIS": "CB CG ND1 CD2 CE1 NE2",
    "ILE": "CB CG1 CG2 CD1", "LEU": "CB CG CD1 CD2", "LYS": "CB CG CD CE NZ",
    "MET": "CB CG SD CE", "PHE": "CB CG CD1 CD2 CE1 CE2 CZ", "PRO": "CB CG CD",
    "SER": "CB OG", "THR": "CB OG1 CG2", "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "CB CG1 CG2",
}
SYMMETRY = {
    "ASP": [("OD1", "OD2")], "GLU": [("OE1", "OE2")], "ARG": [("NH1", "NH2")],
    "PHE": [("CD1", "CD2"), ("CE1", "CE2")],
    "TYR": [("CD1", "CD2"), ("CE1", "CE2")],
    "VAL": [("CG1", "CG2")], "LEU": [("CD1", "CD2")],
}


def identifiers(text):
    return {x.strip() for x in re.split(r"[,;]", text) if x.strip()}


def read_table(path, delimiter=","):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sequence_mapping(reference, sequence):
    """Map 1-based UniProt positions to 1-based polymer positions.

    Align complete polymer sequences, including unresolved residues, rather
    than relying on author numbering. Terminal truncations/tags are free.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, sequence)[0]
    mapping = {}
    for (a, b), (c, d) in zip(*alignment.aligned):
        mapping.update((int(i + 1), int(j + 1)) for i, j in zip(range(a, b), range(c, d)))
    identical = sum(reference[i - 1] == sequence[j - 1] for i, j in mapping.items())
    return mapping, identical / max(1, len(mapping)), len(mapping) / max(1, len(reference))


def fit_coordinates(mobile, fixed):
    """Least-squares proper rotation with row-vector coordinates."""
    mobile, fixed = np.asarray(mobile), np.asarray(fixed)
    mc, fc = mobile.mean(axis=0), fixed.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (fixed - fc))
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    translation = fc - mc @ rotation
    rmsd = float(np.sqrt(np.mean(np.sum((mobile @ rotation + translation - fixed) ** 2, axis=1))))
    return rotation, translation, rmsd


def atom_rmsd(apo, holo, names, rotation, translation, code):
    """Minimize over chemically equivalent atom labels (coupled aromatic swaps)."""
    names = sorted(names)
    if not names:
        return None
    target = np.array([apo[n] for n in names])
    permutations = [dict(zip(names, names))]
    swaps = SYMMETRY.get(code, [])
    if swaps and all(a in names and b in names for a, b in swaps):
        swapped = dict(zip(names, names))
        for a, b in swaps:
            swapped[a], swapped[b] = b, a
        permutations.append(swapped)
    return min(float(np.sqrt(np.mean(np.sum(
        (np.array([holo[p[n]] for n in names]) @ rotation + translation - target) ** 2, axis=1
    )))) for p in permutations)


def load_structure(root, pdb, entity_ids):
    path = Path(root) / pdb[1:3].lower() / (pdb.lower() + ".cif.gz")
    with gzip.open(path, "rt") as handle:
        data = MMCIF2Dict(handle)
    sequences = dict(zip(data["_entity_poly.entity_id"],
                         ["".join(s.split()) for s in data["_entity_poly.pdbx_seq_one_letter_code_can"]]))
    chain_entities = dict(zip(data["_struct_asym.id"], data["_struct_asym.entity_id"]))
    chains = {c: {"entity": e, "sequence": sequences[e], "residues": {}, "auth_chain": ""}
              for c, e in chain_entities.items() if e in entity_ids and e in sequences}
    fields = ["label_asym_id", "label_seq_id", "label_comp_id", "label_atom_id", "type_symbol",
              "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "label_alt_id", "pdbx_PDB_model_num",
              "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code"]
    n = len(data["_atom_site.label_asym_id"])
    columns = [data.get("_atom_site." + f, ["?"] * n) for f in fields]
    first_model = columns[10][0]
    ligands = []
    for c, pos, code, atom, element, x, y, z, occupancy, alt, model, auth, author_pos, ins in zip(*columns):
        if model != first_model or element.upper() in {"H", "D"}:
            continue
        xyz = np.array([float(x), float(y), float(z)])
        if pos in {".", "?"} or chain_entities.get(c) not in sequences:
            ligands.append((code, xyz))
            continue
        if c not in chains:
            continue
        chains[c]["auth_chain"] = auth
        residue = chains[c]["residues"].setdefault(int(pos), {
            "code": code, "author_pos": author_pos + (ins if ins not in {".", "?"} else ""),
            "atoms": {}, "choices": {},
        })
        # Highest occupancy per atom; deterministic blank/A/lexical altloc tie break.
        rank = (-float(occupancy) if occupancy not in {".", "?"} else 0,
                0 if alt in {".", "?"} else 1 if alt == "A" else 2, alt)
        if atom not in residue["choices"] or rank < residue["choices"][atom]:
            residue["choices"][atom] = rank
            residue["atoms"][atom] = xyz
    return chains, ligands, path


def ligand_codes(text):
    codes = set()
    for item in text.split(";"):
        code = item.split(":")[0]
        if code.startswith("oligo["):
            codes.update(code[6:-1].split("+"))
        elif code:
            codes.add(code)
    return codes


def contacts_ligand(chain, ligands, codes, cutoff):
    from scipy.spatial import cKDTree
    protein = [xyz for r in chain["residues"].values() for xyz in r["atoms"].values()]
    ligand = [xyz for code, xyz in ligands if code in codes]
    return bool(protein and ligand and np.min(cKDTree(protein).query(ligand)[0]) <= cutoff)


def compare_group(group, catalytic, reference, states, args):
    uni = group["uniprot_ids"]
    base = {"uniprot_id": uni, "seq_group": group.get("seq_group", ""),
            "enzyme_name": group["enzyme_name"], "apo_pdb": group["best_apo"],
            "holo_pdb": group["best_holo"], "apo_resolution_A": group["best_apo_res_A"],
            "holo_resolution_A": group["best_holo_res_A"],
            "mcsa_ids": ";".join(sorted({r["mcsa_id"] for r in catalytic})),
            "status": "", "reason": ""}
    if not catalytic or not reference:
        base.update(status="unavailable", reason="missing catalytic annotation or reference sequence")
        return base, []
    wanted = {int(r["residue_pos"]): r["residue_code"].upper() for r in catalytic}
    if any(wanted[int(r["residue_pos"])] != r["residue_code"].upper() for r in catalytic):
        base.update(status="unavailable", reason="conflicting M-CSA annotations at one position")
        return base, []
    if any(p < 1 or p > len(reference) or protein_letters_3to1_extended.get(code) != reference[p - 1]
           for p, code in wanted.items()):
        base.update(status="unavailable", reason="M-CSA residue does not match reference sequence")
        return base, []
    structures, candidates = [], []
    for state in ("apo", "holo"):
        pdb = base[state + "_pdb"]
        rows = [r for r in states.get(pdb, []) if uni in identifiers(r["uniprot_ids"])
                and r[args.state_column] == state
                and (not base["seq_group"] or r["seq_group"] == base["seq_group"])]
        if not rows:
            base.update(status="unavailable", reason="no matching " + state + " entity")
            return base, []
        chains, ligands, path = load_structure(args.mmcif_root, pdb, {r["entity_id"] for r in rows})
        structures.append(str(path))
        eligible = []
        for label, chain in sorted(chains.items()):
            allowed_chains = group.get(state + "_chains")
            if allowed_chains and label not in identifiers(allowed_chains):
                continue
            mapping, identity, coverage = sequence_mapping(reference, chain["sequence"])
            if identity < args.min_identity or coverage < args.min_coverage:
                continue
            codes = set().union(*(ligand_codes(r["ligands"] if args.state_column == "state"
                                               else r["organic_ligands"])
                                  for r in rows if r["entity_id"] == chain["entity"]))
            if state == "holo" and not contacts_ligand(chain, ligands, codes, args.contact_cutoff):
                continue
            eligible.append((label, chain, mapping, identity, coverage))
        candidates.append(eligible)
    best = None
    for a in candidates[0]:
        for h in candidates[1]:
            common = sorted((set(a[2]) & set(h[2])) - set(wanted))
            fit_pos = [p for p in common if reference[p - 1] == a[1]["sequence"][a[2][p] - 1]
                       == h[1]["sequence"][h[2][p] - 1]
                       and "CA" in a[1]["residues"].get(a[2][p], {}).get("atoms", {})
                       and "CA" in h[1]["residues"].get(h[2][p], {}).get("atoms", {})]
            if len(fit_pos) < args.min_fit_atoms:
                continue
            rotation, translation, rmsd = fit_coordinates(
                [h[1]["residues"][h[2][p]]["atoms"]["CA"] for p in fit_pos],
                [a[1]["residues"][a[2][p]]["atoms"]["CA"] for p in fit_pos])
            residues = []
            for p, code in sorted(wanted.items()):
                row = {"uniprot_id": uni, "seq_group": base["seq_group"], "uniprot_position": p,
                       "residue_code": code, "apo_pdb": base["apo_pdb"], "holo_pdb": base["holo_pdb"],
                       "apo_chain": a[0], "holo_chain": h[0], "status": "missing", "reason": "unresolved/unmapped"}
                ar = a[1]["residues"].get(a[2].get(p))
                hr = h[1]["residues"].get(h[2].get(p))
                if ar and hr:
                    row.update(apo_author_position=ar["author_pos"], holo_author_position=hr["author_pos"])
                    if ar["code"] != code or hr["code"] != code:
                        row.update(reason="catalytic mutation/modified residue")
                    else:
                        aa, ha = ar["atoms"], hr["atoms"]
                        side = set(SIDECHAIN.get(code, "").split()) or {"CA"}
                        required = side | {"N", "CA", "C", "O"}
                        row.update(status="ok" if required <= aa.keys() & ha.keys() else "partial",
                                   reason="" if required <= aa.keys() & ha.keys() else "missing catalytic heavy atoms",
                                   ca_displacement_A=atom_rmsd(aa, ha, {"CA"} & aa.keys() & ha.keys(),
                                                               rotation, translation, code),
                                   sidechain_rmsd_A=atom_rmsd(aa, ha, side, rotation, translation, code)
                                   if side <= aa.keys() & ha.keys() else None,
                                   heavy_atom_rmsd_A=atom_rmsd(aa, ha, required, rotation, translation, code)
                                   if required <= aa.keys() & ha.keys() else None)
                residues.append(row)
            rank = (-sum(r["status"] == "ok" for r in residues), -len(fit_pos), rmsd, a[0], h[0])
            if best is None or rank < best[0]:
                best = (rank, a, h, rmsd, len(fit_pos), residues)
    if best is None:
        base.update(status="unavailable", reason="no eligible chain pair with enough aligned C-alpha atoms")
        return base, []
    _, a, h, rmsd, nfit, residues = best
    complete = all(r["status"] == "ok" for r in residues)
    base.update(status="complete" if complete else "partial", apo_chain=a[0], holo_chain=h[0],
                apo_auth_chain=a[1]["auth_chain"], holo_auth_chain=h[1]["auth_chain"],
                apo_identity=a[3], holo_identity=h[3], apo_coverage=a[4], holo_coverage=h[4],
                fit_ca_atoms=nfit, fit_ca_rmsd_A=rmsd, n_catalytic_residues=len(wanted),
                n_complete_residues=sum(r["status"] == "ok" for r in residues),
                apo_file=structures[0], holo_file=structures[1])
    for metric in ("ca_displacement_A", "sidechain_rmsd_A", "heavy_atom_rmsd_A"):
        values = [r[metric] for r in residues if r.get(metric) is not None]
        base["max_" + metric] = max(values) if values else None
    return base, residues


def write_table(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fields or ["status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_movement_csv(path, pairs, residues, groups, catalytic, threshold=None):
    """Export all annotated residues, or thresholded residues from complete pairs.

    With no threshold, retain missing measurements as blanks with explicit
    statuses; unavailable pairs still get one row per annotated position.
    """
    key = lambda r: (r.get("uniprot_id", r.get("uniprot_ids")), r.get("seq_group", ""))
    pair_index = {key(r): r for r in pairs}
    group_index = {key(r): r for r in groups}
    fields = ["uniprot_id", "enzyme_name", "ec_numbers", "mcsa_ids", "mcsa_ec_numbers",
              "seq_group", "residue_code", "uniprot_position", "catalytic_roles",
              "apo_pdb", "apo_chain", "apo_auth_chain", "apo_author_position",
              "holo_pdb", "holo_chain", "holo_auth_chain", "holo_author_position",
              "sidechain_rmsd_A", "ca_displacement_A", "heavy_atom_rmsd_A",
              "fit_ca_rmsd_A", "apo_resolution_A", "holo_resolution_A",
              "movement_threshold_A", "pair_status", "residue_status", "reason"]
    export_residues = list(residues)
    if threshold is None:
        represented = {key(r) for r in residues}
        for pair in pairs:
            if key(pair) in represented:
                continue
            annotations = catalytic[pair["uniprot_id"]]
            positions = sorted({(int(r["residue_pos"]), r["residue_code"].upper()) for r in annotations})
            for position, code in positions:
                export_residues.append({"uniprot_id": pair["uniprot_id"],
                                        "seq_group": pair.get("seq_group", ""),
                                        "uniprot_position": position, "residue_code": code,
                                        "status": "unavailable", "reason": pair.get("reason", "")})
    count = 0
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        for residue in export_residues:
            pair = pair_index[key(residue)]
            value = residue.get("sidechain_rmsd_A")
            if threshold is not None and (pair["status"] != "complete" or value in (None, "")
                                          or float(value) < threshold):
                continue
            annotations = [r for r in catalytic[residue["uniprot_id"]]
                           if int(r["residue_pos"]) == int(residue["uniprot_position"])]
            row = {f: residue.get(f, pair.get(f, "")) for f in fields}
            group = group_index[key(residue)]
            row.update(enzyme_name=pair.get("enzyme_name") or group["enzyme_name"],
                       ec_numbers=group["ec_numbers"],
                       mcsa_ids=";".join(sorted({r["mcsa_id"] for r in annotations})),
                       mcsa_ec_numbers=";".join(sorted(set().union(*(
                           identifiers(r["ec_numbers"]) for r in annotations)))),
                       catalytic_roles=";".join(sorted({r["roles"] for r in annotations if r["roles"]})),
                       movement_threshold_A=threshold if threshold is not None else "",
                       pair_status=pair["status"], residue_status=residue["status"])
            writer.writerow(row)
            count += 1
    return count


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=Path, required=True, help="apo_holo_groups[_uniprot/_strict].tsv")
    parser.add_argument("--entity-states", type=Path, required=True)
    parser.add_argument("--mmcif-root", type=Path, required=True)
    parser.add_argument("--mcsa-annotations", type=Path, required=True)
    parser.add_argument("--catalytic-residues", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--thresholds", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    parser.add_argument("--csv-threshold", type=float, default=1.0,
                        help="minimum side-chain RMSD for moving_residues.csv (default: 1 A)")
    parser.add_argument("--state-column", choices=["state", "state_strict"], default="state")
    parser.add_argument("--min-identity", type=float, default=0.95)
    parser.add_argument("--min-coverage", type=float, default=0.5)
    parser.add_argument("--min-fit-atoms", type=int, default=20)
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--limit", type=int, help="optional small validation run")
    parser.add_argument("--uniprot", action="append", help="restrict to accession (repeatable)")
    args = parser.parse_args(argv)
    if (any(t <= 0 for t in args.thresholds) or args.csv_threshold <= 0 or not 0 <= args.min_identity <= 1
            or not 0 <= args.min_coverage <= 1 or args.min_fit_atoms < 3
            or args.contact_cutoff <= 0 or (args.limit is not None and args.limit < 1)):
        parser.error("thresholds/cutoff/limit must be positive, fractions in [0,1], fit atoms >=3")
    annotations = read_table(args.mcsa_annotations)
    mcsa_unis = set().union(*(identifiers(r["uniprot_ids"]) for r in annotations))
    groups = [r for r in read_table(args.groups, "\t") if r["uniprot_ids"] in mcsa_unis]
    if args.uniprot:
        groups = [r for r in groups if r["uniprot_ids"] in args.uniprot]
    groups.sort(key=lambda r: (r["uniprot_ids"], r.get("seq_group", "")))
    if args.limit:
        groups = groups[:args.limit]
    catalytic = defaultdict(list)
    for row in read_table(args.catalytic_residues):
        if row["is_reference"].lower() == "true":
            catalytic[row["uniprot_id"]].append(row)
    references = {}
    for record in SeqIO.parse(args.reference_fasta, "fasta"):
        parts = record.id.split("|")
        references[parts[2] if len(parts) >= 3 and parts[0] == "mcsa"
                   else parts[1] if len(parts) >= 3 else parts[0]] = str(record.seq)
    states = defaultdict(list)
    needed = {r[k] for r in groups for k in ("best_apo", "best_holo")}
    for row in read_table(args.entity_states, "\t"):
        if row["pdb_id"] in needed:
            states[row["pdb_id"]].append(row)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs, residues = [], []
    for i, group in enumerate(groups, 1):
        uni = group["uniprot_ids"]
        try:
            pair, rows = compare_group(group, catalytic[uni], references.get(uni), states, args)
        except (OSError, ValueError, KeyError) as error:
            pair, rows = {"uniprot_id": uni, "seq_group": group.get("seq_group", ""),
                          "apo_pdb": group["best_apo"], "holo_pdb": group["best_holo"],
                          "status": "error", "reason": str(error)}, []
        pairs.append(pair)
        residues.extend(rows)
        if i % 25 == 0 or i == len(groups):
            print("Compared {}/{} groups".format(i, len(groups)), flush=True)
    write_table(args.output_dir / "pairs.tsv", pairs)
    write_table(args.output_dir / "residues.tsv", residues)
    all_rows = write_movement_csv(args.output_dir / "all_catalytic_residues.csv", pairs, residues,
                                 groups, catalytic)
    print("Wrote {} catalytic residues to {}".format(all_rows, args.output_dir / "all_catalytic_residues.csv"))
    csv_rows = write_movement_csv(args.output_dir / "moving_residues.csv", pairs, residues,
                                 groups, catalytic, args.csv_threshold)
    print("Wrote {} moving residues to {}".format(csv_rows, args.output_dir / "moving_residues.csv"))
    summary = {"groups": len(groups), "unique_uniprot": len({r["uniprot_ids"] for r in groups}),
               "status_counts": {s: sum(r["status"] == s for r in pairs)
                                 for s in ("complete", "partial", "unavailable", "error")},
               "movement": {}}
    for metric in ("ca_displacement_A", "sidechain_rmsd_A", "heavy_atom_rmsd_A"):
        summary["movement"][metric] = {}
        for threshold in sorted(set(args.thresholds)):
            complete = [r for r in pairs if r["status"] == "complete"]
            moved = [r for r in complete if r.get("max_" + metric) is not None
                     and r["max_" + metric] >= threshold]
            summary["movement"][metric][str(threshold)] = {
                "complete_groups": len(complete), "moved_groups": len(moved),
                "complete_unique_uniprot": len({r["uniprot_id"] for r in complete}),
                "moved_unique_uniprot": len({r["uniprot_id"] for r in moved}),
            }
    inputs = {name: {"path": str(getattr(args, name).resolve()), "sha256": sha256(getattr(args, name))}
              for name in ("groups", "entity_states", "mcsa_annotations", "catalytic_residues", "reference_fasta")}
    files = sorted({str(args.mmcif_root / pdb[1:3].lower() / (pdb.lower() + ".cif.gz"))
                    for pdb in needed
                    if (args.mmcif_root / pdb[1:3].lower() / (pdb.lower() + ".cif.gz")).exists()})
    manifest = {"parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                "versions": {"python": platform.python_version(), "numpy": np.__version__,
                             "biopython": Bio.__version__, "scipy": scipy.__version__},
                "script_sha256": sha256(__file__), "inputs": inputs,
                "structures": {p: sha256(p) for p in files}}
    for name, value in (("summary.json", summary), ("manifest.json", manifest)):
        with open(args.output_dir / name, "w") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
