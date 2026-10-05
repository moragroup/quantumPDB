"""Compare wild-type substrate-free/bound chains with the same cofactor retained.

Ligand identities must be specified explicitly; absence means absence of a
modeled contacting ligand, rather than proof of absence in the experiment.
"""

import argparse
import csv
import json
import platform
from collections import defaultdict
from itertools import product
from pathlib import Path

import Bio
import numpy as np
import scipy
from Bio import SeqIO

from qp.analyze.catalytic_movement import (
    compare_group, contacts_ligand, identifiers, load_structure, read_table,
    sequence_mapping, sha256, write_table,
)


def classify_chain(chain, ligands, reference, substrate, cofactor, excluded, cutoff):
    """Require full reference coverage, wild-type sequence and cofactor contact."""
    _, identity, coverage = sequence_mapping(reference, chain["sequence"])
    if identity < 1 or coverage < 1:
        return "excluded", "mutation or incomplete reference sequence", identity, coverage
    if not contacts_ligand(chain, ligands, {cofactor}, cutoff):
        return "excluded", "required cofactor not contacting chain", identity, coverage
    if contacts_ligand(chain, ligands, excluded, cutoff):
        return "excluded", "excluded product/analogue contacting chain", identity, coverage
    bound = contacts_ligand(chain, ligands, {substrate}, cutoff)
    return ("substrate_bound" if bound else "substrate_free"), "", identity, coverage


def descriptive_fields(row):
    """Use substrate terminology in outputs, including chain/position fields."""
    return {k.replace("apo_", "substrate_free_").replace("holo_", "substrate_bound_"): v
            for k, v in row.items()}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniprot", required=True)
    parser.add_argument("--substrate", required=True, help="mmCIF component code, e.g. TRP")
    parser.add_argument("--cofactor", required=True, help="exact mmCIF component code, e.g. FAD")
    parser.add_argument("--exclude-ligand", action="append", default=[], help="product/analogue code")
    parser.add_argument("--entity-states", type=Path, required=True)
    parser.add_argument("--mmcif-root", type=Path, required=True)
    parser.add_argument("--catalytic-residues", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    args = parser.parse_args(argv)
    if args.contact_cutoff <= 0:
        parser.error("contact cutoff must be positive")
    references = {}
    for record in SeqIO.parse(args.reference_fasta, "fasta"):
        parts = record.id.split("|")
        references[parts[2] if parts[0] == "mcsa" and len(parts) >= 3
                   else parts[1] if len(parts) >= 3 else parts[0]] = str(record.seq)
    reference = references.get(args.uniprot)
    catalytic = [r for r in read_table(args.catalytic_residues)
                 if r["uniprot_id"] == args.uniprot and r["is_reference"].lower() == "true"]
    if not reference or not catalytic:
        parser.error("accession needs a reference sequence and reference catalytic annotations")
    rows = [r for r in read_table(args.entity_states, "\t")
            if args.uniprot in identifiers(r["uniprot_ids"])]
    by_pdb = defaultdict(list)
    for row in rows:
        by_pdb[row["pdb_id"]].append(row)
    inventory, selected, structures = [], defaultdict(dict), {}
    for pdb, entities in sorted(by_pdb.items()):
        chains, ligands, path = load_structure(args.mmcif_root, pdb, {r["entity_id"] for r in entities})
        structures[str(path)] = sha256(path)
        for label, chain in sorted(chains.items()):
            state, reason, identity, coverage = classify_chain(
                chain, ligands, reference, args.substrate, args.cofactor,
                set(args.exclude_ligand), args.contact_cutoff)
            entity = next(r for r in entities if r["entity_id"] == chain["entity"])
            inventory.append({"pdb_id": pdb, "chain": label, "auth_chain": chain["auth_chain"],
                              "state": state, "reason": reason, "identity": identity,
                              "coverage": coverage, "resolution_A": entity["resolution_A"],
                              "entity_ligands": entity["ligands"]})
            if state != "excluded":
                selected[state].setdefault(pdb, []).append((label, entity))
    # Reuse the same alignment/geometry implementation, with explicit chain
    # restrictions and substrate-specific states in place of global apo/holo.
    args.state_column = "state"
    args.min_identity, args.min_coverage, args.min_fit_atoms = 1.0, 1.0, 20
    pairs, residues = [], []
    for free, bound in product(sorted(selected["substrate_free"]), sorted(selected["substrate_bound"])):
        if free == bound:
            continue
        f, b = selected["substrate_free"][free], selected["substrate_bound"][bound]
        group = {"uniprot_ids": args.uniprot, "seq_group": free + "_" + bound,
                 "enzyme_name": f[0][1]["enzyme_name"], "best_apo": free, "best_holo": bound,
                 "best_apo_res_A": f[0][1]["resolution_A"], "best_holo_res_A": b[0][1]["resolution_A"],
                 "apo_chains": ";".join(x[0] for x in f), "holo_chains": ";".join(x[0] for x in b)}
        states = defaultdict(list)
        for pdb, entries, state in ((free, f, "apo"), (bound, b, "holo")):
            for entity_id in sorted({r["entity_id"] for _, r in entries}):
                row = dict(next(r for _, r in entries if r["entity_id"] == entity_id))
                row.update(state=state, seq_group=group["seq_group"], ligands=args.substrate if state == "holo" else "")
                states[pdb].append(row)
        pair, measured = compare_group(group, catalytic, reference, states, args)
        pairs.append(descriptive_fields(pair))
        residues.extend(descriptive_fields(r) for r in measured)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_table(args.output_dir / "structure_inventory.tsv", inventory)
    write_table(args.output_dir / "pairs.tsv", pairs)
    write_table(args.output_dir / "residues.tsv", residues)
    pair_index = {r["seq_group"]: r for r in pairs}
    export = []
    for residue in residues:
        pair = pair_index[residue["seq_group"]]
        annotations = [r for r in catalytic if int(r["residue_pos"]) == residue["uniprot_position"]]
        row = dict(pair, **residue)
        row.update(pair_status=pair["status"], residue_status=residue["status"],
                   catalytic_roles=";".join(sorted({r["roles"] for r in annotations})),
                   substrate=args.substrate, cofactor=args.cofactor)
        export.append(row)
    fields = list(dict.fromkeys(k for r in export for k in r)) or ["status"]
    with open(args.output_dir / "all_catalytic_residues.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        writer.writerows(export)
    summary = {"uniprot_id": args.uniprot, "substrate": args.substrate, "cofactor": args.cofactor,
               "pairs": len(pairs), "complete_pairs": sum(p["status"] == "complete" for p in pairs),
               "substrate_free_structures": sorted(selected["substrate_free"]),
               "substrate_bound_structures": sorted(selected["substrate_bound"]),
               "interpretation": "Static structural differences; not a measurement of time-dependent dynamics."}
    manifest = {"parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                "versions": {"python": platform.python_version(), "numpy": np.__version__,
                             "biopython": Bio.__version__, "scipy": scipy.__version__},
                "script_sha256": sha256(__file__),
                "geometry_script_sha256": sha256(Path(__file__).with_name("catalytic_movement.py")),
                "inputs": {k: {"path": str(getattr(args, k).resolve()), "sha256": sha256(getattr(args, k))}
                           for k in ("entity_states", "catalytic_residues", "reference_fasta")},
                "structures": structures}
    for name, value in (("summary.json", summary), ("manifest.json", manifest)):
        with open(args.output_dir / name, "w") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
