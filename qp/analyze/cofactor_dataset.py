"""Dataset-wide M-CSA substrate comparisons with retained cofactors.

Discover substrates by exact ChEBI cross-references to M-CSA reactants.
Cache PDBe chemical-component annotations for offline reproducibility.
"""

import argparse
import csv
import json
import platform
import time
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from itertools import product
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError

from Bio import SeqIO
import Bio
import numpy as np
import scipy

from qp.analyze import catalytic_movement as geometry
from qp.analyze.substrate_movement import descriptive_fields


# Exact component identities, not interchangeable oxidation states/analogues.
ORGANIC_COFACTORS = set("FAD FDA FMN NAD NAP NAI NDP HEM HEA HEC PLP PMP TPP B12 COB SF4 FES F3S".split())
METALS = {"FE": "iron", "ZN": "zinc", "MG": "magnesium", "MN": "manganese",
          "CU": "copper", "CO": "cobalt", "NI": "nickel", "CA": "calcium",
          "MO": "molybdenum", "W": "tungsten"}


def fetch_json(url):
    for attempt in range(4):
        try:
            with urlopen(url, timeout=90) as handle:
                return json.load(handle)
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def component_chebi(component):
    return {r["resource_id"].split(":")[-1] for r in (component.get("cross_links") or [])
            if r["resource"] == "ChEBI"}


def ligand_roles(compounds, components, cofactors, chebi=None):
    """Exact reaction mappings; a product-only component is never a substrate."""
    reactants = {str(r["chebi_id"]) for r in compounds if r["type"] == "reactant"}
    products = {str(r["chebi_id"]) for r in compounds if r["type"] == "product"}
    chebi = chebi or {}
    def chemical_key(data):
        key = data.get("default_structure", {}).get("standard_inchi_key") if data.get("default_structure") else None
        return "-".join(key.split("-")[:2]) if key else None
    reactant_keys = {chemical_key(chebi[c]) for c in reactants if c in chebi} - {None}
    product_keys = {chemical_key(chebi[c]) for c in products if c in chebi} - {None}
    substrates, product_codes = set(), set()
    for code, component in components.items():
        ids = component_chebi(component)
        key = component.get("inchi_key")
        key = "-".join(key.split("-")[:2]) if key else None
        if code in cofactors or code in ORGANIC_COFACTORS or code in METALS:
            continue
        if (ids & reactants or key in reactant_keys) and (component.get("weight") or 0) > 50:
            substrates.add(code)
        elif ids & products or key in product_keys:
            product_codes.add(code)
    return substrates, product_codes


def classify_contacts(codes, substrates, cofactors, artifacts):
    """Only cofactor(s), or cofactor(s)+known substrate(s), after artifacts."""
    relevant = codes - artifacts
    retained = relevant & cofactors
    bound = relevant & substrates
    unknown = relevant - cofactors - substrates
    if not retained:
        return "excluded", retained, bound, "no supported cofactor contact"
    if unknown:
        return "excluded", retained, bound, "other ligand: " + ";".join(sorted(unknown))
    return ("substrate_bound" if bound else "substrate_free"), retained, bound, ""


def protein_summary(pairs):
    grouped = defaultdict(list)
    for pair in pairs:
        grouped[pair["uniprot_id"]].append(pair)
    output = []
    for uni, rows in sorted(grouped.items()):
        complete = [r for r in rows if r["status"] == "complete"]
        row = {"uniprot_id": uni, "enzyme_name": rows[0]["enzyme_name"],
               "ec_numbers": rows[0]["ec_numbers"], "pairs": len(rows),
               "complete_pairs": len(complete),
               "substrates": ";".join(sorted({r["substrates"] for r in rows})),
               "cofactors": ";".join(sorted({r["cofactors"] for r in rows}))}
        for metric in ("ca_displacement_A", "sidechain_rmsd_A", "heavy_atom_rmsd_A"):
            values = [float(r["max_" + metric]) for r in complete if r.get("max_" + metric) not in (None, "")]
            row["max_" + metric] = max(values) if values else ""
        output.append(row)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("entity-states", "mmcif-root", "catalytic-residues", "reference-fasta",
                 "mcsa-root", "artifacts", "output-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8, help="annotation downloads only")
    parser.add_argument("--uniprot", action="append")
    parser.add_argument("--discover-only", action="store_true")
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache = args.output_dir / "annotation_cache"
    cache.mkdir(exist_ok=True)
    catalytic = defaultdict(list)
    for row in geometry.read_table(args.catalytic_residues):
        if row["is_reference"].lower() == "true":
            catalytic[row["uniprot_id"]].append(row)
    references = {}
    for record in SeqIO.parse(args.reference_fasta, "fasta"):
        parts = record.id.split("|")
        references[parts[2] if parts[0] == "mcsa" and len(parts) >= 3
                   else parts[1] if len(parts) >= 3 else parts[0]] = str(record.seq)
    rows_by_uni = defaultdict(list)
    for row in geometry.read_table(args.entity_states, "\t"):
        uni = row["uniprot_ids"]
        if uni in catalytic and uni in references and (not args.uniprot or uni in args.uniprot):
            rows_by_uni[uni].append(row)
    all_cofactors = ORGANIC_COFACTORS | set(METALS)
    candidate_unis = {uni for uni, rows in rows_by_uni.items()
                      if any(geometry.ligand_codes(r["ligands"]) & all_cofactors for r in rows)}
    observed = {uni: set().union(*(geometry.ligand_codes(r["ligands"]) | geometry.identifiers(r["additives"])
                                  for r in rows_by_uni[uni]))
                for uni in candidate_unis}
    codes = sorted(set().union(*observed.values()))
    batches = [codes[i:i + 40] for i in range(0, len(codes), 40)]

    def download_batch(batch):
        import hashlib
        path = cache / ("components_" + hashlib.sha256(",".join(batch).encode()).hexdigest()[:16] + ".json")
        if not path.exists():
            try:
                data = fetch_json("https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/" + ",".join(batch))
            except HTTPError as error:
                if error.code != 404:
                    raise
                data = {}
                for code in batch:
                    try:
                        data.update(fetch_json("https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/" + code))
                    except HTTPError as error:
                        if error.code != 404:
                            raise
            path.write_text(json.dumps(data, indent=2))
        return json.loads(path.read_text())

    components = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for data in pool.map(download_batch, batches):
            components.update({k: v[0] for k, v in data.items() if v})
    print("Cached {} component annotations for {} candidate proteins".format(len(components), len(candidate_unis)), flush=True)

    def download_uniprot(uni):
        path = cache / ("uniprot_" + uni + ".json")
        if not path.exists():
            data = fetch_json("https://rest.uniprot.org/uniprotkb/" + uni + ".json")
            path.write_text(json.dumps(data, indent=2))
        return uni, json.loads(path.read_text())

    # Metals require enzyme-specific UniProt cofactor support, rather than
    # treating every crystallization calcium/magnesium ion as a cofactor.
    supported = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for uni, data in pool.map(download_uniprot, sorted(candidate_unis)):
            text = " ".join(c["name"] for r in data.get("comments", [])
                            if r["commentType"] == "COFACTOR" for c in r.get("cofactors", [])).lower()
            supported[uni] = ORGANIC_COFACTORS | {k for k, v in METALS.items()
                                                if v in text or k.lower() + "(" in text}
    reactions = defaultdict(list)
    reaction_hashes = {}
    for path in sorted(args.mcsa_root.glob("*/entry_*_api.json")):
        entry = json.loads(path.read_text())
        for uni in candidate_unis:
            if str(entry["mcsa_id"]) in {r["mcsa_id"] for r in catalytic[uni]}:
                reactions[uni].extend(entry.get("reaction", {}).get("compounds", []))
                reaction_hashes[str(path)] = geometry.sha256(path)
    chebi_ids = sorted({str(r["chebi_id"]) for compounds in reactions.values() for r in compounds})
    def download_chebi(identifier):
        path = cache / ("chebi_" + identifier + ".json")
        if not path.exists():
            try:
                data = fetch_json("https://www.ebi.ac.uk/chebi/backend/api/public/compound/" + identifier + "/?format=json")
            except HTTPError as error:
                if error.code != 404:
                    raise
                data = {"unavailable": True}
            path.write_text(json.dumps(data, indent=2))
        return identifier, json.loads(path.read_text())
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        chebi = dict(pool.map(download_chebi, chebi_ids))
    print("Cached {} reaction compound structures".format(len(chebi)), flush=True)
    with open(args.artifacts) as handle:
        artifacts = {line.split("#")[0].strip().split()[0] for line in handle
                     if line.split("#")[0].strip()}
    specifications, discovery = {}, []
    for uni in sorted(candidate_unis):
        cofactors = observed[uni] & supported[uni]
        substrates, products = ligand_roles(reactions[uni], {c: components[c] for c in observed[uni]
                                                            if c in components}, cofactors, chebi)
        # An explicitly mapped substrate cannot be ignored as an artifact.
        ignored = artifacts - substrates - cofactors - products
        # Entity summaries merge contacts across chains. Do not reject a
        # candidate here because another chain has substrate/other ligands.
        cofactor_rows = [r for r in rows_by_uni[uni] if geometry.ligand_codes(r["ligands"]) & cofactors]
        possible_bound = any((geometry.ligand_codes(r["ligands"]) | geometry.identifiers(r["additives"])) & substrates
                             for r in cofactor_rows)
        possible_free = len({r["pdb_id"] for r in cofactor_rows}) >= 2
        discovery.append({"uniprot_id": uni, "enzyme_name": rows_by_uni[uni][0]["enzyme_name"],
                          "cofactors": ";".join(sorted(cofactors)), "substrates": ";".join(sorted(substrates)),
                          "products": ";".join(sorted(products)), "candidate_pair": possible_free and possible_bound})
        if possible_free and possible_bound:
            specifications[uni] = (substrates, cofactors, ignored)
    geometry.write_table(args.output_dir / "discovery.tsv", discovery)
    print("{} proteins have candidate cofactor-only and substrate+cofactor entities".format(len(specifications)), flush=True)
    if args.discover_only:
        return
    args.state_column, args.min_identity, args.min_coverage = "state", 1.0, 1.0
    args.min_fit_atoms, args.contact_cutoff = 20, 4.0
    original_loader = geometry.load_structure
    structure_hashes = {}

    @lru_cache(maxsize=24)
    def cached_structure(root, pdb, entity_ids):
        result = original_loader(root, pdb, entity_ids)
        structure_hashes[str(result[2])] = geometry.sha256(result[2])
        return result

    def loader(root, pdb, entity_ids):
        return cached_structure(root, pdb, frozenset(entity_ids))

    inventory, pairs, residues, errors = [], [], [], []
    geometry.load_structure = loader
    try:
        for i, uni in enumerate(sorted(specifications), 1):
            substrates, cofactors, ignored = specifications[uni]
            by_pdb = defaultdict(list)
            for row in rows_by_uni[uni]:
                # Only structures whose entity inventory has a supported cofactor.
                if geometry.ligand_codes(row["ligands"]) & cofactors:
                    by_pdb[row["pdb_id"]].append(row)
            selected = defaultdict(lambda: defaultdict(list))
            for pdb, entities in sorted(by_pdb.items()):
                try:
                    chains, ligands, _ = loader(args.mmcif_root, pdb, {r["entity_id"] for r in entities})
                    # Test all non-water component contacts, including ligands
                    # omitted from global apo/holo labels as artifacts.
                    ligand_codes = {code for code, _ in ligands} - {"HOH", "DOD"}
                    for label, chain in sorted(chains.items()):
                        _, identity, coverage = geometry.sequence_mapping(references[uni], chain["sequence"])
                        codes = {c for c in ligand_codes if geometry.contacts_ligand(chain, ligands, {c}, 4.0)}
                        state, retained, bound, reason = classify_contacts(codes, substrates, cofactors, ignored)
                        if identity < 1 or coverage < 1:
                            state, reason = "excluded", "mutation or incomplete reference sequence"
                        entity = next(r for r in entities if r["entity_id"] == chain["entity"])
                        inventory.append({"uniprot_id": uni, "pdb_id": pdb, "chain": label,
                                          "state": state, "reason": reason, "identity": identity,
                                          "coverage": coverage, "cofactors": ";".join(sorted(retained)),
                                          "substrates": ";".join(sorted(bound)), "contacts": ";".join(sorted(codes))})
                        if state != "excluded":
                            selected[tuple(sorted(retained))][state].append((pdb, label, entity, bound))
                except (OSError, ValueError, KeyError) as error:
                    errors.append({"uniprot_id": uni, "pdb_id": pdb, "reason": str(error)})
            for retained, states in sorted(selected.items()):
                free, bound = defaultdict(list), defaultdict(list)
                for pdb, label, entity, codes in states["substrate_free"]:
                    free[pdb].append((label, entity, codes))
                # Separate distinct substrate combinations in the same PDB.
                for pdb, label, entity, codes in states["substrate_bound"]:
                    bound[(pdb, tuple(sorted(codes)))].append((label, entity, codes))
                for f, (b, bound_codes) in product(sorted(free), sorted(bound)):
                    if f == b:
                        continue
                    fs, bs = free[f], bound[(b, bound_codes)]
                    group = {"uniprot_ids": uni, "seq_group": "_".join((uni, f, b, "+".join(retained), "+".join(bound_codes))),
                             "enzyme_name": fs[0][1]["enzyme_name"], "ec_numbers": fs[0][1]["ec_numbers"],
                             "best_apo": f, "best_holo": b, "best_apo_res_A": fs[0][1]["resolution_A"],
                             "best_holo_res_A": bs[0][1]["resolution_A"],
                             "apo_chains": ";".join(r[0] for r in fs), "holo_chains": ";".join(r[0] for r in bs)}
                    synthetic = defaultdict(list)
                    for pdb, entries, state in ((f, fs, "apo"), (b, bs, "holo")):
                        for entity_id in sorted({r[1]["entity_id"] for r in entries}):
                            row = dict(next(r[1] for r in entries if r[1]["entity_id"] == entity_id))
                            row.update(state=state, seq_group=group["seq_group"], ligands=";".join(bound_codes) if state == "holo" else "")
                            synthetic[pdb].append(row)
                    try:
                        pair, measured = geometry.compare_group(group, catalytic[uni], references[uni], synthetic, args)
                    except (OSError, ValueError, KeyError) as error:
                        pair, measured = {"uniprot_id": uni, "seq_group": group["seq_group"],
                                          "apo_pdb": f, "holo_pdb": b, "status": "error", "reason": str(error)}, []
                    pair.update(substrates=";".join(bound_codes), cofactors=";".join(retained), ec_numbers=group["ec_numbers"])
                    pairs.append(descriptive_fields(pair))
                    residues.extend(descriptive_fields(r) for r in measured)
            cached_structure.cache_clear()
            print("Compared {}/{} proteins: {} pairs".format(i, len(specifications), len(pairs)), flush=True)
    finally:
        geometry.load_structure = original_loader
    geometry.write_table(args.output_dir / "structure_inventory.tsv", inventory)
    geometry.write_table(args.output_dir / "pairs.tsv", pairs)
    geometry.write_table(args.output_dir / "proteins.tsv", protein_summary(pairs))
    geometry.write_table(args.output_dir / "residues.tsv", residues)
    geometry.write_table(args.output_dir / "errors.tsv", errors)
    pair_index = {r["seq_group"]: r for r in pairs}
    export = []
    for residue in residues:
        pair = pair_index[residue["seq_group"]]
        annotations = [r for r in catalytic[residue["uniprot_id"]] if int(r["residue_pos"]) == residue["uniprot_position"]]
        row = dict(pair, **residue)
        row.update(pair_status=pair["status"], residue_status=residue["status"],
                   catalytic_roles=";".join(sorted({r["roles"] for r in annotations})))
        export.append(row)
    fields = list(dict.fromkeys(k for r in export for k in r)) or ["status"]
    with open(args.output_dir / "all_catalytic_residues.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        writer.writerows(export)
    complete = [p for p in pairs if p["status"] == "complete"]
    summary = {"cofactor_candidate_proteins": len(candidate_unis), "ligand_candidate_proteins": len(specifications),
               "paired_proteins": len({r["uniprot_id"] for r in pairs}), "pairs": len(pairs),
               "status_counts": dict(Counter(r["status"] for r in pairs)),
               "complete_proteins": len({r["uniprot_id"] for r in complete}), "structure_errors": len(errors), "movement": {}}
    for metric in ("ca_displacement_A", "sidechain_rmsd_A", "heavy_atom_rmsd_A"):
        summary["movement"][metric] = {}
        for threshold in (0.5, 1.0, 2.0):
            moved = [r for r in complete if r.get("max_" + metric) is not None and r["max_" + metric] >= threshold]
            summary["movement"][metric][str(threshold)] = {"pairs": len(moved), "proteins": len({r["uniprot_id"] for r in moved})}
    manifest = {"parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                "script_sha256": geometry.sha256(__file__),
                "geometry_sha256": geometry.sha256(geometry.__file__), "structures": structure_hashes,
                "versions": {"python": platform.python_version(), "numpy": np.__version__,
                             "biopython": Bio.__version__, "scipy": scipy.__version__},
                "reaction_inputs": reaction_hashes,
                "inputs": {str(getattr(args, k)): geometry.sha256(getattr(args, k))
                           for k in ("entity_states", "catalytic_residues", "reference_fasta", "artifacts")},
                "annotations": {str(p): geometry.sha256(p) for p in cache.glob("*.json")}}
    for name, data in (("summary.json", summary), ("manifest.json", manifest)):
        (args.output_dir / name).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
