"""Graph-based endpoint generation; BE interpolation is never Cartesian input."""

import numpy as np

from qp.reaction.mapping import read_xyz, write_xyz
from qp.reaction.schema import HypothesisError, bonds


def rdkit_molecule(record, endpoint):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise HypothesisError("Graph-only geometry generation requires RDKit; alternatively supply mapped endpoint XYZ files") from exc
    states = {a["map_id"]: a for a in record["endpoints"][endpoint]["atom_states"]}
    indices = {}
    mol = Chem.RWMol()
    for atom in record["atoms"]:
        state = states[atom["map_id"]]
        a = Chem.Atom(atom["element"])
        a.SetAtomMapNum(atom["map_id"])
        a.SetFormalCharge(state["formal_charge"])
        a.SetNumRadicalElectrons(state["radical_electrons"])
        a.SetNoImplicit(True)
        indices[atom["map_id"]] = mol.AddAtom(a)
    types = {1.0: Chem.BondType.SINGLE, 2.0: Chem.BondType.DOUBLE, 3.0: Chem.BondType.TRIPLE,
             1.5: Chem.BondType.AROMATIC}
    for (a, b), order in bonds(record["endpoints"][endpoint]).items():
        if order not in types:
            raise HypothesisError("Fractional/radical bonds need manually audited Cartesian geometries")
        mol.AddBond(indices[a], indices[b], types[order])
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    # Stereo must be explicitly realized, rather than lost during generic embedding.
    if any(a.get("stereo") for a in states.values()) or any(b.get("stereo") for b in record["endpoints"][endpoint]["bonds"]):
        raise HypothesisError("Stereo-specified graphs require mapped endpoint geometries or an audited stereochemical embedding")
    return mol


def generate_product(record, model, reactant_xyz, seed=2026):
    """Embed changed topology near reactant coordinates, keeping environment fixed.

    Distance-geometry bounds encode product bonds, then a molecular force field
    regularizes the participant geometry. This is an initial guess, not a QM
    stationary point. Covalent region/environment bonds require supplied product
    geometry; they cannot be silently severed by this isolated-graph generator.
    """
    if model.get("covalent_environment_links"):
        raise HypothesisError("Supply product XYZ for covalently linked QM/environment models")
    if record["provenance"].get("stereochemistry_audit_required"):
        raise HypothesisError("Stereo-bearing BE predictions require an audited mapped product geometry")
    mol = rdkit_molecule(record, "product")
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D
    rows = {r["map_id"]: r["qm_index"] for r in model["atom_mapping"]}
    indices = [rows[a["map_id"]] for a in record["atoms"]]
    result = reactant_xyz.copy()
    from rdkit import Chem
    fragment_indices = []
    fragments = Chem.GetMolFrags(mol, asMols=True, fragsMolAtomMapping=fragment_indices)
    # Embed/align each product fragment separately. Embedding disconnected
    # fragments together can place an isolated proton directly inside a bond.
    for fragment, global_indices in zip(fragments, fragment_indices):
        qm_indices = [indices[i] for i in global_indices]
        anchors = {i: Point3D(*map(float, reactant_xyz[q])) for i, q in enumerate(qm_indices)
                   if q in model.get("frozen_indices", [])}
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.useRandomCoords = True
        if anchors:
            params.SetCoordMap(anchors)
        if AllChem.EmbedMolecule(fragment, params) != 0:
            raise HypothesisError("Product distance-geometry embedding failed")
        ff = AllChem.UFFGetMoleculeForceField(fragment) if AllChem.UFFHasAllMoleculeParams(fragment) else None
        if ff is not None:
            for i in anchors:
                ff.AddFixedPoint(i)
            ff.Minimize(maxIts=500)
        product = np.asarray(fragment.GetConformer().GetPositions(), dtype=float)
        reference = reactant_xyz[qm_indices]
        p, r = product.mean(axis=0), reference.mean(axis=0)
        u, _, vt = np.linalg.svd((product - p).T @ (reference - r))
        correction = np.eye(3)
        correction[-1, -1] = np.linalg.det(u @ vt)
        result[qm_indices] = (product - p) @ (u @ correction @ vt) + r
    for i in model.get("frozen_indices", []):
        result[i] = reactant_xyz[i]
    return result


def prepare_geometries(record, model, source_dir, output_dir):
    from pathlib import Path
    source_dir, output_dir = Path(source_dir), Path(output_dir)
    elements, reactant = read_xyz(source_dir / model["reactant_xyz"])
    if elements != model["elements"]:
        raise HypothesisError("Reactant XYZ ordering differs from model inventory")
    if model.get("product_xyz"):
        other, product = read_xyz(source_dir / model["product_xyz"])
        if other != elements:
            raise HypothesisError("Endpoint atom inventory/order differs")
    else:
        product = generate_product(record, model, reactant, model.get("geometry_seed", 2026))
    frozen = model.get("frozen_indices", [])
    if frozen and not np.allclose(reactant[frozen], product[frozen], atol=1e-5, rtol=0):
        raise HypothesisError("Endpoint boundary coordinates must match")
    for side, xyz in (("reactant", reactant), ("product", product)):
        write_xyz(output_dir / (side + ".xyz"), elements, xyz)
    if model.get("ts_xyz"):
        other, ts = read_xyz(source_dir / model["ts_xyz"])
        if other != elements or (frozen and not np.allclose(ts[frozen], reactant[frozen], atol=1e-5, rtol=0)):
            raise HypothesisError("TS guess atom order/boundaries differ")
        write_xyz(output_dir / "ts_guess.xyz", elements, ts)
