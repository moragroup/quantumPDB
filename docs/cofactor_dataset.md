# Dataset-wide substrate binding with cofactor retained

Run from the repository root:

```bash
python -m qp.analyze.cofactor_dataset \
  --entity-states /mnt/labs/shared/enzymes/apo_holo/entity_states.tsv \
  --mmcif-root /mnt/labs/shared/enzymes/mmCIF \
  --catalytic-residues ../Fluxion/data/processed/mcsa_catalytic_residues_full.csv \
  --reference-fasta ../Fluxion/data/processed/mcsa_catalytic_sequences.fasta \
  --mcsa-root ../Fluxion/data/raw/M-CSA/all_enzymes \
  --artifacts /mnt/labs/shared/enzymes/scripts/ligand_artifacts.txt \
  --output-dir substrate_movement_results/all_mcsa
```

The analysis scans the shared enzyme dataset for proteins with M-CSA catalytic
annotations and reference sequences, including proteins absent from the global
apo/holo groups. It downloads and caches PDBe component, ChEBI compound, and
UniProt annotations; reruns reuse `annotation_cache/`. `--discover-only` stops
before coordinate analysis; `--uniprot` optionally restricts accessions.

## Eligibility

- Substrates match M-CSA reaction **reactants** by exact ChEBI cross-reference
  or the first two InChIKey blocks (preserving stereochemistry, allowing
  protonation differences). Small compounds of molecular weight ≤50 are not
  substrate targets. Products do not qualify unless also annotated reactants.
- Recognized organic cofactors and iron–sulfur clusters are enumerated in
  `ORGANIC_COFACTORS`. Metal ions require an enzyme-specific UniProt cofactor
  annotation. ATP, SAM, and free CoA are not automatically called cofactors.
- Ligand heavy atoms must contact the selected protein chain within 4 Å.
  Each chain must cover the entire M-CSA reference at 100% identity, as in
  the PrnA comparison. Unresolved coordinates are recorded separately.
- The substrate-free chain has only supported cofactors after removal of
  listed artifacts. The bound chain has those cofactors plus mapped substrate(s).
  Other contacting ligands exclude the chain. Reaction-mapped substrates and
  products override the artifact list.
- Both chains must have the **identical set of cofactor component codes**.
  Oxidized/reduced variants are not interchangeable. Cofactor stoichiometry
  and protonation within a single component code are not controlled.
- All eligible **distinct PDB free/bound combinations** are compared. One
  representative chain pair is selected per combination/cofactor/substrate set,
  preferring complete catalytic residues, more fitted Cα atoms, then lower fit RMSD.

Candidate discovery deliberately overincludes entity-level ligand mixtures;
final eligibility is chain-specific. Superposition and symmetry-aware residue
metrics reuse `catalytic_movement.py`: matching noncatalytic Cα atoms, at least
20, with no outlier pruning. Movement thresholds count only complete pairs.

## Outputs and result

Outputs include `discovery.tsv`, `structure_inventory.tsv`, `pairs.tsv`,
`proteins.tsv`, `residues.tsv`, `all_catalytic_residues.csv`, `errors.tsv`,
`summary.json`, and `manifest.json`. The manifest hashes input tables,
M-CSA reactions, cached annotations, analysis code, and inspected structures.
CSV contains catalytic positions, roles, EC numbers, both PDB/chain IDs,
resolutions, fit RMSD, and three displacement metrics.

The current run found 676 cofactor-bearing candidate proteins, 106 requiring
coordinate checks, and **19 proteins with 74 eligible pairs**: 69 complete pairs
covering 18 proteins, and 5 partial pairs. No structure-loading errors occurred.

| Threshold (Å) | Proteins with catalytic side-chain movement | Complete pairs with movement |
|---|---:|---:|
| ≥0.5 | 12 / 18 | 36 / 69 |
| ≥1.0 | 11 / 18 | 21 / 69 |
| ≥2.0 | 6 / 18 | 10 / 69 |

These are conservative **identified matches**, not an exhaustive census of all
cofactor enzymes: proteins without M-CSA mappings, unmapped substrate chemistry,
polymeric or covalently attached substrates, unlisted cofactors, mutants, and
truncated constructs are outside eligibility. The absence of modeled substrate
does not prove experimental absence. Comparisons are static structural
differences; coordinate uncertainty, crystal packing, and global conformational
changes can contribute. Inspect fit RMSD before interpreting large movements.
