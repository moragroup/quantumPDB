# M-CSA catalytic-residue apo/holo movement

The offline script compares the **best-resolution apo and holo PDB entries per
input group**, intersects those groups with M-CSA by exact UniProt accession, and
writes auditable measurements. It requires Python with NumPy, Biopython and SciPy
(`python -m pip install numpy biopython scipy`). No network requests are made.

## Run on the current dataset

From the quantumPDB repository root:

```bash
python -m qp.analyze.catalytic_movement \
  --groups /mnt/labs/shared/enzymes/apo_holo/apo_holo_groups_uniprot.tsv \
  --entity-states /mnt/labs/shared/enzymes/apo_holo/entity_states.tsv \
  --mmcif-root /mnt/labs/shared/enzymes/mmCIF \
  --mcsa-annotations ../Fluxion/data/processed/mcsa_annotations.csv \
  --catalytic-residues ../Fluxion/data/processed/mcsa_catalytic_residues_full.csv \
  --reference-fasta ../Fluxion/data/processed/mcsa_catalytic_sequences.fasta \
  --output-dir results/catalytic_movement \
  --thresholds 0.5 1.0 2.0
```

For a small run, append `--limit 10` or `--uniprot P00698` (repeatable).
To compare identical-sequence groups, replace `--groups` with
`apo_holo_groups.tsv`. For organic-ligand-only labeling use
`apo_holo_groups_strict.tsv` **and** `--state-column state_strict`.
Group input must have a single exact UniProt accession per row; multi-accession
rows are excluded. The UniProt table may pair different constructs or mutants;
the identical-sequence table gives a more controlled sequence comparison.

## Method

1. Read reference catalytic positions and residue identities from M-CSA.
   Positions are **UniProt sequence positions**, not PDB author numbers.
2. Align each entity's complete mmCIF canonical polymer sequence to the supplied
   reference FASTA. Use a deterministic local alignment with unaligned termini
   (match +2, mismatch -1, internal gap open -10, extend -0.5). Default requirements:
   at least 95% identity among aligned positions and 50% reference coverage.
   Unresolved residues retain their polymer positions. Catalytic mutations,
   reference mismatches, and missing atoms are explicitly recorded.
3. Consider chains from the matching enzyme entity and state only. Require a
   holo chain to contact a ligand listed for that entity within 4 A. Contact is
   checked against deposited non-protein heavy atoms, without symmetry mates.
   The source labels are entity-level; the contact check avoids selecting an
   unbound copy of a holo enzyme.
4. Fit holo to apo by ordinary least-squares Kabsch alignment of matched,
   identical, resolved C-alpha atoms, excluding annotated catalytic positions.
   At least 20 fit atoms are required. Fit RMSD is reported; no outlier pruning
   or local active-site fit is performed.
5. Choose the chain pair with most fully resolved catalytic residues, then most
   fit atoms, then lowest fit RMSD, then lexical chain IDs. This does not choose
   pairs by the amount of catalytic movement. Only the first deposited model
   is used. Each atom uses its highest-occupancy conformer, with blank/A/lexical
   altloc tie breaks; conformational ensembles are not modeled.
6. For each annotated residue measure C-alpha displacement, side-chain heavy-atom
   RMSD, and whole-residue heavy-atom RMSD (backbone N/CA/C/O plus side chain).
   Glycine uses CA for the side-chain metric. Equivalent ASP/GLU/ARG/VAL/LEU
   labels and coupled PHE/TYR ring flips are minimized over. All required atoms
   must exist in both states for a full-residue comparison.

For each metric, a group is called **moved** at a threshold when **at least one**
catalytic residue reaches that threshold. Only groups with every annotated
catalytic residue fully resolved in both structures contribute to summary
movement counts. Partial measurements remain available in the tables. When
several sequence groups represent one UniProt, unique-UniProt counts mean any
fully mapped representative group satisfies the threshold. Catalytic positions
are the union across M-CSA entries for that accession, not separate counts of
M-CSA reactions or multisubunit active sites.

## Outputs

- `pairs.tsv`: accession, group, M-CSA IDs, selected PDBs/chains, sequence identity
  and coverage, backbone fit RMSD, mapping status, and maximum movement metrics.
- `residues.tsv`: each catalytic position, author numbering, missing/mutated
  status, and individual movement metrics. Chain IDs are mmCIF **label** IDs;
  author chain IDs are also recorded in `pairs.tsv`.
- `summary.json`: coverage and moved counts at every requested threshold.
- `all_catalytic_residues.csv`: every annotated catalytic residue for the selected
  representative apo/holo pairs, with unrounded displacement values and the same
  enzyme/EC/PDB metadata as the filtered CSV. Includes partial and unavailable
  comparisons; blank measurements mean unmeasurable, not zero movement.
  Enzymes without any catalytic annotation have no residue rows and remain
  listed in `pairs.tsv`.
- `moving_residues.csv`: one row per catalytic residue with side-chain RMSD
  at least 1 A in a complete pair, including enzyme name, source-group EC numbers,
  residue-specific M-CSA IDs/EC numbers/roles, UniProt and PDB author positions,
  apo/holo PDBs and chains, resolution, alignment RMSD, and movement metrics.
  Change this cutoff with `--csv-threshold`.
- `manifest.json`: command parameters, Python/NumPy/Biopython/SciPy versions,
  script hash, input-table/FASTA hashes, and hashes of available selected mmCIF files.

Reruns overwrite these six files in the chosen output directory. Different
inputs should use different directories. Preserve the manifest with results.

These are **apo/holo-associated structural differences**, not causal proof of
binding-induced movement. The source definition includes metal/cofactor binding
and ligands outside the active site. Crystal packing, construct differences,
resolution, and fitting across moving domains can contribute to displacement.
One representative pair cannot establish that an enzyme never moves, nor does
it enumerate all available apo/holo pairs. No significance estimate from atomic
coordinate uncertainty is assigned to the 0.5/1/2 A cutoffs.

## Verify the calculations

```bash
python -m pytest qp/tests/test_catalytic_movement.py -q
```
