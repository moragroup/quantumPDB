# Substrate-binding comparison with cofactor retained

`qp.analyze.substrate_movement` compares substrate-free and substrate-bound
structures independently of their global apo/holo labels. Ligand component
identities are supplied explicitly. Classification uses non-polymer ligand
heavy atoms within 4 Å of each protein chain; polymer tryptophan residues do
not count as bound substrate.

The script requires full reference-sequence coverage and 100% sequence
identity, a contacting cofactor of the specified chemical identity, and no
contacting product/analogue from the supplied exclusion list. It compares all
eligible distinct substrate-free/bound PDB pairs, selecting one chain pair
per PDB pair with the same ranking and geometry as `catalytic_movement.py`.
Chain eligibility is checked separately on both sides. A substrate-free
classification means no modeled contacting substrate or specified analogue;
it does not establish the experimental solution composition.

## PrnA command

Run from the repository root:

```bash
python -m qp.analyze.substrate_movement \
  --uniprot P95480 --substrate TRP --cofactor FAD --exclude-ligand CTE \
  --entity-states /mnt/labs/shared/enzymes/apo_holo/entity_states.tsv \
  --mmcif-root /mnt/labs/shared/enzymes/mmCIF \
  --catalytic-residues ../Fluxion/data/processed/mcsa_catalytic_residues_full.csv \
  --reference-fasta ../Fluxion/data/processed/mcsa_catalytic_sequences.fasta \
  --output-dir substrate_movement_results/prna
```

Outputs: `structure_inventory.tsv` (including exclusion reasons), `pairs.tsv`,
`residues.tsv`, `all_catalytic_residues.csv`, `summary.json`, and `manifest.json`.
The manifest hashes the input tables, both analysis scripts, and every
examined structure, and records parameters and software versions. Output
columns use `substrate_free_` / `substrate_bound_` prefixes.

## PrnA result

One eligible wild-type pair: **2APG chain A → 2AQJ chain A**, at 1.9 and 1.8 Å
resolution. Both retain FAD; 2AQJ additionally contains tryptophan. The fit
uses 515 matching resolved noncatalytic Cα atoms, with RMSD **0.242 Å**.
Catalytic atoms are excluded from the fit; side-chain RMSD accounts for
chemically equivalent atom labels.

| Catalytic residue | Cα displacement (Å) | Side-chain RMSD (Å) | All-heavy-atom RMSD (Å) |
|---|---:|---:|---:|
| Lys79 | 0.199 | 0.137 | 0.182 |
| Glu346 | 0.172 | 0.240 | 0.204 |

Neither catalytic residue exceeds the earlier 0.5, 1, or 2 Å movement
thresholds. This pair therefore shows little catalytic-residue rearrangement
on substrate binding. Sub-Å differences can include coordinate uncertainty;
these results do not establish time-dependent dynamics or a general result
for the halogenase family.

Excluded structures:

- **2AR8:** bound CTE (7-chlorotryptophan product).
- **2ARD:** FDA (reduced flavin), rather than matching FAD.
- **2JKC:** catalytic E346D mutation, despite bound tryptophan.
- **4Z43:** E450K mutation.
- **4Z44:** F454K mutation.
