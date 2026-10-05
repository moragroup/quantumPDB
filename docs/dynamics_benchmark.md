# QuantumPDB catalytic-motion benchmark

## Main result

QuantumPDB's current QM calculations **cannot sample catalytic heavy-atom
rearrangements**. `optimization: true` freezes all heavy atoms in TeraChem
(`qp/manager/create.py:find_heavy`) and uses `optimizehydrogens true` in ORCA
(`qp/manager/job_scripts.py:write_orca`). Single-point calculations do not move
atoms either. This is a capability assessment, not a completed simulation result.

The benchmark prepares both experimental endpoints for six proteins. Differences
between calculations started from different crystals are inherited input
differences; they are not evidence that the calculation predicted binding motion.
PrnA provides a small-motion comparison, but agreement at this scale alone does
not validate conformational sampling.

## Prepared cases

One complete pair per example is chosen to reproduce the previously reported
maximum catalytic sidechain RMSD. This retrospective selection is intentionally
a stress test, not an unbiased estimate of typical motion.

| Protein | Substrate-free → bound | Max sidechain RMSD (Å) |
|---|---|---:|
| P450cam, P00183 | 1GEK A → 3WRI B | 4.98 |
| p-Cresol methylhydroxylase, P09788 | 1DII B → 1DIQ A | 2.97 |
| Seryl-tRNA synthetase, Q46AN5 | 2CIM A → 2CJB B | 3.30 |
| Epi-aristolochene synthase, Q40577 | 5ILD A → 5IK0 A | 2.59 |
| Tryptophan synthase beta, P0A2K1 | 2RHG B → 8B06 B | 2.32 |
| PrnA, P95480 | 2APG A → 2AQJ A | 0.24 |

P450cam and seryl-tRNA synthetase also have appreciable global fit RMSDs
(1.94 and 2.67 Å). Their sidechain values include global/domain deformation.

## Setup and execution

```bash
python -m qp.analyze.dynamics_benchmark prepare
python -m qp.analyze.dynamics_benchmark run --stage run
python -m qp.analyze.dynamics_benchmark run --stage submit
```

`prepare` runs offline against `substrate_movement_results/all_mcsa/`. It writes
12 local first-model PDBs from the original mmCIFs, individual CSV/YAML inputs,
`cases.tsv`, `targets.tsv`, a SHA256 manifest, and an initial evaluation under
`dynamics_benchmark_results/`. Absolute paths make configurations independent
of the launch directory. Retained protein partners remain in the input.

Each cluster is centered on all annotated catalytic residues plus contacting
cofactors, using explicit author-chain/residue keys. Two interaction spheres,
ACE/NME caps, and no atom-count pruning prioritize target coverage. These can
be large models; inspect generated atom counts before choosing resources.
This all-target cluster setup differs from a single metal-centered cluster.

Preparation requires licensed Modeller in the running Python environment and
network access to Protoss. The current environment lacks Modeller, so the
attempted run stopped at dependency preflight and no QM jobs were launched.

Default job generation uses ORCA/B3LYP/def2-SVP, D3BJ, CPCM dielectric 10,
16 CPU processes, 64 GB, and the `cpu` Slurm partition. These are editable
starting settings. The generated script expects the `orca/6.1.1` environment
module. `/usr/bin/orca` on this host reports version 46.1 (the desktop screen
reader), and is **not the quantum-chemistry executable**. Configure the actual
ORCA module/path on the execution nodes.

P450cam and p-cresol methylhydroxylase CSV multiplicities are deliberately blank:
their heme redox/spin assignments must be specified from the intended state.
The job-generation wrapper rejects blank electronic parameters. Other examples
start with singlet multiplicity. `oxidation: 0` is an additive charge correction
to QuantumPDB's Protoss-derived ligand/residue charges; inspect `charge.csv`
and the generated total charge rather than adding metal charges twice. Covalent
PLP/FAD/heme links and the protonation states also need to be represented
correctly in the extracted models.

`--stage submit` invokes `qp submit` to create inputs. Configurations initially
set `submit_jobs: false`; enable that field when the settings are ready to launch.
Logs are saved per endpoint, and the wrapper catches CLI-reported failures even
when the CLI exits with status zero.

## Evaluate actual outputs

```bash
python -m qp.analyze.dynamics_benchmark evaluate
python -m qp.analyze.dynamics_benchmark evaluate --outputs final_outputs.csv
```

The optional output CSV has columns:

```csv
case_id,state,cluster_dir,final_xyz,converged
P95480,substrate_free,/absolute/path/to/cluster,/absolute/path/to/final.xyz,true
```

Choose one cluster per endpoint containing the targets. `final_xyz` must contain
one final frame with the same atom ordering as the numbered sphere PDBs used by
QuantumPDB's XYZ writer. Mark convergence from the QM output, not simply from
the existence of an XYZ file. Element/count mismatch is rejected; element checks
cannot detect swaps of atoms of the same element, so preserve ordering.

`evaluation.tsv` reports the observed target movement, fixed-atom protocol
capability, cluster sidechain coverage, and (when supplied) simulated sidechain
relaxation relative to that endpoint's starting cluster. Relaxation is fitted
on at least three non-collinear other CA atoms and accounts for symmetric
sidechain atom labels. It is **not** a predicted free-to-bound transition score.
Missing residues, incomplete sidechains, failed convergence, and insufficient
fit anchors are explicit statuses. Motion above 0.01 Å triggers a mapping/input
consistency note because the configured heavy atoms should be fixed.

To test recovery of the experimental transitions, a separate heavy-atom-mobile
protocol is needed, with substrate perturbation from a common starting geometry
and matched atom selections/boundaries. Optimizing both experimental endpoints
would test endpoint stability, not whether the transition was predicted.
Geometry minimization also does not establish thermal dynamics or sampling
probabilities.
