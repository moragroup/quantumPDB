# Plan: Fluxion-guided QM reaction searches with QuantumPDB

## Goal

Use Fluxion's reaction hypotheses to reduce the cost of finding and validating
enzyme reaction pathways. QuantumPDB supplies prepared active-site structures;
a new bridge translates atom-mapped reaction hypotheses into targeted QM
calculations. Validated intermediates, transition states, barriers, and failures
are returned as structured evidence for comparing Fluxion hypotheses.

The initial bridge is implemented in `qp/reaction/`, with the Fluxion exporter
in `../Fluxion/scripts/export_qm_hypotheses.py`. See
[`reaction_search.rst`](reaction_search.rst) for the implemented contract,
commands, mapping recipe, validation scope and real ORCA pilot results. This
document retains the broader design and enzyme-benchmark milestones.

## Current capabilities and important boundaries

- Fluxion's `../Fluxion/MODEL_FLOW.md` describes predictions in bond-electron
  (BE) matrices: reaction edits and graph endpoints. These are not, by
  themselves, Cartesian transition-state geometries or physical reaction paths.
- Its start–end models can predict a net substrate-to-product transformation.
  As documented in `../Fluxion/startend_implementation_summary.md`, this can
  collapse multiple mechanism steps. Such predictions must be distinguished
  from elementary-step proposals before a transition-state search.
- Flow interpolation time and fractional BE entries must not be interpreted
  directly as physical time, molecular distances, or a transition state.
- QuantumPDB prepares structures, assigns protonation through Protoss, extracts
  capped QM clusters, and generates TeraChem/ORCA jobs.
- The existing cluster-optimization mode freezes all heavy atoms: TeraChem constraints
  originate in `qp/manager/create.py:find_heavy`; ORCA uses
  `optimizehydrogens true` in `qp/manager/job_scripts.py:write_orca`.
  `qp.reaction` provides explicit movable atoms and frozen boundaries for reaction
  searches, with full-cluster IRC or constrained downhill connectivity validation.
- The six existing examples in `dynamics_benchmark_results/` are structural
  endpoint benchmarks. Catalytic sidechain displacement is not a reaction
  barrier or evidence of a transition state.

## Proposed architecture

```text
Fluxion prediction + checkpoint provenance
                 |
                 v
Validated reaction-hypothesis record
  atom identities, bond edits, electron/proton changes, candidate rank
                 |
                 v
Structure mapping + QuantumPDB active-site preparation
  substrate/cofactor/partners, waters, residues, boundaries
                 |
                 v
3D candidates + inexpensive geometry/path screening
                 |
                 v
DFT endpoint optimization -> TS/path search -> validation
                 |
                 v
Structured results + cost accounting -> hypothesis comparison
```

Keep this as an adapter initially rather than making either repository depend
on the other's training/runtime stack. Use a versioned JSON interchange format
and ordinary structure files. Preserve existing QuantumPDB calculation modes.

## Phase 1 — Define the interchange contract

### Tasks

1. Inspect the actual inference outputs and checkpoint used for the pilot.
   Establish whether it is stepwise or start–end and whether stable atom maps,
   explicit hydrogens, radicals, and intermediate graphs are available.
2. Write an exporter for candidate hypotheses. Retain raw predictions and their
   ranking, then derive graph changes with an explicitly recorded conversion
   policy. Do not silently repair invalid chemistry or turn model scores into
   probabilities unless calibrated.
3. Introduce a schema with these required groups:

| Group | Required information |
|---|---|
| Identity | schema version, reaction ID, hypothesis ID, elementary step ID |
| Provenance | checkpoint/config hashes, sampling settings, raw-output hash |
| Chemical atoms | stable map ID, element, molecular participant, formal charge, radical state |
| Endpoints | reactant/product graphs; explicit hydrogen and atom-mapping policy |
| Reaction edits | forming/breaking bonds, bond-order changes, proton donor/acceptor atom maps |
| Electronic state | total charge, multiplicity, cofactor redox state; alternatives when unresolved |
| Structure references | PDB/mmCIF, chain, ligand instance, geometry file where available |
| Structure mapping | chemical map ID to structure atom key and QM index |
| Interpretation | net transformation or elementary step; candidate rank/score semantics |

4. Permit optional proposed intermediate/TS XYZ geometries and reaction
   coordinates. An absent 3D guess should be explicit, not synthesized from BE
   interpolation without geometry generation.
5. Validate element/atom conservation, electron accounting, charge consistency,
   stereochemistry, map uniqueness, and required participants. For proton or
   electron exchange, represent the donor/acceptor or reservoir explicitly in
   the model bookkeeping.

### Deliverables and acceptance

- Schema, example JSON records, and Fluxion exporter.
- Round-trip tests preserve atom IDs and distinguish a net reaction from a step.
- Invalid/ambiguous hypotheses produce explicit reasons, not QM jobs.
- Mapping tests include symmetry, hydrogen transfer, and changed bond orders.

## Phase 2 — Map hypotheses into prepared active sites

### Tasks

1. Select the exact substrate-bound structure, biological partners, cofactor
   instance, relevant reactants, and waters. Store both author and label chain/
   residue identifiers to avoid ambiguity.
2. Map graph atoms to structure atoms using ligand topology and stereochemistry,
   not atom names alone. Enumerate chemically equivalent mappings where needed.
3. Use QuantumPDB preparation, then verify that remodeling/protonation has
   retained every reaction participant and its intended chemistry.
4. Include substrate, reactive cofactor region, participating waters, and all
   residues involved in proton/electron transfer or metal coordination in QM.
   If catalytic annotations are unavailable, use pocket proximity and chemical
   roles to propose the initial region; missing annotations do not justify
   omitting the protein environment.
5. Define boundary atoms away from the reaction center. Freeze suitable boundary
   anchors while permitting the reactive atoms and relevant sidechains to move.
   Preserve covalent ligand–protein links; do not cap through a bond that changes
   in the proposed step.
6. Generate a persistent atom-index table connecting chemical map IDs,
   prepared structure keys, capped cluster atoms, and QM indices. Document the
   different index conventions expected by each QM engine.
7. Build separate electronic/protonation alternatives where the chemistry is
   uncertain. Verify total charge and spin after cluster extraction.

### Deliverables and acceptance

- Prepared model, atom map, boundary/movable atom sets, and preparation report.
- All reacting atoms have unambiguous mappings and are movable.
- Reactant and product models have identical atoms, boundaries, and an
  internally consistent environmental treatment for energy comparison.
- Tests cover capped clusters, covalent cofactors, and ligand/residue numbering.

## Phase 3 — Add targeted reaction-search job modes

Implement ORCA first, using the actual installed quantum-chemistry version and
its documented syntax. TeraChem support can follow after the workflow works.

### New calculation modes

| Mode | Purpose |
|---|---|
| `endpoint_opt` | Relax reactant/product/intermediate geometries with boundary constraints |
| `constrained_scan` | Explore selected forming/breaking distances or coupled proton-transfer coordinates |
| `path_search` | Find a path between mapped endpoint geometries, e.g. NEB with subsequent TS refinement |
| `ts_opt` | Refine a saddle-point guess using suitable curvature/Hessian information |
| `frequency` | Establish minima or reaction-specific saddle character |
| `irc` | Check the connections to the intended neighboring intermediates |

Treat these as separate from the existing boolean `optimization` option.
Record the exact engine input, constraints, model, and method for every job.

### Search decision tree

1. **Usable 3D TS guess:** validate structure/electronic state, then attempt
   targeted TS refinement. Do not first minimize the guess unconstrained,
   which would normally collapse it into a minimum.
2. **Mapped endpoint geometries:** relax both consistently, generate a path,
   and refine promising saddle candidates.
3. **Graph edits only:** generate plausible endpoint geometries and orientation
   candidates; use a small set of constrained scans or path searches to seed TS
   refinement.
4. **Net multistep transformation:** obtain stepwise hypotheses or explicitly
   propose intermediate graphs before assigning one TS to the transformation.

### Cost reduction

- Screen geometry/path candidates with a lower-cost method suitable for the
  cofactor chemistry; verify it against a small DFT reference subset.
- Prioritize a bounded number of Fluxion candidates, while keeping an alternative
  path so a wrong first hypothesis does not terminate discovery.
- Reuse wavefunctions and Hessians only between compatible atom orderings,
  electronic states, and engine settings. Check state tracking for metals.
- Cache identical endpoint calculations using model/geometry/settings hashes.
- Schedule independent candidates concurrently within an explicit CPU/GPU-hour
  budget. Stop and report failed searches rather than endlessly restarting them.

### Acceptance

- A synthetic known reaction exercises every new mode end to end.
- Atom indexing and frozen/movable selections are checked in generated inputs.
- The existing hydrogen-only and single-point modes retain their behavior.
- Restart logic never reuses files across incompatible charge/spin states.

## Phase 4 — Validate chemistry and return results to Fluxion

For each proposed elementary step:

1. Confirm optimizer and electronic-structure convergence.
2. Verify endpoint minima in the movable degrees of freedom.
3. Check that the candidate TS has one significant imaginary mode in the
   evaluated movable subspace and that its displacement matches the proposed
   reaction. Boundary-constrained or partial-Hessian results must be labeled;
   they are not full-system vibrational validation.
4. Run an IRC where suitable and supported, or a documented equivalent
   connectivity check, and optimize the connected endpoints. Confirm graph and
   atom-mapping agreement with the proposed neighboring intermediates.
5. Compare energies using the same QM atoms, environment, method, and electronic
   state conventions. Report electronic barriers separately from corrected
   free-energy estimates, stating temperature, standard state, and treatment of
   low-frequency modes. Cluster corrections are not full enzyme free energies.
6. Check sensitivity to important protonation/spin states and a larger QM region
   before interpreting close competing barriers.

Return structured records containing:

- Validated geometries and atom maps.
- Endpoint/TS energies, barriers, units, reference states, and corrections.
- Imaginary frequencies and reaction-mode interpretation.
- Connectivity evidence and validation status.
- Charge, spin, constraints, region, and computational methods.
- Failure reason, iterations, elapsed time, allocation, and resource cost.

Use statuses such as `validated_step`, `wrong_connectivity`, `not_converged`,
`invalid_mapping`, and `electronic_state_unresolved`. A failed search is not
proof that a mechanism is impossible. A lower computed barrier is not proof of
the enzyme's dominant mechanism without model/environment sensitivity checks.

## Phase 5 — Pilot and measure the speedup

### Pilot selection

Start with one well-defined elementary step with a reference mechanism and
usable structure. Prefer manageable non-metal chemistry over P450cam or
multicofactor redox chemistry for the first integration test.

Tryptophan synthase beta (`P0A2K1`) is a candidate from the existing examples,
but select it only after checking the exact PLP covalent state, proposed step,
participants, and atom maps. If those are ambiguous, choose a simpler M-CSA
example. PrnA remains a useful preparation/mapping example, but its small
structural displacement does not make its halogenation chemistry a simple pilot.

### Matched comparison

Run the same chemical step with the same model, methods, compute budget, and
validation criteria under:

1. A conventional endpoint-driven search without Fluxion-specific reaction
   guidance.
2. A Fluxion-guided search using graph edits and candidate intermediates.
3. Where available, a search initialized from a proposed 3D TS geometry.

Define the baseline's available chemical information before running it; do not
remove essential participants to manufacture a speedup. Keep reference TS
coordinates out of the prediction inputs and disclose checkpoint training
overlap with M-CSA pilot mechanisms.

Measure:

- Fraction of searches finding the intended validated saddle/path.
- Number of expensive energy/gradient/Hessian evaluations.
- Total CPU/GPU-hours, including preparation and failed candidates.
- Wall time with queue time reported separately.
- Barrier agreement and sensitivity to model changes.
- Performance across multiple initial guesses, including unsuccessful trials.

Report savings in expensive evaluations and resource cost first. Wall-time
speedup also depends on scheduling and parallelism. Do not claim a universal
speedup from a single successful case.

### Pilot success criteria

- At least one fully mapped, validated elementary step is reproducible.
- Guided searches recover the same intended connectivity as the baseline.
- Comparisons account for failed attempts and all screening costs.
- An interpretable improvement in search cost or success rate is demonstrated;
  if not, identify whether mapping, geometry generation, or the hypothesis is
  the bottleneck before scaling.

## Proposed implementation files

These are proposed names, not existing interfaces:

```text
QuantumPDB:
  qp/reaction/schema.py            # interchange validation
  qp/reaction/mapping.py           # graph -> prepared structure -> QM atom IDs
  qp/reaction/geometry.py          # endpoint/TS candidate generation
  qp/reaction/jobs.py              # engine reaction-search modes
  qp/reaction/workflow.py          # scheduling, restarts, budgets
  qp/reaction/results.py           # parsing, validation records, cost reporting
  qp/tests/test_reaction_*.py      # mapping, constraints, parsers, integration

Fluxion:
  scripts/export_qm_hypotheses.py   # checkpoint output -> interchange records

Per reaction:
  hypothesis.json
  atom_mapping.tsv
  model/                          # structures, region, constraints
  candidates/<id>/                # endpoints, scans/paths, TS, validation
  results.json
  cost_report.tsv
  manifest.json
```

Store engine version, package versions, input hashes, checkpoint identity, and
raw output references in the manifest. Initial feedback should be a results
table rather than automatic model retraining; decide on supervision only once
validation quality and bias are understood.

## Estimated effort and compute budget

Planning estimates, not measured timings:

| Work | Approximate effort |
|---|---|
| Export/schema and pilot chemistry audit | 2–4 working days |
| Structure mapping and movable-region implementation | 3–7 working days |
| One-engine job modes, validation parsing, and orchestration | 1–2 weeks |
| First elementary-step pilot and matched cost comparison | 1–3 weeks, overlapping development |

Overall, allow roughly **3–6 weeks for a credible first integrated benchmark**,
depending on mapping quality and chemistry. A tractable individual guided step
may take hours to days of compute; several guesses, frequency/IRC validation,
and electronic-state alternatives can extend it to days or weeks. Metal/redox
systems and multistep mechanisms can take substantially longer.

Before scaling, time representative energy/gradient and Hessian calculations
on the actual cluster and hardware. Use those measurements to set per-candidate
budgets rather than extrapolating solely from atom count.

## Immediate next actions

1. Select one elementary-step Fluxion output and its checkpoint/config.
2. Audit whether it contains stable atom maps and stepwise intermediates.
3. Agree on the pilot's protonation, cofactors, reactants, and electronic state.
4. Export a hypothesis record and map it into one QuantumPDB active-site model.
5. Implement heavy-atom-mobile endpoint and one targeted TS/path-search route.
6. Validate that step, then run the matched search-cost comparison.

The central expected benefit is **fewer expensive searches**, not cheaper DFT
evaluations. Catalytic residues can remain absent from Fluxion's symbolic
prediction only if the bridge explicitly adds chemically required participants
to the QM model; excluding them from the physical model is a separate and often
unreliable approximation.
