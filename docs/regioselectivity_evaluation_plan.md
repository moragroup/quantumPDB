# Progesterone CYP regioselectivity: Fluxion → quantumPDB evaluation

## Recommendation and present status

Start with **human CYP21A2 (P08686), progesterone C21 hydroxylation**, using
the substrate-bound structure **4Y8W**. Then add **human CYP17A1 (P05093),
progesterone C17 hydroxylation**, as the contrasting enzyme. The same substrate
and two experimentally observed sites make this more informative than a single
correctly predicted reaction.

The offline benchmark/structure preflight has been run. **Enzyme QM calculations
have not been run.** The existing reaction bridge passed a small synthetic TS
pilot, but the CYP metal chemistry needs an additional audited adapter before
it can be used here. This document is the execution plan for that evaluation.

## 1. Fixed cases and labels

| Enzyme | Source row, zero-based | Source ID / reaction ID | Observed hydroxylation | Starting structure |
|---|---:|---|---|---|
| CYP21A2 / P08686 | 337 | r506 / 1000337 | C21 | 4Y8W, chain A; B/C sensitivity checks |
| CYP17A1 / P05093 | 333 | r363 / 1000333 | C17, 17α product | 4NKX, **A105L mutant**, requiring a WT model |

Both cases belong to the **publication-packaged substrate-held-out test split**:

- Source: `../Fluxion/CYPs/CYP_benchmark_regioselectivity_train.tsv`.
- Split: `../Fluxion/data/publication/reproducibility/Fig3/input/data/splits/cyp_subholdout/split.json`.
- Packaged split sizes: 1,805 training, 200 validation, 472 test rows.
- Progesterone substrate key: `RJKFOVLPORLFTN`.

Use this packaged split consistently. The current working split has different
training/validation sizes; the older enzyme–substrate-pair split is a separate
benchmark. Audit the exact checkpoint training manifest, including M-CSA
warm-start overlap, before describing the run as held-out inference.

The labels identify observed products, not quantitative product ratios or proof
that every other site is chemically impossible. This case was selected after
inspecting predictions, so it is a **retrospective diagnostic**, not an unbiased
estimate of benchmark-wide improvement.

## 2. Completed preflight

Reproduce with Python containing RDKit, NumPy and Biopython:

```bash
python -m qp.analyze.regioselectivity_preflight \
  --fluxion-root ../Fluxion \
  --mmcif-root /mnt/labs/shared/enzymes/mmCIF \
  --output-dir regioselectivity_evaluation_results/preflight
```

The report `regioselectivity_evaluation_results/preflight/preflight.json` stores
benchmark rows, held-out membership, stereochemistry-aware substrate/CCD atom
mapping, crystal distances, all 200 selected archived prediction rows, RDKit
version, and input/script SHA-256 hashes. Canonical ranks are computed **after
removing atom-map numbers**, matching the benchmark evaluator; they are not
biochemical carbon numbering or QM indices.

| Archived canonical site rank | Progesterone CCD carbon |
|---:|---|
| 0 | C21 |
| 18 | C17 |
| 7 | C6 |
| 9 | C16 |
| 10 | C15 |

For archived SRS `best_acc` runs s00–s09, each with 2,048 draws:

| Enzyme | Observed-site top-1, out of 10 runs | Mean C21 draw fraction | Mean C17 draw fraction |
|---|---:|---:|---:|
| CYP21A2 | 8 | 0.6090 | 0.0740 |
| CYP17A1 | 0 | 0.8531 | 0.0755 |

These are archived sampling frequencies, not calibrated probabilities. For
`best_site_auc`, top-1 counts are 7/10 and 0/10 respectively. Fix checkpoint
selection before new calculations; do not choose a checkpoint by this pair's
outcome. The archive names `cyp_reps/srs/s00/...` checkpoint paths, whose actual
weights have not been located. Do not substitute an enzyme-held-out checkpoint
and call it the same experiment.

### Structural audit

- **4Y8W**: 2.64 Å, progesterone `STR`, heme `HEM`; expression tags/truncation
  require sequence mapping, but no engineered mutation is annotated. Chain A
  has STR author residue 604, HEM 603 and axial Cys 429.
- Chain A Fe–C21 distance is approximately **4.03 Å**, versus **5.42 Å** for
  Fe–C17. Across chains A/B/C: C21 4.03–4.11 Å; C17 5.42–5.60 Å.
- **4NKX**: 2.794 Å, four progesterone-bound chains; the protein is **A105L**.
  The mutation corresponds to author residue 87, not author residue 105.
  Chain A STR is 601, HEM 600, axial Cys 442. Fe–C17 spans 4.72–4.83 Å;
  Fe–C21 spans 4.32–4.94 Å.

These are ferric crystal Fe–carbon distances, **not Compound I oxo–H distances,
activation barriers, or predictions of turnover**. Distance alone would not
resolve the CYP17A1 contrast.

Structure references: Pallan et al. (2015), DOI
`10.1074/jbc.M115.646307` (4Y8W); Petrunak et al. (2014), DOI
`10.1074/jbc.M114.610998` (4NKX). Local sources are
`/mnt/labs/shared/enzymes/mmCIF/y8/4y8w.cif.gz` and
`/mnt/labs/shared/enzymes/mmCIF/nk/4nkx.cif.gz`.

## 3. Question and competing hypotheses

Primary question: **Can enzyme-specific QM ranking distinguish C21 from C17
hydrogen abstraction, and can Fluxion guidance reduce the cost of finding
validated reaction paths?** Treat ranking quality and search speed as separate
outcomes.

Model the elementary **Compound I hydrogen-atom transfer (HAT)** step:

1. Break the selected progesterone C–H bond.
2. Form oxo O–H, yielding Compound II plus the substrate radical.
3. Preserve all atoms, overall charge and the chosen total spin surface.

Do not feed the net hydroxylated product to this step's TS search. Oxygen rebound
is a second step; add it only after HAT is validated, especially if HAT ordering
does not explain the observed product or C17 stereochemistry.

Pilot competitors are C21 and C17. Resolve distinct transferable hydrogens at
C21 instead of treating its methyl hydrogens as one geometry. For the full site
comparison add C6, C16 and C15, covering all five archived candidate sites.
Record carbon identity, individual H identity, approach orientation and pose.

## 4. Required implementation and model preparation

### Gate A — CYP-specific hypothesis conversion

The generic `../Fluxion/scripts/export_qm_hypotheses.py` deliberately rejects
unaudited metal BE conversion. Implement and test a CYP adapter before jobs:

- Reuse Fluxion's M-CSA 699-derived Compound I/HAT template and CYP waypoint
  builder; export HAT predictions before canonical SMILES loses atom identity.
- Preserve substrate, oxo, heme and thiolate stable maps. Convert BE electron
  bookkeeping into an explicitly audited physical graph; heme coordination and
  Fe redox cannot be inferred from ordinary organic valence rules.
- Supply mapped Cartesian reactant/product geometries. Do not ask generic
  RDKit embedding to reconstruct the Fe–porphyrin/thiolate complex.
- Declare porphyrin, axial thiolate, oxo and substrate protonation/charges,
  cluster charge, electron count and multiplicity. Do not equate per-atom BE
  radical parity with the physical global spin state.
- Tests must establish atom/electron conservation, correct C–H/O–H edits,
  stable map-to-QM indices, identical endpoint composition, and refusal of
  ambiguous metal states. Record any manually supplied alternatives distinctly
  from Fluxion-generated hypotheses.

### Gate B — enzyme geometry and electronic state

1. Begin with CYP21A2 chain A. Align the expression construct to UniProt and
   retain persistent author/label/QM mappings.
2. Prepare progesterone, complete heme, axial Cys/thiolate and nearby polar
   residues/waters. Inspect protonation and cluster boundary charges. Estimate
   150–250 atoms as a starting target, then use actual extraction counts;
   complete heme and substrate take priority over an arbitrary atom limit.
3. Place the Compound I oxo on the distal Fe side and optimize the reactant.
   Retain chemically appropriate Fe–O/Fe–S/porphyrin coordination. Review
   broken-symmetry solutions, spin populations and ⟨S²⟩ for **doublet and
   quartet** surfaces rather than selecting a singlet default.
4. Freeze only defined boundary atoms. Keep reacting C/H/O, Fe, coordination
   partners and local catalytic heavy atoms movable. Existing hydrogen-only
   cluster optimizations are unsuitable for this search.
5. For CYP17A1, first seek a WT progesterone-bound model; otherwise repair
   4NKX A105L → WT using sequence-aligned author L87A, sample the repaired
   pocket, and relax it. Preserve the unmodified mutant as a separate structural
   control. Report WT-model uncertainty; never label the mutant as WT.
6. Later repeat the winning/runner-up sites with alternative chains/poses,
   a larger cluster and dielectric sensitivity. Steroid orientation must not be
   constrained into the observed answer.

Gate acceptance: reproducible mapping, intact coordination, stable reactant
minimum and reviewed charge/spin state. If these fail, report model preparation
failure rather than a regioselectivity result.

## 5. Calculation design and budget

Use ORCA 6.1.1 through `python -m qp.reaction`; retain every failed search and
its cost. Starting protocol: unrestricted B3LYP-D3BJ/def2-SVP optimization with
tight SCF, an explicitly selected dielectric (initial ε = 10), then consistent
def2-TZVP single-point refinement on validated stationary points. Inspect
spin/state stability and compare a second functional if site barriers are close.
Charge, spin, cluster, boundaries and energy conventions must match within each
comparison. xTB-only ranking is not accepted as validation of heme redox chemistry.

### Stage 1: smallest informative pilot

- CYP21A2 chain A only; C21 versus C17; doublet and quartet.
- Four site/spin combinations, each with baseline and guided search: **eight
  route-level attempts initially**, with distinct H choices expanded as needed.
- Cache shared reactant optimizations only under exact compatibility. Use
  coupled C–H breaking/O–H forming scans for guidance, followed by TS refinement;
  compare with an endpoint-based NEB-TS baseline on the same chemical problem.
- Initial search allocation per attempt: 16 CPU cores, 64 GB, 24-hour wall limit,
  at most 20 workflow jobs and 384 allocated CPU-hours. These are pilot caps,
  not runtime predictions; set each bundle's `max_cpu_hours`, `timeout_seconds`
  and Slurm resources accordingly. Initially run two attempts concurrently.
- Review actual timings and convergence before expanding. Do not automatically
  spend the full eight-attempt cap if preparation or state diagnostics fail.

### Stage 2: contrasting enzyme and complete candidates

- Repeat C21/C17 on the accepted CYP17A1 WT model.
- If both enzymes have validated paths, expand each to all five candidate
  carbons, both spin surfaces, and distinct H/pose alternatives.
- Five sites × two spins × two enzymes × two search strategies is **40 nominal
  site/spin/strategy combinations**, before H, pose or replicate expansion.
  Approve an updated measured budget after Stage 1 rather than multiplying
  guessed TS runtimes.

Use Fluxion scores only for candidate prioritization and initial search
hypotheses. Because the archived CYP17A1 run often misses C17, include it as a
fixed diagnostic competitor even when absent from top-k; report **candidate
recall separately** from QM reranking of this augmented candidate set.

## 6. Validation and ranking

For every reported path require:

- SCF and geometry convergence; intended reactant state and intact cofactor.
- Endpoint minima, exactly one significant imaginary mode in the declared
  movable subspace, and an imaginary displacement involving the intended
  C–H/O–H exchange.
- Validated downhill connections to the correct reactant/HAT radical pair.
  For frozen-boundary clusters, use opposite imaginary-mode displacements and
  constrained downhill minimizations, explicitly labelled **movable-subspace
  validation**, not full IRC. Full IRC is reserved for compatible unconstrained
  models.
- Independently audited endpoint graphs with geometry hashes, stereo and metal
  state evidence. The observed benchmark label is not endpoint-validation proof.
- Consistent charge/multiplicity, spin diagnostics, method, basis, environment,
  geometry ordering and energy reference for all competing sites.

Report ΔE‡ and site differences first; add ZPE or thermal corrections only with
consistent, explicitly stated treatment. Partial-Hessian cluster estimates are
not enzyme activation free energies. Do not convert draw fractions or cluster
barriers directly into claimed experimental product yields.

Biological success: CYP21A2 favors C21 over C17, and the accepted CYP17A1 WT
model favors C17 over C21, with the ordering surviving relevant pose/spin/model
sensitivity. Near-degenerate or sensitivity-dependent ordering is inconclusive.
If HAT does not separate sites, investigate pose sampling/rebound rather than
forcing the expected ranking.

Algorithmic success: guided search finds the **same validated chemistry** at
lower attributable compute cost than baseline. Use matched settings/seeds,
include inference, failed attempts and cached source costs, and report CPU-hours,
gradient/Hessian counts, queue and wall time separately. A baseline timeout is
a censored failure, not a measured speedup ratio. Existing manually supplied
synthetic TS seeds are not evidence of Fluxion acceleration on this enzyme.

## 7. Run artifacts and execution sequence

Store under `regioselectivity_evaluation_results/`:

```text
preflight/preflight.json
protocol.json                    # fixed cases, checkpoint, split, budget
models/P08686/chain_A/            # prepared cluster, audit, maps, state variants
models/P05093/wt_model/           # repair/pose provenance plus mutant control
hypotheses/<enzyme>/<site>/<H>/   # elementary HAT JSON + mapped endpoints
runs/<enzyme>/<spin>/<site>/<strategy>/
endpoint_audits/                  # geometry-hashed independent chemistry audits
site_barriers.tsv
search_costs.tsv
evaluation_summary.json
```

Order of work:

1. Locate the archived run's exact checkpoint and training/conditioning store;
   freeze the protocol and regenerate its mapped HAT hypotheses.
2. Implement/test the heme adapter; prepare CYP21A2 Compound I models and review
   electronic states.
3. Generate C21/C17 endpoint geometries and run Stage 1 baseline/guided attempts.
4. Validate paths, inspect failure costs and set the expansion budget.
5. Resolve CYP17A1 WT geometry, run the contrast, then broaden candidates/poses.
6. Publish site ranking, validation evidence, model sensitivity and matched cost
   comparison, including negative and inconclusive results.

Planning estimate: preflight is complete; allow several days for the adapter and
heme/WT model audit, then roughly 1–3 weeks for initial validated paths and
sensitivity checks, depending on queue time and convergence. Measure the first
jobs before promising a full paired-enzyme completion date.
