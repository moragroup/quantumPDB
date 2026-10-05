Fluxion-guided reaction searches
===============================

The ``qp.reaction`` package connects atom-mapped reaction hypotheses to
QuantumPDB cluster structures and ORCA 6 calculations. It provides graph-guided
scans, endpoint NEB-TS searches, direct Cartesian TS refinement, frequency and
connectivity validation, resumable execution, and cost reports.

Run its command-line interface from the QuantumPDB environment::

   PY=/mnt/storage01/home/amora/.conda/envs/qp/bin/python
   "$PY" -m qp.reaction --help

NumPy and Biopython are required. RDKit is additionally required when generating
a product geometry from a graph; supplying both endpoint XYZ files also works
in the existing ``qp`` environment. The optional ``reaction`` installation extra
installs these dependencies. ORCA is a separate installation. Here the defaults
are module ``orca/6.1.1`` and executable
``/mnt/nfs/vol8t/software/software/orca/6.1.1/orca``; override these in settings on
another system. Fluxion's inference stack is needed only for checkpoint sampling.

Export a hypothesis
-------------------

The sibling Fluxion repository contains ``scripts/export_qm_hypotheses.py``.
It exports matrices before canonical SMILES conversion removes atom maps.
For an existing prediction snapshot::

   "$PY" ../Fluxion/scripts/export_qm_hypotheses.py \
     --snapshot predictions.json --metadata metadata.json --output exported

``predictions.json`` contains ``elements``, optional persistent ``map_ids``,
``reactant_be``, and ``candidates`` (each with ``product_be`` and integer
``count``). Matrices must already be quantized, symmetric, and electron
conserving. Export records the conversion policy and rejected candidate reasons.
Sampling frequencies are empirical counts, not calibrated probabilities.

Minimum metadata::

   {
     "reaction_id": "enzyme_step_1",
     "interpretation": {"kind": "elementary", "step_id": "1"},
     "electronic_state": {"charge": 0, "multiplicity": 1},
     "provenance": {"training_overlap": "unknown"}
   }

Set the actual charge and multiplicity for your chemistry. Optional
``participants`` entries contain ``id``, ``role``, and ``atom_maps`` and must
partition the atom inventory. Otherwise connected reactant components receive
unassigned participant IDs. Identify substrate, cofactor, reacting residues,
waters, and other reaction partners before preparing the structure mapping.

For checkpoint inference, use a Fluxion environment::

   python ../Fluxion/scripts/export_qm_hypotheses.py \
     --checkpoint /path/to/stepwise_checkpoint.pt \
     --config /path/to/checkpoint_config.yaml \
     --smiles "$REACTANT_SMILES" --sequence "$ENZYME_SEQUENCE" \
     --positions 79,346 --sample-size 32 --seed 2026 --device cuda \
     --metadata metadata.json --output exported

The sequence, catalytic positions, checkpoint, PLM and charge/spin must describe
the selected enzyme and step. Sampling writes ``raw_prediction.json`` and hashes
the checkpoint, configuration, sequence, raw output and exporter. Existing
``build_reactant`` assigns explicit-H maps ``1..N``; the snapshot preserves that
order. BE predictions do not specify Cartesian geometry or stereochemistry.
Stereo-bearing input SMILES require a supplied audited product geometry and
endpoint graph audits.
Start–end/net proposals can be exported with ``kind: net``, but reaction job
preparation requires an explicitly decomposed elementary step.

JSON contract
-------------

Schema version ``1.0`` requires:

* ``reaction_id``, ``hypothesis_id``, ``interpretation`` and ``provenance``;
* ``hydrogen_policy: explicit`` and stable positive chemical atom maps;
* ``participants`` and ``atoms`` (``map_id``, ``element``, ``participant``);
* reactant/product ``endpoints`` containing mapped ``atom_states`` with formal
  charges and radicals, and ``bonds`` with atom pairs and bond orders;
* an integer electronic charge and multiplicity.

The exporter additionally writes bond edits, transferring-H donors/acceptors,
ranking, structure references and BE nonbonding-electron accounting. Validation
checks atom conservation, endpoint charges, electron accounting where supplied,
spin parity, and agreement between declared edits and graphs. A label such as
``explicit`` is a contract: the producer must include every hydrogen and partner
needed by the proposed chemistry.

Map a prepared QuantumPDB cluster
--------------------------------

Prepare and protonate the active-site cluster using QuantumPDB first. Include
every participant involved in the elementary step. Keep caps away from edited
bonds and include catalytic residues that supply coordination, proton transfer
or the essential local environment.

``map-cluster`` matches audited participant graphs by element, topology,
formal charge, radicals and supplied atom/bond stereochemistry. Each prepared
graph atom identifies its actual zero-based position in the concatenated sphere
atom list. Give sphere files in the same order used to write the QM XYZ::

   "$PY" -m qp.reaction map-cluster \
     --hypothesis exported/hypothesis_001.json \
     --spheres cluster/0.pdb cluster/1.pdb cluster/2.pdb \
     --specification mapping_recipe.json --output mapped_model

A minimal mapping recipe for the synthetic H + H2 exchange fixture is::

   {
     "charge": 0, "multiplicity": 2, "environment_charge": 0,
     "frozen_indices": [],
     "graphs": [
       {"participant": "component_1",
        "atoms": [
          {"id": "HA", "qm_index": 0, "element": "H", "formal_charge": 0, "radical_electrons": 0},
          {"id": "HB", "qm_index": 1, "element": "H", "formal_charge": 0, "radical_electrons": 0}],
        "bonds": [{"atoms": ["HA", "HB"], "order": 1}],
        "anchors": {"11": "HA"}},
       {"participant": "component_2",
        "atoms": [{"id": "HC", "qm_index": 2, "element": "H", "formal_charge": 0, "radical_electrons": 1}],
        "bonds": []}
     ],
     "product_xyz": "product.xyz",
     "ts_xyz": "optional_ts.xyz"
   }

Replace these graphs and indices with the prepared chemistry; they are not
inferred from atom names. CCD/SDF topology and an audited prepared protonation
state can supply the graphs. Symmetric mappings require explicit chemical-map
anchors. ``mapping_alternatives.json`` records alternatives; anchors prune the
graph search before symmetry enumeration. Mapping persists sphere, chain,
residue, insertion code, atom name and altloc alongside zero/one-based QM indices.

All chemical atoms must map exactly once. Reactive atoms cannot be frozen or
mapped onto ACE/NME caps. Cluster charge equals chemical plus environment charge;
cluster spin is separately checked. The result contains ``model.json``,
``reactant.xyz`` and ``atom_mapping.tsv``. Supplied product/TS paths are relative
to the recipe location. Covalent links to the unmapped environment require a
supplied product geometry and ``covalent_environment_links`` metadata.

An existing audited map can also be supplied directly as ``model.json``:
``elements``, ``charge``, ``multiplicity``, ``environment_charge``,
``frozen_indices``, ``reactant_xyz``, optional ``product_xyz``/``ts_xyz``, and
``atom_mapping`` rows with ``map_id``, ``qm_index`` and ``structure_key``.
XYZ atom order is persistent throughout the workflow. Element checks cannot
detect a swap between two atoms of the same element; prepared files must retain
the mapped order.

Prepare and execute searches
---------------------------

Example ``settings.json``::

   {
     "method": "B3LYP", "basis": "def2-SVP",
     "keywords": ["D3BJ", "TightSCF"],
     "nprocs": 8, "maxcore_mb": 2000,
     "partition": "cpu", "walltime": "02:00:00",
     "max_jobs": 20, "max_cpu_hours": 80,
     "job_timeout_seconds": 7200,
     "max_opt_iterations": 150, "max_neb_iterations": 150,
     "neb_images": 8, "max_irc_iterations": 60,
     "scan_steps": 9, "imaginary_cutoff": -50,
     "minimum_mode_overlap": 0.2
   }

Optional ``dielectric`` enables CPCM; ``pointcharges`` supplies a copied external
charge file. Optional ``orca_executable``, ``orca_module``, ``cache_dir``,
``poll_seconds`` and ``max_queue_seconds`` control execution. Resources and
chemical methods should be selected for the model being studied.

Create one bundle per candidate/search strategy::

   "$PY" -m qp.reaction prepare \
     --hypothesis exported/hypothesis_001.json --model mapped_model/model.json \
     --settings settings.json --strategy guided --output searches/candidate_001
   "$PY" -m qp.reaction run searches/candidate_001 --scheduler slurm
   "$PY" -m qp.reaction collect searches/candidate_001

Strategies:

``baseline``
   Optimize the two endpoints, then NEB-TS and TS validation.
``guided``
   Optimize endpoints, derive coupled distance scans for up to three changed
   bonds, select the highest-energy converged scan point and refine with OptTS.
   Larger edit sets, unchanged distances or ``auto_scan: false`` use NEB-TS.
   Explicit ``scans`` settings can override distances and step counts.
``ts_guess``
   Direct OptTS from an actual mapped Cartesian TS proposal, plus endpoint
   optimization. The TS seed does not undergo ordinary minimum optimization.
``screen``
   Optimize and check both endpoints at a separately chosen low-cost method.
   Result ``screened_endpoints`` establishes endpoint screening only.

Heavy atoms are movable except for explicitly listed boundary atoms. Unfrozen
models use full frequencies and a two-way mass-weighted IRC, followed by
optimization and frequency checks at both connected endpoints. Frozen models
use numerical partial Hessians over movable atoms, with an explicit check of
ORCA's reported displaced-atom subset. They validate connectivity through
opposite imaginary-mode displacements and boundary-constrained downhill
optimizations. This is reported as movable-subspace validation, not an IRC.

``run --scheduler local`` executes locally; ``--submit-only`` submits the next
ready Slurm job. Repeating ``run`` resumes dependencies and checks existing
submissions. A running submitted job is left for subsequent polling/resumption.
Bundles reject changed input geometries, settings, staged input/scripts, or
engine identity. Prepare a new output directory for a changed search. Execution
stops on a failed job or exhaustion of job-count/reserved CPU-hour budgets;
timeouts terminate the job process group. Failed attempts remain in the reports.

Exact-input caches verify output hashes, executable identity and source cost.
The TS frequency Hessian is reused by IRC only at its matching geometry and
method. Cache reuse contributes original attributable cost to comparisons.

Interpret results and compare costs
----------------------------------

Each bundle contains ``hypothesis.json``, ``model.json``, ``atom_mapping.tsv``,
prepared XYZ files, ``manifest.json``, ``state.json``, candidate job inputs/logs,
``results.json`` and ``cost_report.tsv``. Manifests include input, geometry,
implementation and executable hashes plus dependency versions and budgets.

``validated_step`` requires converged endpoints and TS, endpoint minima,
exactly one significant reaction-aligned imaginary mode, and connectivity to
both proposed endpoint graphs. Normal ORCA termination alone is insufficient.
An unrefined NEB climbing image is not accepted as the refined TS. Other result
statuses include ``not_converged`` and ``wrong_connectivity``. Invalid mapping
or unresolved charge/spin inputs fail preparation with a diagnostic and exit 2.

Automatic distance connectivity is limited to simple nonmetal single-bond,
non-stereochemical hypotheses. Bond-order, metal and stereochemical chemistry
requires ``endpoint_graphs.json`` entries for ``reactant``, ``product``,
``connected_a`` and ``connected_b``. Each entry must include the final XYZ
``geometry_sha256``, graph-assignment ``method``, ``provenance`` and ``graph``
in endpoint format. This is explicit audited evidence, not a geometry-based
determination of formal charge or spin.

Barriers are electronic cluster energies in kcal/mol at consistent settings.
They are not enzyme free energies; partial-Hessian thermochemistry is not
reported as a full free-energy correction. Spin/protonation, cluster size,
environment and alternative elementary steps still require physical assessment.

Compare matched strategies and include separately incurred screening costs::

   "$PY" -m qp.reaction compare searches/baseline searches/guided \
     --additional-cost-bundle searches/cheap_screen --output comparison.json

Comparison checks identical chemistry, mapped environment/geometries, engine
and settings. Reports include CPU/GPU allocation, queue and run time, failures,
and observed energy/gradient/Hessian counters (output-marker lower bounds).
Speedup is emitted only for matched validated runs with complete cost accounting
and no reference TS. A failed search does not establish an impossible mechanism.

Verified integration pilot
--------------------------

Create a fresh, small H + H2 doublet exchange bundle::

   "$PY" -m qp.reaction.pilot --output my_reaction_pilot
   "$PY" ../Fluxion/scripts/export_qm_hypotheses.py \
     --snapshot my_reaction_pilot/inputs/be_snapshot.json \
     --metadata my_reaction_pilot/inputs/export_metadata.json \
     --output my_reaction_pilot/exported
   "$PY" -m qp.reaction run my_reaction_pilot/ts_guess --scheduler slurm

The pilot is manually specified chemistry with a known TS seed. It tests the
bridge and engine validation rather than enzyme mechanism inference. Actual
ORCA 6.1.1 B3LYP/def2-SVP runs are retained under
``reaction_pilot_results/validated/ts_guess``:

* 11 converged jobs, one imaginary mode at -508.11 cm^-1;
* reaction-mode overlap 0.99999984;
* both IRC branches converged and connected to the mapped endpoint minima;
* forward electronic barrier 3.08882 kcal/mol.

``reaction_pilot_results/boundary_final/ts_guess`` additionally verified frozen
boundaries, the correct movable-atom partial Hessian and opposite-mode downhill
connectivity. It has ``validation_scope: movable_subspace``. The endpoint-only
NEB pilot reached a converged band but failed TS refinement; the coupled-scan
pilot failed an internal-coordinate transformation. Both retain their costs
and ``not_converged`` status. ``reaction_pilot_results/validated/comparison.json``
therefore reports no speedup. Checkpoint-predicted enzyme pathways still need an
audited enzyme-specific elementary step, complete participant mapping and actual
matched QM searches.
