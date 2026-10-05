"""Create an honest synthetic chemical integration pilot (H + H2 exchange).

This known three-atom step tests mapping, spin, TS curvature, and IRC cheaply.
It is not an enzyme mechanism and is not a Fluxion checkpoint prediction.
"""

import argparse
from pathlib import Path

from qp.reaction.mapping import write_xyz
from qp.reaction.schema import write_json
from qp.reaction.workflow import prepare


def create(output, frozen_environment=False):
    output = Path(output).resolve()
    source = output / "inputs"
    source.mkdir(parents=True, exist_ok=True)
    maps = (11, 23, 47)
    record = {"schema_version": "1.0", "reaction_id": "synthetic_h_h2_exchange",
              "hypothesis_id": "known_step", "hydrogen_policy": "explicit",
              "interpretation": {"kind": "elementary", "step_id": "exchange"},
              "provenance": {"source": "Manually specified synthetic integration reaction; not Fluxion inference",
                             "training_overlap": "not_applicable", "reference_ts_used": True},
              "participants": [{"id": "h2", "role": "substrate"}, {"id": "h", "role": "reactant"}],
              "atoms": [{"map_id": i, "element": "H", "participant": "h2" if i != 47 else "h"} for i in maps],
              "endpoints": {side: {"atom_states": [{"map_id": i, "formal_charge": 0, "radical_electrons": int(i == radical)} for i in maps],
                                    "bonds": [{"atoms": pair, "order": 1.0}]}
                            for side, pair, radical in (("reactant", [11, 23], 47), ("product", [23, 47], 11))},
              "electronic_state": {"charge": 0, "multiplicity": 2}}
    model = {"elements": ["H"] * 3 + (["He"] if frozen_environment else []), "charge": 0, "multiplicity": 2,
             "environment_charge": 0, "reactant_xyz": "reactant.xyz", "product_xyz": "product.xyz",
             "ts_xyz": "ts.xyz", "frozen_indices": [3] if frozen_environment else [],
             "atom_mapping": [{"map_id": i, "qm_index": j, "structure_key": {"synthetic_atom_map": i}}
                              for j, i in enumerate(maps)]}
    # A small asymmetry permits stable minimization in an association complex.
    for name, xyz in (("reactant", [[0., 0., 0.], [.75, 0., 0.], [3.75, .05, 0.]]),
                      ("product", [[0., .05, 0.], [3., 0., 0.], [3.75, 0., 0.]]),
                      ("ts", [[-.93, 0., 0.], [0., 0., 0.], [.93, 0., 0.]])):
        if frozen_environment:
            xyz = xyz + [[10., 3., 2.]]
        write_xyz(source / (name + ".xyz"), model["elements"], xyz)
    settings = {"method": "B3LYP", "basis": "def2-SVP", "keywords": ["TightSCF"],
                "nprocs": 1, "maxcore_mb": 1000, "max_opt_iterations": 100,
                "max_neb_iterations": 100, "neb_images": 6, "max_irc_iterations": 80,
                "max_jobs": 20, "max_cpu_hours": 2, "job_timeout_seconds": 300,
                "walltime": "00:05:00", "cache_dir": str(output / "cache"),
                "poll_seconds": 2, "max_queue_seconds": 600}
    write_json(source / "hypothesis.json", record)
    write_json(source / "model.json", model)
    write_json(source / "settings.json", settings)
    write_json(source / "be_snapshot.json", {
        "elements": ["H"] * 3, "map_ids": list(maps),
        "reactant_be": [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
        "candidates": [{"product_be": [[1, 0, 0], [0, 0, 1], [0, 1, 0]], "count": 1}],
        "provenance": record["provenance"]})
    write_json(source / "export_metadata.json", {k: record[k] for k in
               ("reaction_id", "interpretation", "electronic_state", "provenance")})
    for strategy in ("baseline", "guided", "ts_guess"):
        prepare(source / "hypothesis.json", source / "model.json", output / strategy,
                source / "settings.json", strategy)
    write_json(output / "pilot.json", {"purpose": "Synthetic known-step integration, not measured Fluxion speedup",
                                      "reference_ts_used_in_ts_guess": True,
                                      "strategies": ["baseline", "guided", "ts_guess"],
                                       "baseline_and_graph_guided_route": "Baseline NEB-TS; guided coupled changed-bond scans then OptTS"})
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="reaction_pilot_results")
    parser.add_argument("--frozen-environment", action="store_true", help="Add a frozen distant He atom to verify boundary-constrained TS/frequency/downhill jobs")
    args = parser.parse_args()
    create(args.output, args.frozen_environment)


if __name__ == "__main__":
    main()
