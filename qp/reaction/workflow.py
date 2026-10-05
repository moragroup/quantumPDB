"""Reproducible staged searches, exact-input caches, and bounded execution."""

import argparse
import json
import os
import signal
import shutil
import subprocess
import sys
import time
from pathlib import Path

from qp.reaction.geometry import prepare_geometries
from qp.reaction.jobs import write_job
from qp.reaction.mapping import build_cluster_model, read_xyz, validate_mapping, write_mapping, write_xyz
from qp.reaction.results import (endpoint_matches, mode_overlap, normal_modes,
                                 parse_orca, validate_step, write_cost_report)
from qp.reaction.schema import (HypothesisError, content_hash, load, sha256,
                                reaction_edits, write_json)


def _engine_identity(settings):
    executable = Path(settings.get("orca_executable", "/mnt/nfs/vol8t/software/software/orca/6.1.1/orca"))
    return {"module": settings.get("orca_module", "orca/6.1.1"),
             "executable": str(executable), "executable_sha256": sha256(executable) if executable.is_file() else None}


def _versions():
    from importlib.metadata import PackageNotFoundError, version
    result = {"python": sys.version}
    for package in ("numpy", "biopython", "rdkit"):
        try:
            result[package] = version(package)
        except PackageNotFoundError:
            result[package] = None
    return result


def prepare(hypothesis, model_path, output, settings_path=None, strategy="guided"):
    record = load(hypothesis)
    model_path, output = Path(model_path).resolve(), Path(output).resolve()
    if (output / "state.json").exists():
        raise HypothesisError("Bundle already exists; resume it with run or choose a new output directory")
    model = validate_mapping(record, json.loads(model_path.read_text()))
    settings = json.loads(Path(settings_path).read_text()) if settings_path else {}
    budgets = {"max_jobs": settings.get("max_jobs", 20), "max_cpu_hours": settings.get("max_cpu_hours", 10),
               "job_timeout_seconds": settings.get("job_timeout_seconds", 3600)}
    if type(budgets["max_jobs"]) is not int or any(not isinstance(v, (int, float)) or v <= 0 for v in budgets.values()):
        raise HypothesisError("Positive job-count, CPU-hour and timeout budgets required")
    if strategy not in ("baseline", "guided", "ts_guess", "screen"):
        raise HypothesisError("Unknown strategy")
    if strategy == "ts_guess" and not model.get("ts_xyz"):
        raise HypothesisError("ts_guess strategy requires an actual Cartesian TS proposal")
    if settings.get("screening"):
        # Screening can be run as a separate, fully costed bundle, not as an
        # unrecorded method change inside a supposedly matched DFT comparison.
        raise HypothesisError("Prepare a separate low-cost screening bundle and record its costs in the comparison")
    model_dir = output / "model"
    model_dir.mkdir(parents=True, exist_ok=True)
    prepare_geometries(record, model, model_path.parent, model_dir)
    write_json(output / "hypothesis.json", record)
    write_json(output / "model.json", model)
    write_json(output / "settings.json", settings)
    write_mapping(output / "atom_mapping.tsv", record, model)
    if settings.get("pointcharges"):
        source = Path(settings["pointcharges"])
        if not source.is_absolute():
            source = model_path.parent / source
        shutil.copyfile(source, model_dir / "pointcharges.pc")
        settings["pointcharges"] = "pointcharges.pc"
        write_json(output / "settings.json", settings)
    manifest = {"schema_version": "1.0", "strategy": strategy,
                "hypothesis_hash": content_hash(record), "model_hash": content_hash(model),
                "settings_hash": content_hash(settings), "input_hashes": {
                    str(Path(hypothesis).resolve()): sha256(hypothesis), str(model_path): sha256(model_path)},
                "prepared_files": {p.name: sha256(p) for p in model_dir.iterdir() if p.is_file()},
                 "implementation_hashes": {p.name: sha256(p) for p in Path(__file__).parent.glob("*.py")},
                  "engine_identity": _engine_identity(settings),
                  "versions": _versions(),
                 "budgets": budgets,
                "created_unix": time.time(), "training_overlap": record["provenance"].get("training_overlap", "unknown")}
    write_json(output / "manifest.json", manifest)
    write_json(output / "state.json", {"jobs": [], "status": "prepared"})
    stage_jobs(output)
    return output


def _context(root):
    root = Path(root).resolve()
    record = load(root / "hypothesis.json")
    model = validate_mapping(record, json.loads((root / "model.json").read_text()))
    settings = json.loads((root / "settings.json").read_text())
    state = json.loads((root / "state.json").read_text())
    manifest = json.loads((root / "manifest.json").read_text())
    for key, value in (("hypothesis_hash", record), ("model_hash", model), ("settings_hash", settings)):
        if content_hash(value) != manifest[key]:
            raise HypothesisError("Bundle inputs changed; prepare a new bundle rather than reuse incompatible jobs")
    for name, digest in manifest["prepared_files"].items():
        if sha256(root / "model" / name) != digest:
            raise HypothesisError("Prepared geometry/environment changed")
    if manifest.get("engine_identity") and manifest["engine_identity"] != _engine_identity(settings):
        raise HypothesisError("QM executable identity changed; prepare a new bundle")
    for job in state["jobs"]:
        for name, digest in job["input_files"].items():
            if sha256(root / job["directory"] / name) != digest:
                raise HypothesisError("Staged job inputs changed: {} / {}".format(job["job"], name))
    return root, record, model, settings, state, manifest


def _job(root, name, mode, model, settings, xyz, product=None, hessian=None, scans=None):
    directory = root / "candidates" / "001" / name
    directory.mkdir(parents=True, exist_ok=True)
    elements, coords = read_xyz(xyz, last=True)
    if elements != model["elements"]:
        raise HypothesisError("Engine geometry does not preserve QM atom inventory")
    write_xyz(directory / "input.xyz", elements, coords)
    kwargs = {"xyz": "input.xyz", "scans": scans}
    if product:
        other, coords = read_xyz(product, last=True)
        if other != elements:
            raise HypothesisError("Engine product ordering differs")
        write_xyz(directory / "product.xyz", elements, coords)
        kwargs["product"] = "product.xyz"
    if hessian:
        shutil.copyfile(hessian, directory / "initial.hess")
        kwargs["hessian"] = "initial.hess"
    if settings.get("pointcharges"):
        shutil.copyfile(root / "model" / "pointcharges.pc", directory / "pointcharges.pc")
    write_job(directory, mode, model, settings, **kwargs)
    inputs = {p.name: sha256(p) for p in directory.iterdir() if p.name in
              ("job.inp", "input.xyz", "product.xyz", "initial.hess", "pointcharges.pc", "run.sh", "submit.sh")}
    fingerprint = content_hash({"inputs": inputs, "model": model, "settings": settings,
                                "engine": _engine_identity(settings)})
    return {"job": name, "mode": mode, "directory": str(directory.relative_to(root)), "input_files": inputs,
            "input_hash": fingerprint, "status": "ready", "cached": False}


def _geometry(root, job):
    directory = root / job["directory"]
    # ORCA's final optimizer geometry, including the TS refinement after NEB.
    names = ("job_NEB-TS_converged.xyz", "job.xyz") if job["mode"] == "path_search" else ("job.xyz",)
    for name in names:
        if (directory / name).exists():
            return directory / name
    raise HypothesisError("Converged job has no final XYZ: {}".format(job["job"]))


def stage_jobs(bundle):
    root, record, model, settings, state, manifest = _context(bundle)
    jobs = {j["job"]: j for j in state["jobs"]}

    def add(name, mode, xyz, **kwargs):
        if name not in jobs:
            jobs[name] = _job(root, name, mode, model, settings, xyz, **kwargs)

    def done(name):
        return name in jobs and jobs[name]["status"] == "converged"

    for side in ("reactant", "product"):
        add(side, "endpoint_opt", root / "model" / (side + ".xyz"))
        if done(side) and manifest["strategy"] != "screen":
            add(side + "_frequency", "frequency", _geometry(root, jobs[side]))
    if manifest["strategy"] == "screen":
        state["jobs"] = list(jobs.values())
        write_json(root / "state.json", state)
        return state
    if manifest["strategy"] == "ts_guess":
        add("ts", "ts_opt", root / "model" / "ts_guess.xyz")
    elif done("reactant") and done("product"):
        # Optional explicitly configured edit-guided scan, absent in baseline.
        scans = settings.get("scans")
        if manifest["strategy"] == "guided" and not scans and settings.get("auto_scan", True):
            scans = _edit_scans(record, model, _geometry(root, jobs["reactant"]), _geometry(root, jobs["product"]), settings)
        if manifest["strategy"] == "guided" and scans:
            add("scan", "constrained_scan", _geometry(root, jobs["reactant"]), scans=scans)
            if done("scan"):
                directory = root / jobs["scan"]["directory"]
                trajectory = directory / "job.allxyz"
                frames = _xyz_frames(trajectory)
                energies = _scan_energies(directory / "job.out")
                if not energies:
                    # ORCA scan trajectory comments carry point energies. Avoid
                    # pairing them with all SCF iterations in the optimizer log.
                    import re
                    energies = [float(e) for e in re.findall(r"(?im)^.*\bE\s+([-+\d.Ee]+)\s*$", trajectory.read_text())]
                if len(frames) != len(energies):
                    raise HypothesisError("Cannot pair scan geometries with energies for TS seeding")
                seed = frames[max(range(len(energies)), key=lambda i: energies[i])]
                write_xyz(directory / "scan_seed.xyz", *seed)
                add("ts", "ts_opt", directory / "scan_seed.xyz")
        else:
            add("ts", "path_search", _geometry(root, jobs["reactant"]), product=_geometry(root, jobs["product"]))
    if done("ts"):
        add("ts_frequency", "frequency", _geometry(root, jobs["ts"]))
    if done("ts_frequency"):
        if model.get("frozen_indices"):
            # ORCA IRC does not retain Cartesian boundary constraints. Instead
            # displace both ways along the calculated imaginary mode and optimize
            # with the same boundaries, then require two endpoint minima/graphs.
            directory = root / jobs["ts_frequency"]["directory"]
            parsed = parse_orca(directory / "job.out", "frequency")
            imaginary = [f for f in parsed["frequencies_cm_1"] if f["frequency"] < settings.get("imaginary_cutoff", -50)]
            if len(imaginary) == 1:
                elements, xyz = read_xyz(_geometry(root, jobs["ts"]))
                vector = normal_modes(directory / "job.hess")[:, imaginary[0]["index"]].reshape(-1, 3)
                vector[model["frozen_indices"]] = 0
                import numpy as np
                size = np.linalg.norm(vector, axis=1).max()
                if size < 1e-8:
                    raise HypothesisError("Imaginary mode has no movable displacement")
                for side, sign in (("connected_a", 1), ("connected_b", -1)):
                    seed = directory / (side + "_seed.xyz")
                    write_xyz(seed, elements, xyz + sign * settings.get("downhill_displacement_A", .15) * vector / size)
                    add(side, "endpoint_opt", seed)
                    if done(side):
                        add(side + "_frequency", "frequency", _geometry(root, jobs[side]))
            else:
                state["status"] = "invalid_ts_curvature"
        else:
            hessian = root / jobs["ts_frequency"]["directory"] / "job.hess"
            add("irc", "irc", _geometry(root, jobs["ts"]), hessian=hessian if hessian.exists() else None)
    if done("irc"):
        directory = root / jobs["irc"]["directory"]
        for side, suffix in (("connected_a", "F"), ("connected_b", "B")):
            endpoint = directory / ("job_IRC_{}.xyz".format(suffix))
            if not endpoint.exists():
                endpoint = directory / ("job_IRC_{}_trj.xyz".format(suffix))
            add(side, "endpoint_opt", endpoint)
            if done(side):
                add(side + "_frequency", "frequency", _geometry(root, jobs[side]))
    state["jobs"] = list(jobs.values())
    write_json(root / "state.json", state)
    return state


def _edit_scans(record, model, reactant, product, settings):
    """Turn at most three changed bonds into coupled endpoint-distance scans."""
    import numpy as np
    edits = reaction_edits(record)
    if not edits or len(edits) > 3:
        return None  # Use an endpoint path when the scan dimension is excessive.
    _, left = read_xyz(reactant)
    _, right = read_xyz(product)
    index = {r["map_id"]: r["qm_index"] for r in model["atom_mapping"]}
    scans = []
    for edit in edits:
        a, b = [index[i] for i in edit["atoms"]]
        start, end = [float(np.linalg.norm(x[b] - x[a])) for x in (left, right)]
        if abs(end - start) < .05:
            continue
        scans.append({"indices": [a, b], "start": round(start, 5), "end": round(end, 5),
                      "steps": settings.get("scan_steps", 9), "chemical_maps": edit["atoms"]})
    return scans or None


def _xyz_frames(path):
    lines = Path(path).read_text().splitlines()
    import numpy as np
    frames, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i])
        block = [line.split() for line in lines[i + 2:i + n + 2]]
        frames.append(([r[0] for r in block], np.asarray([[float(v) for v in r[1:4]] for r in block])))
        i += n + 2
    return frames


def _scan_energies(path):
    import re
    text = Path(path).read_text()
    # ORCA's converged point-energy table, not intermediate SCF energies.
    match = re.search(r"RELAXED SURFACE SCAN RESULTS(.*?)(?:\n\s*\n|\Z)", text, re.S | re.I)
    if not match:
        return []
    return [float(line.split()[-1]) for line in match.group(1).splitlines()
            if re.fullmatch(r"\s*[-+\d.Ee]+(?:\s+[-+\d.Ee]+)+\s*", line)]


def execute(bundle, scheduler="local", submit_only=False):
    """Advance dependencies; repeated invocation resumes without duplicating jobs."""
    root, _, _, settings, state, manifest = _context(bundle)
    cache = Path(settings.get("cache_dir", root / "cache")).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    while True:
        state = stage_jobs(root)
        _collect_jobs(root, state, settings)
        write_json(root / "state.json", state)
        if any(j["status"] in ("failed", "timeout", "not_converged") for j in state["jobs"]):
            break
        if any(j["status"] == "submitted" for j in state["jobs"]):
            break
        ready = [j for j in state["jobs"] if j["status"] == "ready"]
        if not ready:
            # Collection may have unlocked a new dependency; stage once more.
            before = len(state["jobs"])
            state = stage_jobs(root)
            if len(state["jobs"]) > before:
                continue
            break
        job = ready[0]
        directory = root / job["directory"]
        cached = cache / job["input_hash"]
        if (cached / "evidence.json").exists():
            cached_evidence = json.loads((cached / "evidence.json").read_text())
            digests = cached_evidence.get("output_files", {})
            if (cached_evidence.get("converged") and _engine_identity(settings)["executable_sha256"] and digests and
                    "allocated_cpu_hours" in cached_evidence.get("source_cost", {}) and
                    all((cached / n).exists() and sha256(cached / n) == h for n, h in digests.items())):
                for p in cached.iterdir():
                    if p.name != "evidence.json":
                        shutil.copyfile(p, directory / p.name)
                job.update(status="converged", cached=True, wall_seconds=0,
                            allocated_cpu_hours=0, allocated_gpu_hours=0,
                            cache_source_cost=cached_evidence.get("source_cost", {}),
                           cost={"energy_evaluations_observed": 0, "gradient_evaluations_observed": 0,
                                 "hessian_evaluations_observed": 0})
                write_json(root / "state.json", state)
                continue
        charged = [j for j in state["jobs"] if j.get("started_unix") and not j.get("cached")]
        budget = manifest["budgets"]
        used_cpu = sum(j.get("allocated_cpu_hours", 0) for j in charged)
        reserved = settings.get("nprocs", 1) * budget["job_timeout_seconds"] / 3600.
        if len(charged) >= budget["max_jobs"] or used_cpu + reserved > budget["max_cpu_hours"]:
            state["status"] = "budget_exhausted"
            state["reason"] = "Next job would exceed job count or reserved CPU-hour budget"
            write_json(root / "state.json", state)
            break
        job["started_unix"] = time.time()
        if scheduler == "slurm":
            response = subprocess.run(["sbatch", "--parsable", "submit.sh"], cwd=str(directory),
                                      check=True, text=True, capture_output=True)
            job.update(status="submitted", slurm_job_id=response.stdout.strip().split(";")[0],
                       reserved_cpu_hours=reserved)
            write_json(root / "state.json", state)
            if submit_only:
                return state
            while _slurm_active(job["slurm_job_id"]):
                if time.time() - job["started_unix"] > budget["job_timeout_seconds"] + settings.get("max_queue_seconds", 3600):
                    subprocess.run(["scancel", job["slurm_job_id"]], check=True)
                    job["status"] = "timeout"
                    break
                time.sleep(settings.get("poll_seconds", 10))
        else:
            process = subprocess.Popen(["bash", "run.sh"], cwd=str(directory), start_new_session=True)
            try:
                returncode = process.wait(timeout=budget["job_timeout_seconds"])
                job["returncode"] = returncode
                job["status"] = "finished" if returncode == 0 else "failed"
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                job["status"] = "timeout"
        job["finished_unix"] = time.time()
        job["wall_seconds"] = job["finished_unix"] - job["started_unix"]
        _collect_jobs(root, state, settings)
        if job["status"] == "converged":
            cached.mkdir(parents=True, exist_ok=True)
            for p in directory.iterdir():
                if p.suffix in (".out", ".xyz", ".hess", ".allxyz"):
                    shutil.copyfile(p, cached / p.name)
            cache_evidence = parse_orca(directory / "job.out", job["mode"])
            cache_evidence["output_files"] = {p.name: sha256(p) for p in cached.iterdir() if p.name != "evidence.json"}
            cache_evidence["source_cost"] = {key: job.get(key, 0) for key in
                                            ("allocated_cpu_hours", "allocated_gpu_hours", "wall_seconds")}
            cache_evidence["source_cost"]["cost"] = job.get("cost", {})
            write_json(cached / "evidence.json", cache_evidence)
        write_json(root / "state.json", state)
    collect(root)
    return json.loads((root / "state.json").read_text())


def _slurm_active(job_id):
    return bool(subprocess.run(["squeue", "-h", "-j", job_id], check=True,
                               text=True, capture_output=True).stdout.strip())


def _collect_jobs(root, state, settings):
    for job in state["jobs"]:
        if job["status"] == "submitted" and _slurm_active(job["slurm_job_id"]):
            continue
        path = root / job["directory"] / "job.out"
        if job["status"] == "submitted" and not path.exists():
            accounting = subprocess.run(["sacct", "-n", "-P", "-j", job["slurm_job_id"],
                                         "--format=JobID,State,ExitCode,ElapsedRaw"],
                                        check=True, text=True, capture_output=True).stdout
            for line in accounting.splitlines():
                fields = line.split("|")
                if fields[0] == job["slurm_job_id"] and fields[1] not in ("PENDING", "RUNNING", "CONFIGURING", "COMPLETING"):
                    job.update(status="failed", slurm_accounting=accounting,
                               allocated_cpu_hours=float(fields[3]) * settings.get("nprocs", 1) / 3600,
                               wall_seconds=float(fields[3]))
        if not path.exists() or job["status"] == "ready":
            continue
        evidence = parse_orca(path, job["mode"])
        if job["status"] not in ("timeout", "failed"):
            job["status"] = "converged" if evidence["converged"] else "not_converged"
        if not job.get("cached"):
            job["cost"] = evidence["cost"]
            seconds = evidence["engine_wall_seconds"] or job.get("wall_seconds", 0)
            job["allocated_cpu_hours"] = seconds * settings.get("nprocs", 1) / 3600
            job["allocated_gpu_hours"] = 0
            job["wall_seconds"] = seconds
            if job.get("slurm_job_id"):
                accounting = subprocess.run(["sacct", "-n", "-P", "-j", job["slurm_job_id"],
                                             "--format=JobID,ElapsedRaw,Submit,Start,State"],
                                            check=True, text=True, capture_output=True).stdout
                job["slurm_accounting"] = accounting
                # Keep queue timing absent if scheduler records are not yet available.
                for line in accounting.splitlines():
                    fields = line.split("|")
                    if fields[0] == job["slurm_job_id"] and len(fields) >= 5:
                        try:
                            from datetime import datetime
                            job["queue_seconds"] = (datetime.fromisoformat(fields[3]) - datetime.fromisoformat(fields[2])).total_seconds()
                            job["allocated_cpu_hours"] = float(fields[1]) * settings.get("nprocs", 1) / 3600
                            job["wall_seconds"] = float(fields[1])
                        except ValueError:
                            pass


def collect(bundle):
    root, record, model, settings, state, manifest = _context(bundle)
    _collect_jobs(root, state, settings)
    evidence = {}
    for job in state["jobs"]:
        directory = root / job["directory"]
        if not (directory / "job.out").exists():
            continue
        item = parse_orca(directory / "job.out", job["mode"])
        item["converged"] = item["converged"] and job["status"] == "converged"
        if job["mode"] == "frequency" and model.get("frozen_indices") and item["converged"]:
            expected = set(range(len(model["elements"]))) - set(model["frozen_indices"])
            if set(item.get("displaced_atoms_0based") or []) != expected:
                item["converged"] = False
                item["subspace_error"] = "Numerical Hessian displaced a different atom subset from the movable model"
        item["output_sha256"] = sha256(directory / "job.out")
        if item["converged"] and job["mode"] in ("endpoint_opt", "ts_opt", "path_search"):
            geometry = _geometry(root, job)
            elements, xyz = read_xyz(geometry)
            if elements != model["elements"]:
                raise HypothesisError("Final engine atom inventory/order differs")
            if model.get("frozen_indices"):
                import numpy as np
                _, initial = read_xyz(directory / "input.xyz")
                frozen = model["frozen_indices"]
                if not np.allclose(xyz[frozen], initial[frozen], atol=1e-4, rtol=0):
                    raise HypothesisError("Engine moved frozen boundary atoms")
            item["geometry_sha256"] = sha256(geometry)
        if job["job"] == "ts_frequency" and item["converged"]:
            item["hessian_sha256"] = sha256(directory / "job.hess") if (directory / "job.hess").exists() else None
            imaginary = [f for f in item["frequencies_cm_1"] if f["frequency"] < settings.get("imaginary_cutoff", -50)]
            if len(imaginary) == 1 and (directory / "job.hess").exists():
                _, xyz = read_xyz(directory / "input.xyz")
                modes = normal_modes(directory / "job.hess")
                item["reaction_mode_overlap"] = mode_overlap(record, model, xyz, modes[:, imaginary[0]["index"]].reshape(-1, 3))
        if job["job"] in ("reactant", "product", "connected_a", "connected_b") and item["converged"]:
            _, xyz = read_xyz(_geometry(root, job))
            item["matches"] = []
            audit = root / "endpoint_graphs.json"
            graphs = json.loads(audit.read_text()) if audit.exists() else {}
            audited = graphs.get(job["job"])
            if audited is not None:
                geometry_hash = sha256(_geometry(root, job))
                if (audited.get("geometry_sha256") != geometry_hash or not audited.get("provenance") or not audited.get("method")):
                    raise HypothesisError("Endpoint graph audits require the final geometry hash, method and provenance")
                item["graph_audit"] = audited
            for side in ("reactant", "product"):
                try:
                    if endpoint_matches(record, model, xyz, side, audited.get("graph") if audited else None):
                        item["matches"].append(side)
                except HypothesisError as exc:
                    item["audit_required"] = str(exc)
        evidence[job["job"]] = item
    if model.get("frozen_indices"):
        evidence["irc"] = {
            "converged": all(evidence.get(s, {}).get("converged") for s in ("connected_a", "connected_b")),
            "method": "opposite imaginary-mode displacements and boundary-constrained endpoint optimization",
            "scope": "downhill connectivity evidence; not a mass-weighted IRC trajectory",
            "displacement_A": settings.get("downhill_displacement_A", .15)}
    # Connected endpoints must themselves be minima, not just completed optimizations.
    for side in ("connected_a", "connected_b"):
        freq = evidence.get(side + "_frequency")
        if side in evidence:
            if not freq or not freq["converged"] or any(f["frequency"] < settings.get("imaginary_cutoff", -50) for f in freq["frequencies_cm_1"]):
                evidence[side]["converged"] = False
    result = validate_step(record, model, evidence, settings.get("imaginary_cutoff", -50), settings.get("minimum_mode_overlap", .2))
    if manifest["strategy"] == "screen":
        accepted = all(evidence.get(s, {}).get("converged") and s in evidence[s].get("matches", []) for s in ("reactant", "product"))
        result.update(status="screened_endpoints" if accepted else "screening_failed",
                      validation_scope="endpoint convergence/connectivity only; no stationary-point curvature or TS evidence",
                      reasons=[] if accepted else ["Endpoints failed convergence/connectivity screening"])
    result.update(hypothesis_id=record["hypothesis_id"], reaction_id=record["reaction_id"],
                  strategy=manifest["strategy"], model_hash=manifest["model_hash"],
                   settings_hash=manifest["settings_hash"], training_overlap=manifest["training_overlap"],
                   reference_ts_used=record["provenance"].get("reference_ts_used", False),
                    runtime_implementation_hashes={p.name: sha256(p) for p in Path(__file__).parent.glob("*.py")})
    result["problem_hash"] = content_hash({
        "chemistry": {k: record[k] for k in ("atoms", "endpoints", "electronic_state", "interpretation")},
        "prepared_files": manifest["prepared_files"], "model_hash": manifest["model_hash"],
        "engine": manifest.get("engine_identity")})
    if state.get("status") == "budget_exhausted":
        result["reasons"].append(state["reason"])
    state["validation_status"] = result["status"]
    write_json(root / "state.json", state)
    write_json(root / "results.json", result)
    write_cost_report(root / "cost_report.tsv", state["jobs"])
    return result


def compare(bundles, output, additional_cost_bundles=()):
    results, rows = [], []
    for path in list(bundles) + list(additional_cost_bundles):
        path = Path(path)
        result = collect(path)
        jobs = json.loads((path / "state.json").read_text())["jobs"]
        rows.append({"bundle": str(path.resolve()), "strategy": result["strategy"],
                     "status": result["status"], "cpu_hours": sum(j.get("allocated_cpu_hours", 0) for j in jobs),
                      "gpu_hours": sum(j.get("allocated_gpu_hours", 0) for j in jobs),
                      "attributable_cpu_hours": sum(j.get("cache_source_cost", {}).get("allocated_cpu_hours", 0)
                                                    if j.get("cached") else j.get("allocated_cpu_hours", 0) for j in jobs),
                      "attributable_gpu_hours": sum(j.get("cache_source_cost", {}).get("allocated_gpu_hours", 0)
                                                    if j.get("cached") else j.get("allocated_gpu_hours", 0) for j in jobs),
                     "failed_jobs": sum(j["status"] in ("failed", "timeout", "not_converged") for j in jobs),
                     "wall_seconds_sum": sum(j.get("wall_seconds", 0) for j in jobs),
                     "queue_seconds_sum": sum(j.get("queue_seconds", 0) for j in jobs),
                       "additional_screening_cost": path in [Path(p) for p in additional_cost_bundles]})
        rows[-1]["cost_complete"] = all(not j.get("cached") or "allocated_cpu_hours" in j.get("cache_source_cost", {})
                                        for j in jobs)
        for counter in ("energy_evaluations_observed", "gradient_evaluations_observed", "hessian_evaluations_observed"):
            rows[-1][counter] = sum(j.get("cache_source_cost", {}).get("cost", {}).get(counter, 0)
                                   if j.get("cached") else j.get("cost", {}).get(counter, 0) for j in jobs)
        if not rows[-1]["additional_screening_cost"]:
            results.append(result)
    compatible = len({(r["reaction_id"], r["problem_hash"], r["settings_hash"]) for r in results}) == 1
    report = {"matched_settings": compatible, "runs": rows, "speedup": None,
              "counter_semantics": "Observed output markers are lower bounds; allocated CPU hours include failed jobs",
               "interpretation": "No speedup claim without matched settings and validated runs; shared cached jobs are charged at original source cost",
               "reference_ts_used": any(r["reference_ts_used"] for r in results)}
    if (compatible and results and all(r["status"] == "validated_step" for r in results)
            and all(r["cost_complete"] for r in rows) and not report["reference_ts_used"]):
        baseline = [r for r in rows if r["strategy"] == "baseline"]
        guided = [r for r in rows if r["strategy"] != "baseline" and not r["additional_screening_cost"]]
        screening_cpu = sum(r["attributable_cpu_hours"] for r in rows if r["additional_screening_cost"])
        if len(baseline) == 1:
            report["speedup"] = {r["bundle"]: baseline[0]["attributable_cpu_hours"] / (r["attributable_cpu_hours"] + screening_cpu)
                                   for r in guided if r["attributable_cpu_hours"] + screening_cpu > 0}
    write_json(output, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare")
    p.add_argument("--hypothesis", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--settings")
    p.add_argument("--output", required=True)
    p.add_argument("--strategy", choices=("baseline", "guided", "ts_guess", "screen"), default="guided")
    p = commands.add_parser("map-cluster")
    p.add_argument("--hypothesis", required=True)
    p.add_argument("--spheres", nargs="+", required=True, help="Numbered sphere PDBs in exact XYZ concatenation order")
    p.add_argument("--specification", required=True, help="Prepared participant graphs, anchors, electronic state and boundary indices")
    p.add_argument("--output", required=True)
    p = commands.add_parser("run")
    p.add_argument("bundle")
    p.add_argument("--scheduler", choices=("local", "slurm"), default="local")
    p.add_argument("--submit-only", action="store_true")
    p = commands.add_parser("collect")
    p.add_argument("bundle")
    p = commands.add_parser("compare")
    p.add_argument("bundles", nargs="+")
    p.add_argument("--additional-cost-bundle", action="append", default=[])
    p.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "prepare":
            prepare(args.hypothesis, args.model, args.output, args.settings, args.strategy)
        elif args.command == "map-cluster":
            path = Path(args.specification).resolve()
            build_cluster_model(load(args.hypothesis), args.spheres, json.loads(path.read_text()), args.output, path.parent)
        elif args.command == "run":
            execute(args.bundle, args.scheduler, args.submit_only)
        elif args.command == "collect":
            collect(args.bundle)
        else:
            compare(args.bundles, args.output, args.additional_cost_bundle)
    except (HypothesisError, KeyError, ValueError, OSError) as exc:
        parser.exit(2, "Reaction workflow failed: {}\n".format(exc))


if __name__ == "__main__":
    main()
