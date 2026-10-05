"""ORCA 6 reaction modes with explicit movable heavy atoms and boundaries."""

import re
import shlex
from pathlib import Path

from qp.reaction.schema import HypothesisError


MODES = ("endpoint_opt", "constrained_scan", "path_search", "ts_opt", "frequency", "irc")


def _safe_keyword(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_+*/().,\-]+", value):
        raise HypothesisError("Unsafe/invalid engine keyword")
    return value


def _file(value):
    if any(c in str(value) for c in ('"', "\n", "\r")):
        raise HypothesisError("Invalid engine filename")
    return '"{}"'.format(value)


def orca_input(mode, model, settings, xyz="input.xyz", product="product.xyz",
               scans=None, hessian=None):
    if mode not in MODES:
        raise HypothesisError("Unknown reaction calculation mode")
    method = _safe_keyword(settings.get("method", "B3LYP"))
    basis = _safe_keyword(settings.get("basis", "def2-SVP"))
    extras = [_safe_keyword(v) for v in settings.get("keywords", ["D3BJ", "TightSCF"])]
    frozen = model.get("frozen_indices", [])
    if mode == "irc" and frozen:
        raise HypothesisError("Frozen-boundary IRC is not supported: use a released-boundary model or documented downhill connectivity validation")
    nprocs = settings.get("nprocs", 1)
    maxcore = settings.get("maxcore_mb", 2000)
    if type(nprocs) is not int or nprocs < 1 or type(maxcore) is not int or maxcore < 100:
        raise HypothesisError("Positive nprocs and maxcore_mb >=100 required")
    keyword = {"endpoint_opt": "Opt", "constrained_scan": "Opt",
               "path_search": "NEB-TS", "ts_opt": "OptTS",
               "frequency": "Freq" if not frozen else "NumFreq", "irc": "IRC"}[mode]
    lines = ["! {} {} {} {}".format(method, basis, " ".join(extras), keyword),
             "%pal nprocs {} end".format(nprocs), "%maxcore {}".format(maxcore)]
    epsilon = settings.get("dielectric")
    if epsilon is not None:
        if not isinstance(epsilon, (int, float)) or epsilon <= 1:
            raise HypothesisError("Dielectric must be >1")
        lines[0] += " CPCM"
        lines += ["%cpcm", "  epsilon {}".format(epsilon), "end"]
    if settings.get("pointcharges"):
        lines.append("%pointcharges {}".format(_file(settings["pointcharges"])))
    if mode in ("endpoint_opt", "constrained_scan", "path_search", "ts_opt"):
        lines += ["%geom", "  MaxIter {}".format(int(settings.get("max_opt_iterations", 150)))]
        if frozen:
            if mode != "ts_opt":
                lines.append("  CoordSys Cartesian")
            lines.append("  Constraints")
            lines.extend("    {{ C {} C }}".format(i) for i in frozen)
            lines.append("  end")
        if mode == "ts_opt":
            lines += ["  Calc_Hess true", "  Recalc_Hess {}".format(int(settings.get("recalc_hess", 5)))]
        if mode == "constrained_scan":
            if not scans or len(scans) > 3:
                raise HypothesisError("Supply one to three explicit scan coordinates")
            lines.append("  Scan")
            for scan in scans:
                pair = scan["indices"]
                if len(pair) != 2 or any(type(i) is not int or not 0 <= i < len(model["elements"]) for i in pair) or pair[0] == pair[1]:
                    raise HypothesisError("Invalid zero-based scan atom indices")
                if any(i in frozen for i in pair):
                    raise HypothesisError("Scan atoms must be movable")
                if scan["start"] <= 0 or scan["end"] <= 0 or type(scan["steps"]) is not int or scan["steps"] < 2:
                    raise HypothesisError("Invalid scan distance/steps")
                lines.append("    B {} {} = {}, {}, {}".format(*pair, scan["start"], scan["end"], scan["steps"]))
            lines.append("  end")
            if len(scans) > 1:
                if len({s["steps"] for s in scans}) != 1:
                    raise HypothesisError("Coupled scans require equal step counts")
                lines.append("  Simul_Scan true")
        lines.append("end")
    if mode == "path_search":
        lines += ["%neb", "  NEB_End_XYZFile {}".format(_file(product)),
                   "  NImages {}".format(int(settings.get("neb_images", 8))),
                   "  MaxIter {}".format(int(settings.get("max_neb_iterations", 150)))]
        if frozen:
            lines += ["  Quatern no", "  Fix_center false", "  Remove_extern_Force false"]
        lines.append("end")
    if mode == "frequency" and frozen:
        # ORCA's Partial_Hess list EXCLUDES these atoms from displacement.
        # Listing the movable atoms instead silently computes the wrong subspace.
        lines += ["%freq", "  Partial_Hess {{ {} }} end".format(" ".join(map(str, frozen))),
                  "  ProjectTR false", "  TransInvar false", "  NumHessTransInvar false", "end"]
    if mode == "irc":
        lines += ["%irc", "  Direction both", "  MaxIter {}".format(int(settings.get("max_irc_iterations", 60)))]
        if hessian:
            lines += ["  InitHess read", "  Hess_Filename {}".format(_file(hessian))]
        else:
            lines.append("  InitHess calc_anfreq")
        lines.append("end")
    if not re.fullmatch(r"[A-Za-z0-9_./\-]+", str(xyz)):
        raise HypothesisError("Use a simple whitespace-free XYZ basename")
    # Unlike %neb filenames, ORCA's xyzfile directive reads quotes literally.
    lines.append("* xyzfile {} {} {}".format(model["charge"], model["multiplicity"], xyz))
    return "\n".join(lines) + "\n"


def write_job(directory, mode, model, settings, **kwargs):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "job.inp").write_text(orca_input(mode, model, settings, **kwargs))
    executable = settings.get("orca_executable", "/mnt/nfs/vol8t/software/software/orca/6.1.1/orca")
    module = settings.get("orca_module", "orca/6.1.1")
    # Site profiles often reference unset locale variables; enable nounset only
    # after module initialization rather than breaking before ORCA is launched.
    script = ["#!/bin/bash"]
    script += ["source /etc/profile"] if module else []
    script += ["set -eo pipefail"]
    script += ["module load " + shlex.quote(module)] if module else []
    script += ["set -u", "exec {} job.inp > job.out 2>&1".format(shlex.quote(executable))]
    (directory / "run.sh").write_text("\n".join(script) + "\n")
    partition = settings.get("partition", "cpu")
    for value in (partition, settings.get("walltime", "02:00:00")):
        if not re.fullmatch(r"[A-Za-z0-9_:\-]+", value):
            raise HypothesisError("Invalid Slurm resource value")
    batch = ["#!/bin/bash", "#SBATCH --job-name=qp-reaction", "#SBATCH --partition=" + partition,
             "#SBATCH --cpus-per-task={}".format(settings.get("nprocs", 1)),
             "#SBATCH --mem={}M".format(settings.get("nprocs", 1) * settings.get("maxcore_mb", 2000) + 500),
              "#SBATCH --time=" + settings.get("walltime", "02:00:00"),
              "timeout --signal=TERM --kill-after=5 {}s bash run.sh".format(settings.get("job_timeout_seconds", 3600))]
    (directory / "submit.sh").write_text("\n".join(batch) + "\n")
