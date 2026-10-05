"""Stored job submission scripts"""

def write_qm(optimization, coord_file, basis, method, total_charge, multiplicity, guess, pcm_radii_file, constraint_freeze, dielectric, use_charge_embedding, use_implicit_solvent=True):
    """Generate a TeraChem input file (qmscript.in).

    Creates the input file content for TeraChem with the specified
    calculation parameters. Supports both single-point energy and
    geometry optimization calculations. PCM implicit solvent and
    point charge embedding can be enabled independently or together.

    Parameters
    ----------
    optimization : bool
        If True, include geometry optimization keywords.
    coord_file : str
        Name of the XYZ coordinate file.
    basis : str
        Basis set name (e.g., ``'lacvps_ecp'``).
    method : str
        DFT functional (e.g., ``'wpbeh'``, ``'ub3lyp'``).
    total_charge : int
        Total system charge.
    multiplicity : int
        Spin multiplicity (1 = singlet, 2 = doublet, etc.).
    guess : str
        Initial guess method (e.g., ``'generate'``).
    pcm_radii_file : str
        Path to PCM radii file for cavity construction.
    constraint_freeze : str
        TeraChem constraint block for frozen atoms.
    dielectric : float
        Dielectric constant for PCM solvent.
    use_charge_embedding : bool
        If True, include MM point charges from ``ptchrges.xyz``.
    use_implicit_solvent : bool, optional
        If True, include PCM implicit solvent (COSMO) block. Can be
        enabled alongside ``use_charge_embedding`` for combined
        QM/MM + implicit solvent calculations. Default is True.

    Returns
    -------
    str
        Complete TeraChem input file content.
    """

    minimization_keywords = """new_minimizer yes\nrun minimize\n""" if optimization else ""

    if use_implicit_solvent:
        pcm_section = f"""pcm cosmo
epsilon {dielectric}
pcm_radii read
pcm_radii_file {pcm_radii_file}
pcm_matrix no
"""
    else:
        pcm_section = ""

    if use_charge_embedding:
        pointcharges_section = """pointcharges ptchrges.xyz
pointcharges_self_interaction true
"""
    else:
        pointcharges_section = ""


    qmscript_content = f"""levelshift yes
levelshiftvala 0.25
levelshiftvalb 0.25
{minimization_keywords}coordinates {coord_file}
basis {basis}
method {method}
charge {total_charge}
spinmult {multiplicity}
guess {guess}
maxit 500 
dftd d3
scrdir ./scr
scf diis+a
{pcm_section}{pointcharges_section}ml_prop yes
end

{constraint_freeze}
"""
    return qmscript_content



def _orca_method_keywords(method):
    """Translate a TeraChem-style method name into ORCA simple-input keywords.

    TeraChem marks unrestricted calculations with a ``u`` prefix (e.g.,
    ``ub3lyp``). ORCA uses the separate ``UKS`` keyword instead, so the prefix
    is stripped and ``UKS`` is added. Any other method string is passed
    through unchanged.

    Parameters
    ----------
    method : str
        Method/functional name (e.g., ``'ub3lyp'``, ``'b3lyp'``, ``'r2SCAN-3c'``).

    Returns
    -------
    str
        ORCA keywords for the method (e.g., ``'UKS B3LYP'``).
    """
    if method.lower().startswith("u") and len(method) > 1 and method[1].isalpha() and not method.upper().startswith("UKS"):
        return f"UKS {method[1:].upper()}"
    return method


def _memory_to_mb(memory):
    """Convert a memory string such as ``'16G'`` or ``'500M'`` to megabytes."""
    memory = str(memory).strip().upper().rstrip("B")
    units = {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 ** 2}
    if memory and memory[-1] in units:
        return int(float(memory[:-1]) * units[memory[-1]])
    return int(float(memory))


def write_orca(optimization, coord_file, basis, method, total_charge, multiplicity, guess, dielectric, use_charge_embedding, use_implicit_solvent=True, nprocs=16, memory="64G", dispersion="D3BJ", aux_basis="def2/J", extra_keywords=None):
    """Generate an ORCA input file (qmscript.inp).

    Builds an ORCA input mirroring the TeraChem settings used by
    :func:`write_qm`: the same functional, charge, multiplicity, dispersion,
    level shifting, implicit solvent (CPCM), and point-charge embedding. For
    geometry optimizations, only hydrogen positions are relaxed, matching the
    TeraChem heavy-atom freeze.

    Parameters
    ----------
    optimization : bool
        If True, optimize hydrogen positions (heavy atoms fixed).
    coord_file : str
        Name of the XYZ coordinate file.
    basis : str
        ORCA basis set name (e.g., ``'def2-SVP'``).
    method : str
        DFT functional. A leading ``u`` (TeraChem convention) becomes ``UKS``.
    total_charge : int
        Total system charge.
    multiplicity : int
        Spin multiplicity.
    guess : str
        Initial guess. ``'generate'`` (TeraChem default) uses the ORCA
        default; ORCA guess names (``PModel``, ``Hueckel``, ``HCore``,
        ``PAtom``) are passed through.
    dielectric : float
        Dielectric constant for CPCM.
    use_charge_embedding : bool
        If True, read MM point charges from ``ptchrges.pc``.
    use_implicit_solvent : bool, optional
        If True, include CPCM implicit solvent. Default is True.
    nprocs : int, optional
        Number of MPI processes for ORCA (``%pal``). Default is 16.
    memory : str, optional
        Total job memory (e.g., ``'64G'``). 75% is divided across processes
        for ``%maxcore`` to leave headroom. Default is ``'64G'``.
    dispersion : str or None, optional
        Dispersion correction keyword (e.g., ``'D3BJ'``, ``'D3ZERO'``,
        ``'D4'``). Use None to disable. Default is ``'D3BJ'``.
    aux_basis : str or None, optional
        Auxiliary basis for RIJCOSX. Use None to disable RIJCOSX.
        Default is ``'def2/J'``.
    extra_keywords : str, list of str, or None, optional
        Additional ORCA simple-input keywords appended to the ``!`` line
        (e.g., ``'EnGrad'`` to also compute nuclear gradients, or
        ``'TightSCF'``). Default is None.

    Returns
    -------
    str
        Complete ORCA input file content.
    """
    keywords = [_orca_method_keywords(method)]
    if dispersion:
        keywords.append(dispersion)
    keywords.append(basis)
    if aux_basis:
        keywords += [aux_basis, "RIJCOSX"]
    if use_implicit_solvent:
        keywords.append("CPCM")
    keywords.append("SlowConv")
    if optimization:
        keywords.append("Opt")
    if extra_keywords:
        if isinstance(extra_keywords, str):
            extra_keywords = extra_keywords.split()
        keywords += [k for k in extra_keywords if k not in keywords]

    maxcore = max(_memory_to_mb(memory) * 3 // (4 * nprocs), 500)

    orca_guesses = {"pmodel": "PModel", "hueckel": "Hueckel", "hcore": "HCore", "patom": "PAtom"}
    guess_line = f"  Guess {orca_guesses[guess.lower()]}\n" if guess and guess.lower() in orca_guesses else ""

    blocks = f"""%pal
  nprocs {nprocs}
end
%maxcore {maxcore}
%scf
  MaxIter 500
{guess_line}  Shift Shift 0.25 ErrOff 0.10 end
end
"""
    if use_implicit_solvent:
        blocks += f"""%cpcm
  epsilon {dielectric}
end
"""
    if use_charge_embedding:
        blocks += '%pointcharges "ptchrges.pc"\n'
    if optimization:
        blocks += """%geom
  optimizehydrogens true
end
"""

    return f"""! {' '.join(keywords)}
{blocks}
* xyzfile {total_charge} {multiplicity} {coord_file}
"""


def write_slurm_orca_job(job_name, nprocs, memory, partition="cpu", account=None, time_limit=None, orca_module="orca/6.1.1", structure_name=None):
    """Generate a SLURM submission script for ORCA.

    Requests ``nprocs`` tasks on a single node, loads the ORCA module, and
    calls ORCA by its full path (required for its internal MPI launch).
    ORCA runs in node-local scratch (``$TMPDIR``, falling back to ``/tmp``)
    and only results (``qmscript.out``, ``.gbw``, ``.engrad``, property
    files) are copied back, so large temporary files never reach the job
    directory. A run counts as successful only if ORCA prints its
    normal-termination banner. On success, the wavefunction is converted to
    Molden format at ``scr/<structure_name>.molden`` so that ``qp analyze``
    can find it.

    Parameters
    ----------
    job_name : str
        Name for the SLURM job.
    nprocs : int
        Number of MPI processes (SLURM tasks).
    memory : str
        Total memory for the job (e.g., ``'64G'``).
    partition : str, optional
        SLURM partition. Default is ``'cpu'``.
    account : str or None, optional
        SLURM account to charge.
    time_limit : str or None, optional
        Wall time limit (e.g., ``'2-00:00:00'``).
    orca_module : str or None, optional
        Environment module that provides ORCA. Default is ``'orca/6.1.1'``.
    structure_name : str or None, optional
        Cluster name used for the Molden file. Defaults to ``job_name``.

    Returns
    -------
    str
        Complete SLURM submission script content.
    """
    structure_name = structure_name or job_name
    account_line = f"#SBATCH --account={account}\n" if account else ""
    time_line = f"#SBATCH --time={time_limit}\n" if time_limit else ""
    module_line = f"module load {orca_module}\n" if orca_module else ""

    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
{account_line}#SBATCH --nodes=1
#SBATCH --ntasks={nprocs}
#SBATCH --cpus-per-task=1
#SBATCH --mem={memory}
{time_line}#SBATCH --output={job_name}.log

source /etc/profile
{module_line}
# ORCA must be called by its full path for parallel runs
ORCA_EXE="${{ORCA_BIN:-$(command -v orca)}}"

# Run in node-local scratch so ORCA's large temporary files never land on
# the (shared, quota-limited) job directory; copy back only the results.
JOBDIR="${{SLURM_SUBMIT_DIR:-$(pwd)}}"
WORK="${{TMPDIR:-/tmp}}/orca_${{SLURM_JOB_ID:-$$}}"
mkdir -p "$WORK"
trap 'cd "$JOBDIR"; rm -rf "$WORK"' EXIT
cp "$JOBDIR"/qmscript.inp "$JOBDIR"/*.xyz "$WORK"/
[ -f "$JOBDIR/ptchrges.pc" ] && cp "$JOBDIR/ptchrges.pc" "$WORK"/
cd "$WORK"

echo "Run Start Time: $(date '+%Y-%m-%d %H:%M:%S') on $(hostname)" >> "$JOBDIR/.submit_record"
"$ORCA_EXE" qmscript.inp > "$JOBDIR/qmscript.out"
status=$?
# ORCA can exit 0 after an error termination; trust only the normal-termination banner
if [ $status -eq 0 ] && ! grep -q "ORCA TERMINATED NORMALLY" "$JOBDIR/qmscript.out"; then
    status=1
fi
for f in qmscript.gbw qmscript.engrad qmscript.property.txt qmscript_property.txt qmscript.xyz qmscript_trj.xyz; do
    [ -f "$f" ] && cp "$f" "$JOBDIR"/
done
if [ $status -eq 0 ] && [ -f qmscript.gbw ]; then
    mkdir -p "$JOBDIR/scr"
    orca_2mkl qmscript -molden > /dev/null && mv qmscript.molden.input "$JOBDIR/scr/{structure_name}.molden"
fi
echo "Run End Time: $(date '+%Y-%m-%d %H:%M:%S')" >> "$JOBDIR/.submit_record"
exit $status
"""


def write_slurm_job(job_name, gpus, memory, partition="xeon-g6-volta", account=None, time_limit=None):
    """Generate a SLURM submission script for TeraChem.

    Creates a bash script with SLURM directives for GPU job submission.
    Configured for systems with NVIDIA Volta GPUs and the TeraChem module.

    Parameters
    ----------
    job_name : str
        Name for the SLURM job (used in output filenames).
    gpus : int
        Number of GPUs to request.
    memory : str
        Memory allocation (currently unused but kept for API consistency).
    partition : str, optional
        SLURM partition. Default is ``'xeon-g6-volta'``.
    account : str or None, optional
        SLURM account to charge.
    time_limit : str or None, optional
        Wall time limit (e.g., ``'2-00:00:00'``).

    Returns
    -------
    str
        Complete SLURM submission script content.
    """
    account_line = f"#SBATCH --account={account}\n" if account else ""
    time_line = f"#SBATCH --time={time_limit}\n" if time_limit else ""

    jobscript_content = f"""#! /bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
{account_line}{time_line}#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:{gpus}
#SBATCH --cpus-per-task={gpus * 20}
#SBATCH --output={job_name}.log

source /etc/profile

#---TC setup---
module load terachem/1.9-2023.11-dev

# your command to run terachem
echo "Run Start Time: $(date '+%Y-%m-%d %H:%M:%S')" >> .submit_record
terachem qmscript.in > qmscript.out
echo "Run End Time: $(date '+%Y-%m-%d %H:%M:%S')" >> .submit_record
"""

    return jobscript_content

def write_sge_job(job_name, gpus, memory):
    """Generate a Sun Grid Engine (SGE) submission script for TeraChem.

    Creates a bash script with SGE directives for GPU job submission.
    Includes a minimum runtime enforcement (10 minutes) to prevent
    scheduler issues with very short jobs.

    Parameters
    ----------
    job_name : str
        Name for the SGE job (prefixed with 'Z' in the script).
    gpus : int
        Number of GPUs/parallel threads to request.
    memory : str
        Memory allocation string (e.g., ``'8G'``).

    Returns
    -------
    str
        Complete SGE submission script content.
    """

    jobscript_content = f"""#!/bin/bash
#$ -N Z{job_name}
#$ -cwd
#$ -l h_rt=300:00:00
#$ -l h_rss={memory}
#$ -q (gpusnew|gpus|gpusbig)
#$ -l gpus=1
#$ -pe smp {gpus}
# -fin qmscript.in
# -fin *.xyz
# -fout scr/

module load cuda/10.0
module load terachem/071920-cuda10.0-intel16
module load intel/16.0.109

export OMP_NUM_THREADS={gpus}

# Start time
SECONDS=0

echo "Run Start Time: $(date '+%Y-%m-%d %H:%M:%S')" >> .submit_record
terachem qmscript.in > $SGE_O_WORKDIR/qmscript.out
echo "Run End Time: $(date '+%Y-%m-%d %H:%M:%S')" >> .submit_record

# Calculate elapsed time in seconds
ELAPSED_TIME=$SECONDS

# If elapsed time is less than 600 seconds (10 minutes), sleep for the remainder
# This is just for Gibraltar which can't handle short jobs
if [ $ELAPSED_TIME -lt 600 ]; then
    SLEEP_TIME=$((600 - ELAPSED_TIME))
    sleep $SLEEP_TIME
fi
"""

    return jobscript_content