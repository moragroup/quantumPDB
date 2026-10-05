"""
Unit and regression test for the qp package.
"""

# Import package, test suite, and other packages as needed
import sys
import os
import json
import tempfile
from unittest import mock

import pytest

import qp
from qp.manager.charge_embedding import load_custom_charges, get_charges, parse_pdb_to_xyz
from qp.manager.job_scripts import write_qm, write_orca, write_slurm_orca_job, write_slurm_job


def test_qp_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "qp" in sys.modules


def test_load_custom_charges_valid_json():
    """Test that load_custom_charges correctly loads a valid JSON file."""
    charges = {
        "ALA": {"N": -0.4157, "H": 0.2719},
        "GLY": {"N": -0.4157}
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(charges, f)
        tmppath = f.name
    try:
        result = load_custom_charges(tmppath)
        assert result == charges
        assert result["ALA"]["N"] == -0.4157
        assert result["ALA"]["H"] == 0.2719
        assert result["GLY"]["N"] == -0.4157
    finally:
        os.unlink(tmppath)


def test_load_custom_charges_file_not_found():
    """Test that load_custom_charges raises FileNotFoundError for missing files."""
    with pytest.raises(FileNotFoundError):
        load_custom_charges("/nonexistent/path.json")


def test_load_custom_charges_invalid_json():
    """Test that load_custom_charges raises an error for invalid JSON."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        f.write("not valid json {{")
        tmppath = f.name
    try:
        with pytest.raises(json.JSONDecodeError):
            load_custom_charges(tmppath)
    finally:
        os.unlink(tmppath)


def test_get_charges_default_uses_ff14sb():
    """Test that get_charges uses ff14SB dict when no custom file is provided."""
    with mock.patch('qp.manager.charge_embedding.ff14SB_dict.get_ff14SB_dict') as mock_ff, \
         mock.patch('qp.manager.charge_embedding.rename_and_clean_resnames'), \
         mock.patch('qp.manager.charge_embedding.parse_pdb'), \
         mock.patch('qp.manager.charge_embedding.remove_qm_atoms'), \
         mock.patch('qp.manager.charge_embedding.read_xyz', return_value=__import__('numpy').array([[0, 0, 0]])), \
         mock.patch('qp.manager.charge_embedding.parse_pdb_to_xyz'), \
         mock.patch('os.path.exists', return_value=False), \
         mock.patch('os.mkdir'), \
         mock.patch('os.getcwd', return_value='/fake/output/pdb1/A200/method'), \
         mock.patch('shutil.rmtree'):
        mock_ff.return_value = {"ALA": {"N": -0.4157}}
        get_charges(20)
        mock_ff.assert_called_once()


def test_get_charges_custom_file_used():
    """Test that get_charges uses a custom charges file when provided."""
    charges = {"ALA": {"N": -0.4157, "H": 0.2719}}
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(charges, f)
        tmppath = f.name

    _real_exists = os.path.exists

    def _exists_side_effect(path):
        """Allow the custom charges file to exist, mock others as absent."""
        if path == tmppath:
            return _real_exists(path)
        return False

    try:
        with mock.patch('qp.manager.charge_embedding.ff14SB_dict.get_ff14SB_dict') as mock_ff, \
             mock.patch('qp.manager.charge_embedding.rename_and_clean_resnames'), \
             mock.patch('qp.manager.charge_embedding.parse_pdb'), \
             mock.patch('qp.manager.charge_embedding.remove_qm_atoms'), \
             mock.patch('qp.manager.charge_embedding.read_xyz', return_value=__import__('numpy').array([[0, 0, 0]])), \
             mock.patch('qp.manager.charge_embedding.parse_pdb_to_xyz'), \
             mock.patch('os.path.exists', side_effect=_exists_side_effect), \
             mock.patch('os.mkdir'), \
             mock.patch('os.getcwd', return_value='/fake/output/pdb1/A200/method'), \
             mock.patch('shutil.rmtree'):
            get_charges(20, charge_embedding_charges=tmppath)
            mock_ff.assert_not_called()
    finally:
        os.unlink(tmppath)


def test_parse_pdb_to_xyz_residue_based_selection():
    """Test that parse_pdb_to_xyz includes complete residues.

    If any atom of a residue is within the cutoff distance, all atoms
    of that residue must be included in the output.  Residues entirely
    outside the cutoff must be excluded.
    """
    import numpy as np

    def make_pdb_line(serial, atom_name, res_name, chain, res_seq, x, y, z, charge):
        """Build a fixed-width PDB ATOM line."""
        return (
            f"ATOM  {serial:>5d} {atom_name:<4s} {res_name:>3s} {chain}{res_seq:>4d}"
            f"    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"  1.00{charge:>6.2f}           \n"
        )

    # ALA atom 1 (N)  at (0,0,0)  -> distance 0  <= 10 -> triggers residue inclusion
    # ALA atom 2 (CA) at (12,0,0) -> distance 12 > 10  -> excluded by old code, included by new
    # GLY atom   (N)  at (50,0,0) -> distance 50 > 10  -> residue excluded entirely
    pdb_content = (
        make_pdb_line(1, "N",  "ALA", "A", 1,  0.0, 0.0, 0.0, -0.42)
        + make_pdb_line(2, "CA", "ALA", "A", 1, 12.0, 0.0, 0.0,  0.42)
        + make_pdb_line(3, "N",  "GLY", "A", 2, 50.0, 0.0, 0.0, -0.30)
    )

    qm_centroid = np.array([0.0, 0.0, 0.0])
    cutoff = 10.0

    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as pdb_f:
        pdb_f.write(pdb_content)
        pdb_path = pdb_f.name

    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as out_f:
        out_path = out_f.name

    try:
        parse_pdb_to_xyz(pdb_path, out_path, qm_centroid, cutoff)

        with open(out_path, 'r') as f:
            result_lines = f.readlines()

        # Header: atom count and comment
        assert result_lines[0].strip() == "2", (
            "Expected 2 atoms (full ALA residue), got: " + result_lines[0].strip()
        )
        assert result_lines[1].strip() == "Generated from PDB file"

        # Both ALA atoms present, GLY excluded
        data_lines = result_lines[2:]
        assert len(data_lines) == 2, f"Expected 2 data lines, got {len(data_lines)}"

        # First atom: charge=-0.42, coords=(0,0,0)
        parts_0 = data_lines[0].split()
        assert float(parts_0[0]) == pytest.approx(-0.42)
        assert float(parts_0[1]) == pytest.approx(0.0)

        # Second atom: charge=0.42, coords=(12,0,0)
        parts_1 = data_lines[1].split()
        assert float(parts_1[0]) == pytest.approx(0.42)
        assert float(parts_1[1]) == pytest.approx(12.0)
    finally:
        os.unlink(pdb_path)
        os.unlink(out_path)


# --- Common kwargs for write_qm tests ---

_BASE_QM_KWARGS = dict(
    optimization=False,
    coord_file="cluster.xyz",
    basis="lacvps_ecp",
    method="wpbeh",
    total_charge=-1,
    multiplicity=1,
    guess="generate",
    pcm_radii_file="/path/to/pcm_radii",
    constraint_freeze="",
    dielectric=10,
)


def test_write_qm_charge_embedding_only():
    """Test that charge embedding without implicit solvent produces pointcharges but no PCM."""
    result = write_qm(**_BASE_QM_KWARGS, use_charge_embedding=True, use_implicit_solvent=False)
    assert "pointcharges ptchrges.xyz" in result
    assert "pointcharges_self_interaction true" in result
    assert "pcm cosmo" not in result
    assert "epsilon" not in result


def test_write_qm_implicit_solvent_only():
    """Test that implicit solvent without charge embedding produces PCM but no pointcharges."""
    result = write_qm(**_BASE_QM_KWARGS, use_charge_embedding=False, use_implicit_solvent=True)
    assert "pcm cosmo" in result
    assert "epsilon 10" in result
    assert "pcm_radii_file /path/to/pcm_radii" in result
    assert "pointcharges ptchrges.xyz" not in result


def test_write_qm_both_enabled():
    """Test that both charge embedding and implicit solvent can be enabled together."""
    result = write_qm(**_BASE_QM_KWARGS, use_charge_embedding=True, use_implicit_solvent=True)
    assert "pointcharges ptchrges.xyz" in result
    assert "pointcharges_self_interaction true" in result
    assert "pcm cosmo" in result
    assert "epsilon 10" in result
    assert "pcm_radii_file /path/to/pcm_radii" in result


def test_write_qm_neither_enabled():
    """Test that disabling both produces neither PCM nor pointcharges blocks."""
    result = write_qm(**_BASE_QM_KWARGS, use_charge_embedding=False, use_implicit_solvent=False)
    assert "pointcharges ptchrges.xyz" not in result
    assert "pcm cosmo" not in result
    # Core keywords should still be present
    assert "method wpbeh" in result
    assert "basis lacvps_ecp" in result


def test_write_qm_backward_compatible_default():
    """Test that omitting use_implicit_solvent defaults to True (PCM on).

    This preserves backward compatibility for callers that do not pass
    the new parameter.
    """
    result = write_qm(**_BASE_QM_KWARGS, use_charge_embedding=False)
    assert "pcm cosmo" in result
    assert "pointcharges ptchrges.xyz" not in result


_BASE_ORCA_KWARGS = dict(
    optimization=False,
    coord_file="A302.xyz",
    basis="def2-SVP",
    method="ub3lyp",
    total_charge=3,
    multiplicity=5,
    guess="generate",
    dielectric=10,
    nprocs=32,
    memory="128G",
)


def test_write_orca_implicit_solvent():
    """Test ORCA input with CPCM implicit solvent and unrestricted functional."""
    result = write_orca(**_BASE_ORCA_KWARGS, use_charge_embedding=False, use_implicit_solvent=True)
    keywords = result.splitlines()[0]
    assert keywords.startswith("! UKS B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM")
    assert "Opt" not in keywords
    assert "nprocs 32" in result
    assert "%maxcore 3072" in result
    assert "epsilon 10" in result
    assert "%pointcharges" not in result
    assert "* xyzfile 3 5 A302.xyz" in result


def test_write_orca_charge_embedding_only():
    """Test ORCA input with point charges and no implicit solvent."""
    result = write_orca(**_BASE_ORCA_KWARGS, use_charge_embedding=True, use_implicit_solvent=False)
    assert '%pointcharges "ptchrges.pc"' in result
    assert "CPCM" not in result
    assert "%cpcm" not in result


def test_write_orca_optimization_hydrogens_only():
    """Test ORCA optimization relaxes only hydrogens, mirroring TeraChem heavy-atom freeze."""
    kwargs = dict(_BASE_ORCA_KWARGS, optimization=True)
    result = write_orca(**kwargs, use_charge_embedding=False)
    assert " Opt" in result.splitlines()[0]
    assert "optimizehydrogens true" in result


def test_write_orca_restricted_method_passthrough():
    """Test that methods without a 'u' prefix are passed through unchanged."""
    kwargs = dict(_BASE_ORCA_KWARGS, method="r2SCAN-3c", dispersion=None, aux_basis=None, basis="")
    result = write_orca(**kwargs, use_charge_embedding=False, use_implicit_solvent=False)
    assert result.splitlines()[0].split() == ["!", "r2SCAN-3c", "SlowConv"]


def test_write_orca_extra_keywords():
    """Test that extra ORCA keywords (string or list) are appended once."""
    as_str = write_orca(**_BASE_ORCA_KWARGS, use_charge_embedding=False, extra_keywords="EnGrad TightSCF")
    as_list = write_orca(**_BASE_ORCA_KWARGS, use_charge_embedding=False, extra_keywords=["EnGrad", "TightSCF", "CPCM"])
    assert as_str.splitlines()[0].endswith("SlowConv EnGrad TightSCF")
    assert as_list.splitlines()[0] == as_str.splitlines()[0]


def test_write_slurm_orca_job():
    """Test ORCA SLURM script requests MPI tasks and runs ORCA by full path."""
    result = write_slurm_orca_job("6nieA302", 32, "128G", partition="cpu", account="mora", time_limit="1-00:00:00", structure_name="A302")
    assert "#SBATCH --partition=cpu" in result
    assert "#SBATCH --account=mora" in result
    assert "#SBATCH --ntasks=32" in result
    assert "#SBATCH --mem=128G" in result
    assert "#SBATCH --time=1-00:00:00" in result
    assert "module load orca/6.1.1" in result
    assert '"$ORCA_EXE" qmscript.inp > "$JOBDIR/qmscript.out"' in result
    assert "scr/A302.molden" in result
    assert 'grep -q "ORCA TERMINATED NORMALLY"' in result
    # ORCA runs in node-local scratch and results are copied back
    assert 'WORK="${TMPDIR:-/tmp}/orca_${SLURM_JOB_ID:-$$}"' in result
    assert 'cp "$f" "$JOBDIR"/' in result
    assert "mpirun" not in result


def test_write_slurm_job_defaults_unchanged():
    """Test TeraChem SLURM script keeps its original partition when no overrides are given."""
    result = write_slurm_job("6nieA302", 1, "8G")
    assert "#SBATCH --partition=xeon-g6-volta" in result
    assert "--account" not in result
    assert "--time" not in result
