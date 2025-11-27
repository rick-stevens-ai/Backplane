#!/usr/bin/env python3
"""Run methane single-point energy tests across multiple simulation codes.

Each backend gets its own scratch directory under APPS/test_runs/methane/<code>.
The script prepares tiny input decks, invokes the relevant executable, and
parses the reported total potential energy. Results are printed at the end.

All tests assume the repositories were cloned under ./APPS as set up earlier.
Quantum ESPRESSO/QE relies on existing pseudopotentials in APPS/q-e/pseudo.
CP2K/GROMACS/LAMMPS/GPAW depend on their usual Homebrew locations.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional

REPO_ROOT = Path(__file__).resolve().parent
APPS_DIR = REPO_ROOT / "APPS"
WORK_ROOT = APPS_DIR / "test_runs" / "methane"
WORK_ROOT.mkdir(parents=True, exist_ok=True)

RESULTS: Dict[str, Dict[str, Optional[float]]] = {}

RYDBERG_TO_EV = 13.605693009
HARTREE_TO_EV = 27.211386018
KCALMOL_TO_EV = 0.0433641153087705
KJMOL_TO_EV = 0.0103642723013312


def run(cmd, workdir: Path, env: Optional[dict] = None, timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    """Wrapper around subprocess.run that raises helpful errors."""
    default_env = os.environ.copy()
    if env:
        default_env.update(env)
    return subprocess.run(
        cmd,
        cwd=str(workdir),
        env=default_env,
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )


def record_result(tool: str, success: bool, energy: Optional[float], log: str, notes: str = "", unit: Optional[str] = None) -> None:
    RESULTS[tool] = {
        "success": success,
        "energy": energy,
        "unit": unit,
        "log": log,
        "notes": notes,
    }


# ---------------- Quantum ESPRESSO ---------------- #

def run_quantum_espresso() -> None:
    tool = "quantum_espresso"
    bin_path = APPS_DIR / "q-e" / "bin" / "pw.x"
    pseudo_dir = APPS_DIR / "q-e" / "pseudo"
    if not bin_path.exists():
        record_result(tool, False, None, "", "pw.x not found at expected path")
        return
    pseudo_c = pseudo_dir / "C.UPF"
    pseudo_h = pseudo_dir / "H_US.van"
    if not (pseudo_c.exists() and pseudo_h.exists()):
        record_result(tool, False, None, "", "Missing C.UPF or H_US.van pseudopotential")
        return
    workdir = WORK_ROOT / tool
    workdir.mkdir(parents=True, exist_ok=True)
    tmp_dir = workdir / "tmp"
    tmp_dir.mkdir(exist_ok=True)
    input_path = workdir / "ch4.scf.in"
    output_path = workdir / "pw.out"
    input_path.write_text(
        f"""
&control
    calculation='scf',
    prefix='ch4',
    pseudo_dir='{pseudo_dir}',
    outdir='{tmp_dir}',
    tprnfor=.false.,
    tstress=.false.,
/
&system
    ibrav=1, celldm(1)=20.0, nat=5, ntyp=2,
    ecutwfc=50.0,
    occupations='smearing', smearing='gaussian', degauss=0.02,
    input_dft='PBE'
/
&electrons
    conv_thr=1.0d-7,
    mixing_beta=0.7,
    diagonalization='david'
/
ATOMIC_SPECIES
 C 12.011 C.UPF
 H 1.008 H_US.van
ATOMIC_POSITIONS (angstrom)
 C 0.000 0.000 0.000
 H 0.629 0.629 0.629
 H -0.629 -0.629 0.629
 H -0.629 0.629 -0.629
 H 0.629 -0.629 -0.629
K_POINTS gamma
""".strip()
    )
    env = {"OMPI_MCA_btl": "self", "ESPRESSO_TMPDIR": str(tmp_dir)}
    result = run([str(bin_path), "-in", input_path.name], workdir, env)
    output_path.write_text(result.stdout + "\n" + result.stderr)
    if result.returncode != 0:
        record_result(tool, False, None, str(output_path), "pw.x exited with non-zero status")
        return
    energy = None
    for line in result.stdout.splitlines():
        if line.strip().startswith("!") and "total energy" in line:
            # Format: !    total energy              =   -xxx Ry
            parts = line.split("=", 1)
            if len(parts) == 2:
                try:
                    energy = float(parts[1].split()[0])
                except ValueError:
                    pass
    if energy is None:
        record_result(tool, False, None, str(output_path), "total energy line not found")
    else:
        energy_ev = energy * RYDBERG_TO_EV
        record_result(tool, True, energy_ev, str(output_path), unit="eV")


# ---------------- CP2K ---------------- #

def run_cp2k() -> None:
    tool = "cp2k"
    exe = shutil.which("cp2k.ssmp")
    if not exe:
        record_result(tool, False, None, "", "cp2k.ssmp not found in PATH")
        return
    workdir = WORK_ROOT / tool
    workdir.mkdir(parents=True, exist_ok=True)
    input_path = workdir / "ch4_cp2k.inp"
    log_path = workdir / "cp2k.out"
    input_path.write_text(
        """
&GLOBAL
  PROJECT CH4
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 400
    &END MGRID
    &SCF
      MAX_SCF 50
      EPS_SCF 1.0E-6
      CHOLESKY OFF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 15.0 15.0 15.0
    &END CELL
    &COORD
      C    0.000    0.000    0.000
      H    0.629    0.629    0.629
      H   -0.629   -0.629    0.629
      H   -0.629    0.629   -0.629
      H    0.629   -0.629   -0.629
    &END COORD
    &KIND C
      ELEMENT C
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND H
      ELEMENT H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
""".strip()
    )
    env = {"CP2K_DATA_DIR": "/opt/homebrew/opt/cp2k/share/cp2k/data"}
    result = run([exe, "-i", input_path.name], workdir, env)
    log_path.write_text(result.stdout + "\n" + result.stderr)
    if result.returncode != 0:
        record_result(tool, False, None, str(log_path), "cp2k run failed")
        return
    energy = None
    for line in result.stdout.splitlines():
        line = line.strip()
        if line.startswith("ENERGY| Total FORCE_EVAL"):
            try:
                energy = float(line.rsplit(None, 1)[-1])
            except ValueError:
                pass
    if energy is None:
        record_result(tool, False, None, str(log_path), "ENERGY line not found")
    else:
        energy_ev = energy * HARTREE_TO_EV
        record_result(tool, True, energy_ev, str(log_path), unit="eV")


# ---------------- GPAW ---------------- #

def run_gpaw() -> None:
    tool = "gpaw"
    python_bin = shutil.which("python3.14") or shutil.which("python3")
    gpaw_script = WORK_ROOT / tool / "run_gpaw.py"
    workdir = gpaw_script.parent
    workdir.mkdir(parents=True, exist_ok=True)
    script_body = """
from ase.build import molecule
from gpaw import GPAW, FermiDirac
atoms = molecule('CH4')
atoms.center(vacuum=4.0)
calc = GPAW(mode='lcao', basis='dzp', xc='PBE', occupations=FermiDirac(0.1), txt='gpaw.log')
atoms.calc = calc
energy = atoms.get_potential_energy()
print(f'TOTAL_ENERGY={energy:.8f} eV')
""".strip()
    gpaw_script.write_text(script_body)
    env = {"OMPI_MCA_btl": "self"}
    result = run([python_bin, gpaw_script.name], workdir, env)
    log_path = workdir / "gpaw_stdout.log"
    log_path.write_text(result.stdout + "\n" + result.stderr)
    if result.returncode != 0:
        record_result(tool, False, None, str(log_path), "gpaw Python run failed")
        return
    energy = None
    for line in result.stdout.splitlines():
        if line.startswith("TOTAL_ENERGY="):
            try:
                energy = float(line.split("=", 1)[1].split()[0])
            except ValueError:
                pass
    if energy is None:
        record_result(tool, False, None, str(log_path), "TOTAL_ENERGY marker missing")
    else:
        record_result(tool, True, energy, str(log_path), unit="eV")


# ---------------- LAMMPS ---------------- #

def run_lammps() -> None:
    tool = "lammps"
    exe = shutil.which("lmp_serial") or "/opt/homebrew/opt/lammps/bin/lmp_serial"
    if not Path(exe).exists():
        record_result(tool, False, None, "", "lmp_serial not found")
        return
    workdir = WORK_ROOT / tool
    workdir.mkdir(parents=True, exist_ok=True)
    data_path = workdir / "data.methane"
    input_path = workdir / "in.methane"
    log_path = workdir / "log.lammps"
    data_path.write_text(
        """
LAMMPS methane

5 atoms
4 bonds
6 angles

2 atom types
1 bond types
1 angle types

-5.0 5.0 xlo xhi
-5.0 5.0 ylo yhi
-5.0 5.0 zlo zhi

Masses

1 12.011
2 1.008

Atoms

1 1 1 -0.4 0.000 0.000 0.000
2 1 2 0.1  0.629 0.629 0.629
3 1 2 0.1 -0.629 -0.629 0.629
4 1 2 0.1 -0.629 0.629 -0.629
5 1 2 0.1  0.629 -0.629 -0.629

Bonds

1 1 1 2
2 1 1 3
3 1 1 4
4 1 1 5

Angles

1 1 2 1 3
2 1 2 1 4
3 1 2 1 5
4 1 3 1 4
5 1 3 1 5
6 1 4 1 5
""".strip()
    )
    input_path.write_text(
        """
units real
atom_style full
boundary p p p
read_data data.methane

pair_style lj/cut 12.0
pair_modify mix arithmetic
pair_coeff 1 1 0.1094 3.5
pair_coeff 2 2 0.0157 2.5
pair_coeff 1 2 0.04 3.0

bond_style harmonic
bond_coeff 1 340.0 1.09

angle_style harmonic
angle_coeff 1 33.0 109.5

special_bonds lj/coul 0.0 0.0 0.5
thermo_style custom step pe ke etotal temp
thermo 1
run 0
""".strip()
    )
    result = run([exe, "-in", input_path.name], workdir)
    log_path.write_text(result.stdout + "\n" + result.stderr)
    if result.returncode != 0:
        record_result(tool, False, None, str(log_path), "LAMMPS run failed")
        return
    energy = None
    for line in result.stdout.splitlines():
        if line.strip().startswith("0 ") and len(line.split()) >= 4:
            # step pe ke etotal temp
            try:
                energy = float(line.split()[1])
            except ValueError:
                pass
    if energy is None:
        record_result(tool, False, None, str(log_path), "thermo output missing")
    else:
        energy_ev = energy * KCALMOL_TO_EV
        record_result(tool, True, energy_ev, str(log_path), unit="eV")


# ---------------- GROMACS ---------------- #

def run_gromacs() -> None:
    tool = "gromacs"
    exe = shutil.which("gmx")
    if not exe:
        record_result(tool, False, None, "", "gmx command not found")
        return
    gmxlib = Path("/opt/homebrew/opt/gromacs/share/gromacs/top")
    if not gmxlib.exists():
        record_result(tool, False, None, "", "GROMACS force-field directory not found")
        return
    workdir = WORK_ROOT / tool
    workdir.mkdir(parents=True, exist_ok=True)
    gro_path = workdir / "methane.gro"
    top_path = workdir / "methane.top"
    mdp_path = workdir / "methane.mdp"
    gro_path.write_text(
        """Methane
5
    1METH    C     1   0.000   0.000   0.000
    1METH    H1    2   0.068   0.068   0.068
    1METH    H2    3  -0.068  -0.068   0.068
    1METH    H3    4  -0.068   0.068  -0.068
    1METH    H4    5   0.068  -0.068  -0.068
   3.0   3.0   3.0
""".strip()
    )
    top_path.write_text(
        """
#include "oplsaa.ff/forcefield.itp"

[ moleculetype ]
METH 3

[ atoms ]
1 opls_135 1 METH C 1 -0.24 12.011
2 opls_140 1 METH H1 1 0.06 1.008
3 opls_140 1 METH H2 1 0.06 1.008
4 opls_140 1 METH H3 1 0.06 1.008
5 opls_140 1 METH H4 1 0.06 1.008

[ bonds ]
1 2 1 0.109 284512
1 3 1 0.109 284512
1 4 1 0.109 284512
1 5 1 0.109 284512

[ angles ]
2 1 3 1 109.5 376.56
2 1 4 1 109.5 376.56
2 1 5 1 109.5 376.56
3 1 4 1 109.5 376.56
3 1 5 1 109.5 376.56
4 1 5 1 109.5 376.56

[ system ]
Methane

[ molecules ]
METH 1
""".strip()
    )
    mdp_path.write_text(
        """
integrator       = steep
nsteps           = 0
cutoff-scheme    = verlet
coulombtype      = Cut-off
rcoulomb         = 0.5
rvdw             = 0.5
rlist            = 0.5
constraints      = none
ld_seed          = 12345
""".strip()
    )
    env = {"GMXLIB": str(gmxlib)}
    grompp = run([exe, "grompp", "-f", mdp_path.name, "-c", gro_path.name, "-p", top_path.name, "-o", "methane.tpr"], workdir, env)
    (workdir / "grompp.log").write_text(grompp.stdout + "\n" + grompp.stderr)
    if grompp.returncode != 0:
        record_result(tool, False, None, str(workdir / "grompp.log"), "grompp failed")
        return
    mdrun = run([exe, "mdrun", "-s", "methane.tpr", "-ntmpi", "1", "-ntomp", "1", "-nsteps", "0", "-deffnm", "methane"], workdir, env)
    mdlog = workdir / "methane.log"
    mdlog.write_text(mdrun.stdout + "\n" + mdrun.stderr)
    if mdrun.returncode != 0:
        record_result(tool, False, None, str(mdlog), "mdrun failed")
        return
    energy = None
    for line in (mdrun.stdout + "\n" + mdrun.stderr).splitlines():
        if "Potential Energy" in line:
            try:
                energy = float(line.split("=", 1)[1].split()[0])
            except ValueError:
                pass
    if energy is None:
        record_result(tool, False, None, str(mdlog), "Potential Energy not parsed")
    else:
        energy_ev = energy * KJMOL_TO_EV
        record_result(tool, True, energy_ev, str(mdlog), unit="eV")


# ---------------- main ---------------- #

def main() -> None:
    run_quantum_espresso()
    run_cp2k()
    run_gpaw()
    run_lammps()
    run_gromacs()
    print(json.dumps(RESULTS, indent=2))


if __name__ == "__main__":
    main()
