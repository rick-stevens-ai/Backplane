#!/usr/bin/env python3
"""
CP2K wrapper for mixed Gaussian/plane-wave DFT and QM/MM calculations
"""
from pathlib import Path
from typing import Dict, Any, List, Optional
import subprocess
import re
from .base_wrapper import SimulationWrapper


class CP2KWrapper(SimulationWrapper):
    """Wrapper for CP2K DFT/QM-MM package"""

    def _detect_app_path(self) -> Path:
        """Detect CP2K installation path"""
        # Check in APPS directory first
        apps_path = Path(__file__).parent.parent / "APPS" / "cp2k"
        if apps_path.exists():
            return apps_path

        # Check if cp2k.ssmp is in PATH
        try:
            result = subprocess.run(['which', 'cp2k.ssmp'], capture_output=True, text=True)
            if result.returncode == 0:
                return Path(result.stdout.strip()).parent.parent
        except:
            pass

        # Check for cp2k.popt (parallel version)
        try:
            result = subprocess.run(['which', 'cp2k.popt'], capture_output=True, text=True)
            if result.returncode == 0:
                return Path(result.stdout.strip()).parent.parent
        except:
            pass

        return apps_path

    def _validate_installation(self) -> bool:
        """Validate that CP2K executable exists"""
        cp2k_exe = self._find_executable()
        if cp2k_exe and cp2k_exe.exists():
            try:
                result = subprocess.run([str(cp2k_exe), '--version'], capture_output=True, text=True, timeout=5)
                if result.returncode == 0:
                    version_match = re.search(r'CP2K version (\d+\.\d+(?:\.\d+)?)', result.stdout)
                    if version_match:
                        self.version = version_match.group(1)
                    else:
                        self.version = "unknown"
                    return True
            except:
                pass

        self.version = "not_installed"
        return False

    def _find_executable(self) -> Optional[Path]:
        """Find CP2K executable (try different variants)"""
        # Try different CP2K executable names
        for exe_name in ['cp2k.ssmp', 'cp2k.popt', 'cp2k']:
            # Check in exe directory
            exe_path = self.app_path / "exe" / "local" / exe_name
            if exe_path.exists():
                return exe_path

            # Check in bin directory
            bin_path = self.app_path / "bin" / exe_name
            if bin_path.exists():
                return bin_path

            # Check system PATH
            try:
                result = subprocess.run(['which', exe_name], capture_output=True, text=True)
                if result.returncode == 0:
                    return Path(result.stdout.strip())
            except:
                pass

        return None

    def _generate_input_file(self, job_params: Dict[str, Any], input_file: Path) -> None:
        """
        Generate CP2K input file

        Expected job_params:
        - run_type: 'ENERGY', 'GEO_OPT', 'MD', 'CELL_OPT' (default: 'ENERGY')
        - project_name: Name for the calculation (default: 'CP2K_calc')
        - atomic_structure: Dict with 'atoms' and 'cell' or 'molecule_smiles'
        - cutoff: Plane wave cutoff in Ry (default: 400)
        - method: 'DFT', 'QMMM' (default: 'DFT')
        - functional: 'PBE', 'BLYP', etc. (default: 'PBE')
        """
        run_type = job_params.get('run_type', 'ENERGY')
        project_name = job_params.get('project_name', 'CP2K_calc')
        cutoff = job_params.get('cutoff', 400.0)
        functional = job_params.get('functional', 'PBE')

        # Get atomic structure (priority: atomic_structure > xyz_structure > molecule_smiles > default)
        if 'atomic_structure' in job_params:
            structure = job_params['atomic_structure']
        elif 'cif_structure' in job_params:
            # Parse CIF format string
            structure = self._parse_cif(job_params['cif_structure'])
        elif 'xyz_structure' in job_params:
            structure = self._parse_xyz(job_params['xyz_structure'])
        elif 'molecule_smiles' in job_params:
            structure = self._smiles_to_structure(job_params['molecule_smiles'])
        else:
            structure = self._default_structure()

        atoms = structure['atoms']
        cell = structure.get('cell', [[20.0, 0, 0], [0, 20.0, 0], [0, 0, 20.0]])

        # Build input file
        input_content = f"""&GLOBAL
  PROJECT {project_name}
  RUN_TYPE {run_type}
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep

  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS

    &MGRID
      CUTOFF {cutoff}
      REL_CUTOFF 60
    &END MGRID

    &QS
      EPS_DEFAULT 1.0E-10
    &END QS

    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 50
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF

    &XC
      &XC_FUNCTIONAL {functional}
      &END XC_FUNCTIONAL
    &END XC

    &PRINT
      &MO
        EIGENVECTORS .FALSE.
      &END MO
    &END PRINT
  &END DFT

  &SUBSYS
    &CELL
      A {cell[0][0]:.6f} {cell[0][1]:.6f} {cell[0][2]:.6f}
      B {cell[1][0]:.6f} {cell[1][1]:.6f} {cell[1][2]:.6f}
      C {cell[2][0]:.6f} {cell[2][1]:.6f} {cell[2][2]:.6f}
      PERIODIC NONE
    &END CELL

    &COORD
"""

        # Add atomic positions
        for atom in atoms:
            input_content += f"      {atom['element']:3s}  {atom['x']:12.8f}  {atom['y']:12.8f}  {atom['z']:12.8f}\n"

        input_content += """    &END COORD

"""

        # Add basis sets and potentials for each element
        unique_elements = sorted(set(atom['element'] for atom in atoms))
        for element in unique_elements:
            input_content += f"""    &KIND {element}
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-{functional}
    &END KIND

"""

        input_content += """  &END SUBSYS
&END FORCE_EVAL
"""

        # Write input file
        with open(input_file, 'w') as f:
            f.write(input_content)

    def _parse_xyz(self, xyz_string: str) -> Dict[str, Any]:
        """
        Parse XYZ format string into atomic structure dict

        XYZ format:
        <number of atoms>
        <comment line>
        <element> <x> <y> <z>
        <element> <x> <y> <z>
        ...
        """
        lines = xyz_string.strip().split('\n')

        if len(lines) < 3:
            raise ValueError("Invalid XYZ format: need at least 3 lines (n_atoms, comment, coordinates)")

        try:
            n_atoms = int(lines[0].strip())
        except ValueError:
            raise ValueError(f"Invalid XYZ format: first line must be number of atoms, got '{lines[0]}'")

        # Comment line is lines[1], skip it

        atoms = []
        for i, line in enumerate(lines[2:], start=2):
            if not line.strip():
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            try:
                element = parts[0]
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])

                atoms.append({
                    'element': element,
                    'x': x,
                    'y': y,
                    'z': z
                })
            except (ValueError, IndexError) as e:
                raise ValueError(f"Invalid XYZ format at line {i+1}: {line.strip()}")

        if len(atoms) != n_atoms:
            print(f"Warning: XYZ file claims {n_atoms} atoms but found {len(atoms)} coordinate lines")

        if not atoms:
            raise ValueError("No atoms found in XYZ structure")

        print(f"Successfully parsed XYZ structure with {len(atoms)} atoms")
        return {'atoms': atoms}

    def _parse_cif(self, cif_string: str) -> Dict[str, Any]:
        """
        Parse CIF format string into atomic structure dict.
        Uses the CIF parser module.
        """
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from cif_parser import create_cif_parser

        parser = create_cif_parser()
        structure = parser.parse_cif(cif_string)

        # Return in the format expected by the wrapper
        return {
            'atoms': structure['atoms'],
            'cell': structure['cell_matrix']
        }


    def _smiles_to_structure(self, smiles: str) -> Dict[str, Any]:
        """Convert SMILES to atomic structure (simplified)"""
        test_structures = {
            'CCO': {  # Ethanol
                'atoms': [
                    {'element': 'C', 'x': 0.0, 'y': 0.0, 'z': 0.0},
                    {'element': 'C', 'x': 1.5, 'y': 0.0, 'z': 0.0},
                    {'element': 'O', 'x': 2.0, 'y': 1.3, 'z': 0.0},
                    {'element': 'H', 'x': -0.5, 'y': -0.9, 'z': 0.0},
                    {'element': 'H', 'x': -0.5, 'y': 0.5, 'z': 0.9},
                    {'element': 'H', 'x': -0.5, 'y': 0.5, 'z': -0.9},
                    {'element': 'H', 'x': 2.0, 'y': -0.5, 'z': 0.9},
                    {'element': 'H', 'x': 2.0, 'y': -0.5, 'z': -0.9},
                    {'element': 'H', 'x': 2.9, 'y': 1.3, 'z': 0.0},
                ]
            },
            'O': {  # Water
                'atoms': [
                    {'element': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.0},
                    {'element': 'H', 'x': 0.76, 'y': 0.59, 'z': 0.0},
                    {'element': 'H', 'x': -0.76, 'y': 0.59, 'z': 0.0},
                ]
            }
        }
        return test_structures.get(smiles, self._default_structure())

    def _default_structure(self) -> Dict[str, Any]:
        """Default structure: water molecule"""
        return {
            'atoms': [
                {'element': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.0},
                {'element': 'H', 'x': 0.76, 'y': 0.59, 'z': 0.0},
                {'element': 'H', 'x': -0.76, 'y': 0.59, 'z': 0.0},
            ]
        }

    def _parse_output_file(self, output_file: Path, job_dir: Path) -> Dict[str, Any]:
        """
        Parse CP2K output file

        Extracts:
        - Final energy
        - SCF convergence
        - Forces (if available)
        """
        if not output_file.exists():
            return {"error": "Output file not found"}

        with open(output_file, 'r') as f:
            content = f.read()

        results = {}

        # Extract final energy
        energy_pattern = r'ENERGY\| Total FORCE_EVAL.*?:\s+([-\d.]+)'
        energy_matches = re.findall(energy_pattern, content)
        if energy_matches:
            energy_au = float(energy_matches[-1])
            results['energy_hartree'] = energy_au
            results['energy_ev'] = energy_au * 27.2114  # Convert to eV

        # Check SCF convergence
        if 'SCF run converged' in content:
            results['scf_converged'] = True
            scf_pattern = r'outer SCF iter =\s+(\d+)'
            scf_match = re.search(scf_pattern, content)
            if scf_match:
                results['scf_iterations'] = int(scf_match.group(1))
        else:
            results['scf_converged'] = False

        # Extract forces (for optimization runs)
        forces_pattern = r'ATOMIC FORCES in \[a\.u\.\].*?\n(.*?)\n\s+SUM OF ATOMIC FORCES'
        forces_match = re.search(forces_pattern, content, re.DOTALL)
        if forces_match:
            forces = []
            for line in forces_match.group(1).strip().split('\n'):
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        forces.append({
                            'atom': int(parts[0]),
                            'element': parts[1],
                            'fx': float(parts[3]),
                            'fy': float(parts[4]),
                            'fz': float(parts[5])
                        })
                    except ValueError:
                        continue
            if forces:
                results['forces'] = forces

        # Extract timing information
        timing_pattern = r'CP2K\s+\|\s+Total Program Execution\s+([\d.]+)'
        timing_match = re.search(timing_pattern, content)
        if timing_match:
            results['total_time_seconds'] = float(timing_match.group(1))

        return results

    def _get_run_command(self, input_file: Path, output_file: Path, job_dir: Path) -> List[str]:
        """Get command to run CP2K"""
        cp2k_exe = self._find_executable()
        if not cp2k_exe:
            raise FileNotFoundError("CP2K executable not found")

        return [str(cp2k_exe), '-i', input_file.name, '-o', output_file.name]

    def _get_input_filename(self) -> str:
        return "cp2k.inp"

    def _get_output_filename(self) -> str:
        return "cp2k.out"
