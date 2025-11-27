#!/usr/bin/env python3
"""
Quantum ESPRESSO (QE) wrapper for plane-wave DFT calculations
"""
from pathlib import Path
from typing import Dict, Any, List, Optional
import subprocess
import re
from .base_wrapper import SimulationWrapper


class QuantumEspressoWrapper(SimulationWrapper):
    """Wrapper for Quantum ESPRESSO DFT package"""

    def _detect_app_path(self) -> Path:
        """Detect Quantum ESPRESSO installation path"""
        # Check in APPS directory first
        apps_path = Path(__file__).parent.parent / "APPS" / "q-e"
        if apps_path.exists():
            return apps_path

        # Check if pw.x is in PATH
        try:
            result = subprocess.run(['which', 'pw.x'], capture_output=True, text=True)
            if result.returncode == 0:
                return Path(result.stdout.strip()).parent.parent
        except:
            pass

        # Return APPS path even if not fully installed yet
        return apps_path

    def _validate_installation(self) -> bool:
        """Validate that pw.x executable exists"""
        pw_x = self._find_executable("pw.x")
        if pw_x and pw_x.exists():
            try:
                result = subprocess.run([str(pw_x), '-v'], capture_output=True, text=True, timeout=5)
                # pw.x returns error code but still prints version, so check stdout
                if 'Program PWSCF' in result.stdout:
                    # Extract version from output
                    version_match = re.search(r'v\.(\d+\.\d+(?:\.\d+)?)', result.stdout)
                    if version_match:
                        self.version = version_match.group(1)
                    else:
                        self.version = "unknown"
                    return True
            except:
                pass

        self.version = "not_installed"
        return False

    def _find_executable(self, name: str) -> Optional[Path]:
        """Find executable in app_path or system PATH"""
        # Check in bin directory
        bin_path = self.app_path / "bin" / name
        if bin_path.exists():
            return bin_path

        # Check system PATH
        try:
            result = subprocess.run(['which', name], capture_output=True, text=True)
            if result.returncode == 0:
                return Path(result.stdout.strip())
        except:
            pass

        return None

    def _generate_input_file(self, job_params: Dict[str, Any], input_file: Path) -> None:
        """
        Generate Quantum ESPRESSO PWscf input file

        Expected job_params:
        - calculation: 'scf', 'relax', 'vc-relax', 'md' (default: 'scf')
        - system_name: Name for the calculation (default: 'QE_calculation')
        - atomic_structure: Dict with 'atoms' and 'cell' or 'molecule_smiles'
        - cutoff_wfc: Wavefunction cutoff in Ry (default: 50)
        - cutoff_rho: Charge density cutoff in Ry (default: 400)
        - k_points: K-point grid [nk1, nk2, nk3] (default: [1,1,1] for molecules)
        - pseudopotentials: Dict mapping element to pseudopotential file
        """
        calculation = job_params.get('calculation', 'scf')
        system_name = job_params.get('system_name', 'QE_calculation')
        cutoff_wfc = job_params.get('cutoff_wfc', 50.0)
        cutoff_rho = job_params.get('cutoff_rho', 400.0)
        k_points = job_params.get('k_points', [1, 1, 1])

        # Get atomic structure - supports atomic_structure dict, CIF, XYZ, or SMILES
        if 'atomic_structure' in job_params:
            structure = job_params['atomic_structure']
        elif 'cif_structure' in job_params:
            # Parse CIF format string
            structure = self._parse_cif(job_params['cif_structure'])
        elif 'xyz_structure' in job_params:
            # Parse XYZ format string
            structure = self._parse_xyz(job_params['xyz_structure'])
        elif 'molecule_smiles' in job_params:
            # For small test: create simple structure from SMILES
            structure = self._smiles_to_structure(job_params['molecule_smiles'])
        else:
            # Default: water molecule for testing
            structure = self._default_structure()

        atoms = structure['atoms']
        cell = structure.get('cell', [[20.0, 0, 0], [0, 20.0, 0], [0, 0, 20.0]])  # Default: 20 Å cube for molecules

        # Count atoms and get unique elements
        nat = len(atoms)
        ntyp = len(set(atom['element'] for atom in atoms))

        # Build input file
        pseudo_dir = str(Path(__file__).parent.parent / "pseudo")
        input_content = f"""&CONTROL
  calculation = '{calculation}'
  prefix = '{system_name}'
  outdir = './tmp'
  pseudo_dir = '{pseudo_dir}'
  tprnfor = .true.
  tstress = .true.
/

&SYSTEM
  ibrav = 0
  nat = {nat}
  ntyp = {ntyp}
  ecutwfc = {cutoff_wfc}
  ecutrho = {cutoff_rho}
/

&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.7
/

"""

        # Add IONS section for relaxation
        if calculation in ['relax', 'vc-relax', 'md']:
            input_content += """&IONS
  ion_dynamics = 'bfgs'
/

"""

        # Add CELL section for variable-cell relaxation
        if calculation == 'vc-relax':
            input_content += """&CELL
  cell_dynamics = 'bfgs'
/

"""

        # Atomic species
        input_content += "ATOMIC_SPECIES\n"
        unique_elements = sorted(set(atom['element'] for atom in atoms))

        # Map elements to pseudopotential files
        pseudo_files = {
            'H': 'H.pbe-kjpaw_psl.1.0.0.UPF',
            'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
            'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
            'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF'
        }

        for element in unique_elements:
            mass = self._get_atomic_mass(element)
            pseudo = pseudo_files.get(element, f"{element}.UPF")
            input_content += f"  {element}  {mass}  {pseudo}\n"

        # Cell parameters
        input_content += "\nCELL_PARAMETERS angstrom\n"
        for row in cell:
            input_content += f"  {row[0]:12.8f}  {row[1]:12.8f}  {row[2]:12.8f}\n"

        # Atomic positions
        input_content += "\nATOMIC_POSITIONS angstrom\n"
        for atom in atoms:
            input_content += f"  {atom['element']:3s}  {atom['x']:12.8f}  {atom['y']:12.8f}  {atom['z']:12.8f}\n"

        # K-points
        input_content += f"\nK_POINTS automatic\n"
        input_content += f"  {k_points[0]} {k_points[1]} {k_points[2]}  0 0 0\n"

        # Write input file
        with open(input_file, 'w') as f:
            f.write(input_content)

    def _smiles_to_structure(self, smiles: str) -> Dict[str, Any]:
        """
        Convert SMILES to atomic structure (simplified for small molecules)
        For production, use RDKit or similar library
        """
        # Simple test structure for common molecules
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

    def _parse_xyz(self, xyz_string: str) -> Dict[str, Any]:
        """
        Parse XYZ format string into atomic structure dict.

        XYZ format:
        <number of atoms>
        <comment line>
        <element> <x> <y> <z>
        ...
        """
        lines = xyz_string.strip().split('\n')

        if len(lines) < 3:
            raise ValueError("Invalid XYZ format: need at least 3 lines (n_atoms, comment, atom lines)")

        try:
            n_atoms = int(lines[0].strip())
        except ValueError:
            raise ValueError(f"Invalid XYZ format: first line must be number of atoms, got: {lines[0]}")

        # Parse atoms starting from line 3 (skip n_atoms and comment)
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
                raise ValueError(f"Invalid XYZ format at line {i+1}: {line.strip()} - {str(e)}")

        if len(atoms) != n_atoms:
            # This is a warning, not an error - some XYZ files may have extra blank lines
            pass

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

    def _get_atomic_mass(self, element: str) -> float:
        """Get atomic mass for element"""
        masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
            'Fe': 55.845, 'Cu': 63.546, 'Pt': 195.084, 'Ru': 101.07,
            'Mo': 95.95, 'Co': 58.933, 'Ni': 58.693
        }
        return masses.get(element, 1.0)

    def _parse_output_file(self, output_file: Path, job_dir: Path) -> Dict[str, Any]:
        """
        Parse Quantum ESPRESSO output file

        Extracts:
        - Final energy
        - Forces
        - Stress
        - Convergence status
        """
        if not output_file.exists():
            return {"error": "Output file not found"}

        with open(output_file, 'r') as f:
            content = f.read()

        results = {}

        # Extract final energy
        energy_pattern = r'!\s+total energy\s+=\s+([-\d.]+)\s+Ry'
        energy_matches = re.findall(energy_pattern, content)
        if energy_matches:
            energy_ry = float(energy_matches[-1])
            results['energy_ry'] = energy_ry
            results['energy_ev'] = energy_ry * 13.6057  # Convert to eV

        # Check convergence
        if 'convergence has been achieved' in content or 'JOB DONE' in content:
            results['converged'] = True
        else:
            results['converged'] = False

        # Extract forces
        forces_section = re.search(r'Forces acting on atoms.*?\n(.*?)\n\s*Total force', content, re.DOTALL)
        if forces_section:
            forces = []
            for line in forces_section.group(1).strip().split('\n'):
                if 'atom' in line:
                    parts = line.split()
                    if len(parts) >= 9:
                        forces.append({
                            'atom': int(parts[1]),
                            'fx': float(parts[6]),
                            'fy': float(parts[7]),
                            'fz': float(parts[8])
                        })
            if forces:
                results['forces'] = forces

        # Extract stress (if present)
        stress_pattern = r'total\s+stress.*?=\s+([-\d.]+)'
        stress_match = re.search(stress_pattern, content)
        if stress_match:
            results['total_stress_kbar'] = float(stress_match.group(1))

        # Extract number of iterations
        scf_pattern = r'convergence has been achieved in\s+(\d+)\s+iterations'
        scf_match = re.search(scf_pattern, content)
        if scf_match:
            results['scf_iterations'] = int(scf_match.group(1))

        return results

    def _get_run_command(self, input_file: Path, output_file: Path, job_dir: Path) -> List[str]:
        """Get command to run pw.x"""
        pw_x = self._find_executable("pw.x")
        if not pw_x:
            raise FileNotFoundError("pw.x executable not found")

        return [str(pw_x), '-in', input_file.name]

    def _get_input_filename(self) -> str:
        return "pw.in"

    def _get_output_filename(self) -> str:
        return "pw.out"
