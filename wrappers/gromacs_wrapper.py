#!/usr/bin/env python3
"""
GROMACS wrapper for molecular dynamics of biochemical systems
"""
from pathlib import Path
from typing import Dict, Any, List, Optional
import subprocess
import re
from .base_wrapper import SimulationWrapper


class GROMACSWrapper(SimulationWrapper):
    """Wrapper for GROMACS molecular dynamics package"""

    def _detect_app_path(self) -> Path:
        """Detect GROMACS installation path"""
        # Check if gmx is in PATH
        try:
            result = subprocess.run(['which', 'gmx'], capture_output=True, text=True)
            if result.returncode == 0:
                gmx_path = Path(result.stdout.strip())
                return gmx_path.parent.parent
        except:
            pass

        # Check APPS directory
        apps_path = Path(__file__).parent.parent / "APPS" / "gromacs"
        return apps_path

    def _validate_installation(self) -> bool:
        """Validate that GROMACS is properly installed"""
        try:
            result = subprocess.run(['gmx', '--version'], capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                version_match = re.search(r'GROMACS version:\s+(\S+)', result.stdout)
                if version_match:
                    self.version = version_match.group(1)
                else:
                    self.version = "unknown"
                return True
        except:
            pass

        self.version = "not_installed"
        return False

    def _generate_input_file(self, job_params: Dict[str, Any], input_file: Path) -> None:
        """
        Generate GROMACS input files (.mdp and structure files)

        Expected job_params:
        - simulation_type: 'em' (energy minimization), 'nvt', 'npt', 'md' (default: 'em')
        - system_name: Name for the calculation (default: 'GROMACS_sim')
        - atomic_structure: Dict with 'atoms' and 'cell' or 'molecule_smiles'
        - temperature: Temperature in K (default: 300)
        - pressure: Pressure in bar (for NPT, default: 1.0)
        - timestep: Timestep in ps (default: 0.002)
        - steps: Number of MD steps (default: 50000)
        - force_field: Force field name (default: 'oplsaa')
        """
        sim_type = job_params.get('simulation_type', 'em')
        system_name = job_params.get('system_name', 'GROMACS_sim')
        temperature = job_params.get('temperature', 300.0)
        pressure = job_params.get('pressure', 1.0)
        timestep = job_params.get('timestep', 0.002)
        steps = job_params.get('steps', 50000)
        force_field = job_params.get('force_field', 'oplsaa')

        # Get atomic structure - supports atomic_structure dict, XYZ, or SMILES
        if 'atomic_structure' in job_params:
            structure = job_params['atomic_structure']
        elif 'cif_structure' in job_params:
            # Parse CIF format string
            structure = self._parse_cif(job_params['cif_structure'])
        elif 'xyz_structure' in job_params:
            # Parse XYZ format string
            structure = self._parse_xyz(job_params['xyz_structure'])
        elif 'molecule_smiles' in job_params:
            structure = self._smiles_to_structure(job_params['molecule_smiles'])
        else:
            structure = self._default_structure()

        atoms = structure['atoms']
        cell = structure.get('cell', [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]])  # Default: 3 nm cube

        # Store structure for topology generation
        self._current_structure = structure

        # Generate structure file (.gro)
        gro_file = input_file.parent / "system.gro"
        self._generate_gro_file(atoms, cell, gro_file, system_name)

        # Generate topology file (.top)
        top_file = input_file.parent / "topol.top"
        self._generate_topology_file(atoms, structure.get('bonds', []), top_file)

        # Generate MDP file based on simulation type
        if sim_type == 'em':
            self._generate_em_mdp(input_file, steps)
        elif sim_type == 'nvt':
            self._generate_nvt_mdp(input_file, temperature, timestep, steps)
        elif sim_type == 'npt':
            self._generate_npt_mdp(input_file, temperature, pressure, timestep, steps)
        else:  # generic MD
            self._generate_md_mdp(input_file, temperature, timestep, steps)

    def _generate_gro_file(self, atoms: List[Dict], cell: List[List[float]],
                          gro_file: Path, system_name: str) -> None:
        """Generate GROMACS structure file (.gro)"""
        nat = len(atoms)

        content = f"{system_name}\n"
        content += f"{nat}\n"

        for i, atom in enumerate(atoms, 1):
            # GROMACS .gro format: %5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f
            # residue_number, residue_name, atom_name, atom_number, x, y, z (in nm)
            x_nm = atom['x'] / 10.0  # Convert Angstrom to nm
            y_nm = atom['y'] / 10.0
            z_nm = atom['z'] / 10.0
            content += f"{1:5d}{'MOL':<5s}{atom['element']:>5s}{i:5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n"

        # Box vectors (in nm)
        box_x = cell[0][0] / 10.0
        box_y = cell[1][1] / 10.0
        box_z = cell[2][2] / 10.0
        content += f"   {box_x:.5f}   {box_y:.5f}   {box_z:.5f}\n"

        with open(gro_file, 'w') as f:
            f.write(content)

    def _generate_em_mdp(self, mdp_file: Path, steps: int) -> None:
        """Generate energy minimization MDP file"""
        content = f"""; Energy Minimization
integrator  = steep
nsteps      = {steps}
emtol       = 100.0
emstep      = 0.01

; Output control
nstlog      = 100
nstenergy   = 100

; Neighbor searching
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 10
rlist       = 1.2

; Electrostatics
coulombtype = PME
rcoulomb    = 1.2

; van der Waals
vdwtype     = Cut-off
rvdw        = 1.2

; Temperature and pressure
tcoupl      = no
pcoupl      = no

; Periodic boundary conditions
pbc         = xyz

; Constraints
constraints = none
"""
        with open(mdp_file, 'w') as f:
            f.write(content)

    def _generate_nvt_mdp(self, mdp_file: Path, temperature: float,
                         timestep: float, steps: int) -> None:
        """Generate NVT (constant volume) MDP file"""
        content = f"""; NVT ensemble
integrator  = md
dt          = {timestep}
nsteps      = {steps}

; Output control
nstlog      = 1000
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000

; Neighbor searching
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 10
rlist       = 1.2

; Electrostatics
coulombtype = PME
rcoulomb    = 1.2

; van der Waals
vdwtype     = Cut-off
rvdw        = 1.2

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = {temperature}

; Pressure coupling
pcoupl      = no

; Velocity generation
gen_vel     = yes
gen_temp    = {temperature}
gen_seed    = -1

; Periodic boundary conditions
pbc         = xyz

; Constraints
constraints = h-bonds
"""
        with open(mdp_file, 'w') as f:
            f.write(content)

    def _generate_npt_mdp(self, mdp_file: Path, temperature: float,
                         pressure: float, timestep: float, steps: int) -> None:
        """Generate NPT (constant pressure) MDP file"""
        content = f"""; NPT ensemble
integrator  = md
dt          = {timestep}
nsteps      = {steps}

; Output control
nstlog      = 1000
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000

; Neighbor searching
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 10
rlist       = 1.2

; Electrostatics
coulombtype = PME
rcoulomb    = 1.2

; van der Waals
vdwtype     = Cut-off
rvdw        = 1.2

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = {temperature}

; Pressure coupling
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = {pressure}
compressibility = 4.5e-5

; Velocity generation
gen_vel     = yes
gen_temp    = {temperature}
gen_seed    = -1

; Periodic boundary conditions
pbc         = xyz

; Constraints
constraints = h-bonds
"""
        with open(mdp_file, 'w') as f:
            f.write(content)

    def _generate_md_mdp(self, mdp_file: Path, temperature: float,
                        timestep: float, steps: int) -> None:
        """Generate generic MD MDP file"""
        self._generate_nvt_mdp(mdp_file, temperature, timestep, steps)

    def _generate_topology_file(self, atoms: List[Dict], bonds: List[Dict], top_file: Path) -> None:
        """Generate minimal GROMACS topology file"""
        nat = len(atoms)
        unique_elements = sorted(set(atom['element'] for atom in atoms))

        # Create atom types
        atom_types = []
        for elem in unique_elements:
            mass = self._get_atomic_mass(elem)
            # Minimal LJ parameters (epsilon in kJ/mol, sigma in nm)
            lj_params = {
                'C': (0.439, 0.3431),
                'H': (0.126, 0.2500),
                'N': (0.289, 0.3660),
                'O': (0.586, 0.3118),
            }
            epsilon, sigma = lj_params.get(elem, (0.4, 0.3))
            atom_types.append((elem, mass, 0.0, epsilon, sigma))

        # Build topology content
        content = """; Topology file generated by Backplane
; Minimal topology for molecule simulation

[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
  1             2               yes             0.5     0.8333

[ atomtypes ]
; name  mass      charge   ptype   sigma      epsilon
"""
        for elem, mass, charge, epsilon, sigma in atom_types:
            content += f"{elem:6s}  {mass:8.4f}  {charge:8.4f}  A  {sigma:10.5f}  {epsilon:10.5f}\n"

        content += """
[ moleculetype ]
; Name            nrexcl
MOL              3

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
"""
        elem_to_type = {elem: elem for elem in unique_elements}
        for i, atom in enumerate(atoms, 1):
            elem = atom['element']
            mass = self._get_atomic_mass(elem)
            content += f"{i:5d}  {elem:4s}  1  MOL  {elem:4s}  1  0.000  {mass:8.4f}\n"

        # Add bonds if present
        if bonds:
            content += """
[ bonds ]
; ai  aj  funct  c0         c1
"""
            for bond in bonds:
                # GROMACS bond: type 1 (harmonic), length in nm, force constant in kJ/mol/nm^2
                content += f"{bond['atom1']+1:5d}  {bond['atom2']+1:5d}  1  0.15  250000.0\n"

        content += """
[ system ]
; Name
Molecule in water

[ molecules ]
; Compound  nmols
MOL         1
"""

        with open(top_file, 'w') as f:
            f.write(content)

    def _get_atomic_mass(self, element: str) -> float:
        """Get atomic mass for element"""
        masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
        }
        return masses.get(element, 12.0)

    def _smiles_to_structure(self, smiles: str) -> Dict[str, Any]:
        """Convert SMILES to atomic structure using RDKit"""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return self._default_structure()

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

            # Get conformer
            conf = mol.GetConformer()

            # Extract atoms and positions
            atoms = []
            bonds = []

            # Center molecule at (15, 15, 15)
            offset_x, offset_y, offset_z = 15.0, 15.0, 15.0

            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                atoms.append({
                    'element': atom.GetSymbol(),
                    'x': pos.x + offset_x,
                    'y': pos.y + offset_y,
                    'z': pos.z + offset_z,
                })

            # Extract bonds
            for bond in mol.GetBonds():
                bonds.append({
                    'atom1': bond.GetBeginAtomIdx(),
                    'atom2': bond.GetEndAtomIdx(),
                    'type': int(bond.GetBondTypeAsDouble())
                })

            # Determine box size (molecule extent + 10 Å padding)
            coords = [(a['x'], a['y'], a['z']) for a in atoms]
            max_x = max(c[0] for c in coords) + 5.0
            max_y = max(c[1] for c in coords) + 5.0
            max_z = max(c[2] for c in coords) + 5.0
            min_x = min(c[0] for c in coords) - 5.0
            min_y = min(c[1] for c in coords) - 5.0
            min_z = min(c[2] for c in coords) - 5.0

            box_size = max(max_x - min_x, max_y - min_y, max_z - min_z, 30.0)

            return {
                'atoms': atoms,
                'bonds': bonds,
                'cell': [[box_size, 0, 0], [0, box_size, 0], [0, 0, box_size]]
            }

        except ImportError:
            # RDKit not available, use default
            return self._default_structure()
        except Exception as e:
            print(f"Error converting SMILES: {e}")
            return self._default_structure()

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
            raise ValueError("Invalid XYZ format: need at least 3 lines")

        try:
            n_atoms = int(lines[0].strip())
        except ValueError:
            raise ValueError(f"Invalid XYZ format: first line must be number of atoms")

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

        # Determine box size based on molecule extent (convert to nm for GROMACS)
        if atoms:
            coords = [(a['x'], a['y'], a['z']) for a in atoms]
            max_x = max(c[0] for c in coords) + 10.0  # Add 10 Angstrom padding
            max_y = max(c[1] for c in coords) + 10.0
            max_z = max(c[2] for c in coords) + 10.0
            min_x = min(c[0] for c in coords) - 10.0
            min_y = min(c[1] for c in coords) - 10.0
            min_z = min(c[2] for c in coords) - 10.0

            # Convert to Angstrom, ensure minimum 20 Angstrom
            box_x = max(max_x - min_x, 20.0)
            box_y = max(max_y - min_y, 20.0)
            box_z = max(max_z - min_z, 20.0)

            return {
                'atoms': atoms,
                'cell': [[box_x, 0, 0], [0, box_y, 0], [0, 0, box_z]]
            }

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


    def _default_structure(self) -> Dict[str, Any]:
        """Default structure: water molecule"""
        return {
            'atoms': [
                {'element': 'O', 'x': 15.0, 'y': 15.0, 'z': 15.0},
                {'element': 'H', 'x': 15.76, 'y': 15.59, 'z': 15.0},
                {'element': 'H', 'x': 14.24, 'y': 15.59, 'z': 15.0},
            ],
            'cell': [[30.0, 0, 0], [0, 30.0, 0], [0, 0, 30.0]]
        }

    def _parse_output_file(self, output_file: Path, job_dir: Path) -> Dict[str, Any]:
        """
        Parse GROMACS output

        Extracts:
        - Final energy
        - Temperature, pressure (if applicable)
        - Convergence status
        """
        if not output_file.exists():
            return {"error": "Output file not found"}

        with open(output_file, 'r') as f:
            content = f.read()

        results = {}

        # Extract final potential energy
        energy_pattern = r'Potential Energy\s+=\s+([-\d.eE+]+)'
        energy_matches = re.findall(energy_pattern, content)
        if energy_matches:
            energy_kj = float(energy_matches[-1])
            results['potential_energy_kj_mol'] = energy_kj
            results['potential_energy_kcal_mol'] = energy_kj / 4.184

        # Extract temperature
        temp_pattern = r'Temperature\s+=\s+([-\d.]+)'
        temp_matches = re.findall(temp_pattern, content)
        if temp_matches:
            results['temperature_K'] = float(temp_matches[-1])

        # Extract pressure
        press_pattern = r'Pressure\s+=\s+([-\d.eE+]+)'
        press_matches = re.findall(press_pattern, content)
        if press_matches:
            results['pressure_bar'] = float(press_matches[-1])

        # Check for successful completion
        if 'Finished mdrun' in content:
            results['completed'] = True
        else:
            results['completed'] = False

        # Extract timing
        time_pattern = r'Time:\s+([\d.]+)'
        time_match = re.search(time_pattern, content)
        if time_match:
            results['simulation_time_ps'] = float(time_match.group(1))

        # For energy minimization, check convergence
        if 'Steepest Descents converged' in content:
            results['converged'] = True
            results['minimization_type'] = 'steepest_descents'
        elif 'Energy minimization' in content and 'reached the maximum number' in content:
            results['converged'] = False
            results['minimization_type'] = 'steepest_descents'

        return results

    def _get_run_command(self, input_file: Path, output_file: Path, job_dir: Path) -> List[str]:
        """
        Get command to run GROMACS
        This is a simplified version - real GROMACS workflows require grompp first
        """
        # In practice, GROMACS requires:
        # 1. gmx grompp -f input.mdp -c system.gro -o topol.tpr
        # 2. gmx mdrun -s topol.tpr -o output.trr -e output.edr -g output.log

        # For this simplified version, we'll create a shell script
        script_file = job_dir / "run_gromacs.sh"
        gro_file = job_dir / "system.gro"

        script_content = f"""#!/bin/bash
# GROMACS run script

# Preprocessing (grompp)
gmx grompp -f {input_file.name} -c system.gro -p topol.top -o topol.tpr -maxwarn 10 2>&1

# Run simulation
gmx mdrun -s topol.tpr -o traj.trr -e ener.edr -g {output_file.name} -c confout.gro 2>&1
"""

        with open(script_file, 'w') as f:
            f.write(script_content)
        script_file.chmod(0o755)

        return ['bash', script_file.name]

    def _get_input_filename(self) -> str:
        return "gromacs.mdp"

    def _get_output_filename(self) -> str:
        return "md.log"
