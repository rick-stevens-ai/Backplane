#!/usr/bin/env python3
"""
CIF (Crystallographic Information File) Parser
Converts CIF files to internal structure representation for computational chemistry wrappers
"""
import re
import math
from typing import Dict, List, Tuple, Optional, Any


class CIFParser:
    """Parser for CIF (Crystallographic Information File) format"""

    def __init__(self):
        """Initialize CIF parser"""
        pass

    def parse_cif(self, cif_string: str) -> Dict[str, Any]:
        """
        Parse CIF file content into internal structure representation

        Args:
            cif_string: Content of CIF file as string

        Returns:
            Dictionary containing:
                - atoms: List of atom dictionaries with element, fractional coords
                - cell: Unit cell parameters [a, b, c, alpha, beta, gamma]
                - cell_matrix: 3x3 cell matrix in Cartesian coordinates
                - space_group: Space group symbol (if available)
                - formula: Chemical formula
        """
        lines = cif_string.strip().split('\n')

        # Extract cell parameters
        cell_params = self._extract_cell_parameters(lines)

        # Extract space group
        space_group = self._extract_space_group(lines)

        # Extract atomic positions
        atoms = self._extract_atoms(lines)

        # Convert cell parameters to matrix
        cell_matrix = self._cell_params_to_matrix(cell_params)

        # Convert fractional to Cartesian coordinates
        atoms_cartesian = self._fractional_to_cartesian(atoms, cell_matrix)

        # Determine chemical formula
        formula = self._determine_formula(atoms)

        return {
            'atoms': atoms_cartesian,
            'atoms_fractional': atoms,
            'cell': cell_params,
            'cell_matrix': cell_matrix,
            'space_group': space_group,
            'formula': formula
        }

    def _extract_cell_parameters(self, lines: List[str]) -> List[float]:
        """Extract cell parameters a, b, c, alpha, beta, gamma from CIF"""
        params = {
            '_cell_length_a': None,
            '_cell_length_b': None,
            '_cell_length_c': None,
            '_cell_angle_alpha': None,
            '_cell_angle_beta': None,
            '_cell_angle_gamma': None
        }

        for line in lines:
            line = line.strip()
            for key in params.keys():
                if line.startswith(key):
                    # Extract numeric value (handle uncertainties in parentheses)
                    value_str = line.split()[1]
                    # Remove uncertainty notation like 8.396(3)
                    value_str = re.sub(r'\([^)]*\)', '', value_str)
                    try:
                        params[key] = float(value_str)
                    except ValueError:
                        pass

        # Check if all parameters found
        if None in params.values():
            raise ValueError("Could not extract all cell parameters from CIF")

        return [
            params['_cell_length_a'],
            params['_cell_length_b'],
            params['_cell_length_c'],
            params['_cell_angle_alpha'],
            params['_cell_angle_beta'],
            params['_cell_angle_gamma']
        ]

    def _extract_space_group(self, lines: List[str]) -> Optional[str]:
        """Extract space group from CIF"""
        for line in lines:
            line = line.strip()
            if line.startswith('_space_group_name_H-M') or \
               line.startswith('_symmetry_space_group_name_H-M'):
                # Extract space group name (may be quoted)
                match = re.search(r"'([^']+)'", line)
                if match:
                    return match.group(1)
                # Try without quotes
                parts = line.split()
                if len(parts) > 1:
                    return ' '.join(parts[1:])
        return None

    def _extract_atoms(self, lines: List[str]) -> List[Dict[str, Any]]:
        """Extract atomic positions from CIF"""
        atoms = []

        # Find atom site loop
        in_loop = False
        loop_columns = []

        for i, line in enumerate(lines):
            line = line.strip()

            # Start of atom site loop
            if line.startswith('loop_'):
                in_loop = True
                loop_columns = []
                continue

            # Column headers in loop
            if in_loop and line.startswith('_atom_site'):
                loop_columns.append(line)
                continue

            # End of loop or start of new section
            if in_loop and (line.startswith('_') and not line.startswith('_atom_site')):
                in_loop = False
                loop_columns = []
                continue

            # Parse atom data
            if in_loop and loop_columns and not line.startswith('_'):
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) >= len(loop_columns):
                    atom = self._parse_atom_line(parts, loop_columns)
                    if atom:
                        atoms.append(atom)

        if not atoms:
            raise ValueError("Could not extract atomic positions from CIF")

        return atoms

    def _parse_atom_line(self, parts: List[str], columns: List[str]) -> Optional[Dict[str, Any]]:
        """Parse a single atom line based on column headers"""
        atom = {}

        # Map column indices
        col_map = {}
        for i, col in enumerate(columns):
            col_map[col] = i

        try:
            # Extract element symbol
            if '_atom_site_type_symbol' in col_map:
                idx = col_map['_atom_site_type_symbol']
                # Remove oxidation state (e.g., Fe3+ -> Fe)
                element = re.sub(r'[0-9+-]', '', parts[idx])
                atom['element'] = element
            elif '_atom_site_label' in col_map:
                idx = col_map['_atom_site_label']
                # Extract element from label (e.g., Fe1 -> Fe)
                element = re.sub(r'[0-9+-]', '', parts[idx])
                atom['element'] = element
            else:
                return None

            # Extract fractional coordinates
            if '_atom_site_fract_x' in col_map:
                x_idx = col_map['_atom_site_fract_x']
                y_idx = col_map['_atom_site_fract_y']
                z_idx = col_map['_atom_site_fract_z']

                # Remove uncertainties
                x_str = re.sub(r'\([^)]*\)', '', parts[x_idx])
                y_str = re.sub(r'\([^)]*\)', '', parts[y_idx])
                z_str = re.sub(r'\([^)]*\)', '', parts[z_idx])

                atom['fract_x'] = float(x_str)
                atom['fract_y'] = float(y_str)
                atom['fract_z'] = float(z_str)
            else:
                return None

            # Optional: occupancy
            if '_atom_site_occupancy' in col_map:
                occ_idx = col_map['_atom_site_occupancy']
                occ_str = re.sub(r'\([^)]*\)', '', parts[occ_idx])
                atom['occupancy'] = float(occ_str)
            else:
                atom['occupancy'] = 1.0

            return atom

        except (ValueError, IndexError) as e:
            return None

    def _cell_params_to_matrix(self, params: List[float]) -> List[List[float]]:
        """
        Convert cell parameters to 3x3 Cartesian cell matrix

        Args:
            params: [a, b, c, alpha, beta, gamma] in Angstroms and degrees

        Returns:
            3x3 matrix where rows are lattice vectors in Cartesian coordinates
        """
        a, b, c, alpha, beta, gamma = params

        # Convert angles to radians
        alpha_rad = math.radians(alpha)
        beta_rad = math.radians(beta)
        gamma_rad = math.radians(gamma)

        # Calculate cell matrix
        # a vector along x-axis
        ax = a
        ay = 0.0
        az = 0.0

        # b vector in xy plane
        bx = b * math.cos(gamma_rad)
        by = b * math.sin(gamma_rad)
        bz = 0.0

        # c vector
        cx = c * math.cos(beta_rad)
        cy = c * (math.cos(alpha_rad) - math.cos(beta_rad) * math.cos(gamma_rad)) / math.sin(gamma_rad)
        cz = c * math.sqrt(1.0 - math.cos(beta_rad)**2 -
                          ((math.cos(alpha_rad) - math.cos(beta_rad) * math.cos(gamma_rad)) /
                           math.sin(gamma_rad))**2)

        return [
            [ax, ay, az],
            [bx, by, bz],
            [cx, cy, cz]
        ]

    def _fractional_to_cartesian(self, atoms: List[Dict], cell_matrix: List[List[float]]) -> List[Dict]:
        """Convert fractional coordinates to Cartesian"""
        cartesian_atoms = []

        for atom in atoms:
            fract = [atom['fract_x'], atom['fract_y'], atom['fract_z']]

            # Matrix multiplication: cart = cell_matrix^T * fract
            x = (cell_matrix[0][0] * fract[0] +
                 cell_matrix[1][0] * fract[1] +
                 cell_matrix[2][0] * fract[2])
            y = (cell_matrix[0][1] * fract[0] +
                 cell_matrix[1][1] * fract[1] +
                 cell_matrix[2][1] * fract[2])
            z = (cell_matrix[0][2] * fract[0] +
                 cell_matrix[1][2] * fract[1] +
                 cell_matrix[2][2] * fract[2])

            cartesian_atoms.append({
                'element': atom['element'],
                'x': x,
                'y': y,
                'z': z,
                'fract_x': atom['fract_x'],
                'fract_y': atom['fract_y'],
                'fract_z': atom['fract_z'],
                'occupancy': atom.get('occupancy', 1.0)
            })

        return cartesian_atoms

    def _determine_formula(self, atoms: List[Dict]) -> str:
        """Determine chemical formula from atoms"""
        element_counts = {}

        for atom in atoms:
            element = atom['element']
            occupancy = atom.get('occupancy', 1.0)

            if element not in element_counts:
                element_counts[element] = 0.0
            element_counts[element] += occupancy

        # Build formula string
        formula_parts = []
        for element in sorted(element_counts.keys()):
            count = element_counts[element]
            if abs(count - round(count)) < 0.01:
                count = int(round(count))
                if count == 1:
                    formula_parts.append(element)
                else:
                    formula_parts.append(f"{element}{count}")
            else:
                formula_parts.append(f"{element}{count:.2f}")

        return ''.join(formula_parts)

    def to_xyz(self, structure: Dict[str, Any]) -> str:
        """
        Convert parsed CIF structure to XYZ format string

        Args:
            structure: Parsed CIF structure from parse_cif()

        Returns:
            XYZ format string with extended format including cell
        """
        atoms = structure['atoms']
        cell_matrix = structure['cell_matrix']

        # Flatten cell matrix for XYZ extended format
        cell_flat = []
        for row in cell_matrix:
            cell_flat.extend(row)

        cell_str = ' '.join([f"{x:.6f}" for x in cell_flat])

        # Build XYZ string
        xyz_lines = [
            str(len(atoms)),
            f'Lattice="{cell_str}" Properties=species:S:1:pos:R:3'
        ]

        for atom in atoms:
            xyz_lines.append(
                f"{atom['element']:3s} {atom['x']:12.6f} {atom['y']:12.6f} {atom['z']:12.6f}"
            )

        return '\n'.join(xyz_lines)


def create_cif_parser() -> CIFParser:
    """Factory function to create CIF parser"""
    return CIFParser()


# Example usage and testing
if __name__ == "__main__":
    # Example: Fe3O4 magnetite CIF
    example_cif = """
data_Fe3O4
_chemical_name_mineral 'Magnetite'
_chemical_formula_sum 'Fe3 O4'
_cell_length_a    8.396
_cell_length_b    8.396
_cell_length_c    8.396
_cell_angle_alpha 90.000
_cell_angle_beta  90.000
_cell_angle_gamma 90.000
_space_group_name_H-M_alt 'F d -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0.125 0.125 0.125 1.0
Fe2 Fe 0.5 0.5 0.5 1.0
O1 O 0.255 0.255 0.255 1.0
"""

    parser = create_cif_parser()

    try:
        structure = parser.parse_cif(example_cif)

        print("Parsed CIF Structure:")
        print(f"Formula: {structure['formula']}")
        print(f"Space Group: {structure['space_group']}")
        print(f"Cell Parameters: {structure['cell']}")
        print(f"\nAtoms ({len(structure['atoms'])}):")
        for atom in structure['atoms']:
            print(f"  {atom['element']:3s} ({atom['x']:8.4f}, {atom['y']:8.4f}, {atom['z']:8.4f})")

        print("\nXYZ Format:")
        print(parser.to_xyz(structure))

    except Exception as e:
        print(f"Error parsing CIF: {e}")
