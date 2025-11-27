#!/usr/bin/env python3
"""
Extract FeMo-cofactor from PDB structure
Extracts the [Fe7MoS9C] cluster with homocitrate ligand
"""
import sys

def extract_femo_cofactor(pdb_file, output_xyz):
    """
    Extract FeMo-cofactor from PDB file
    Looking for residue name 'FMO' or atoms Fe, Mo, S near each other
    """

    # Read PDB file
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    # Find HETATM lines with metal cluster components
    # FeMo-cofactor is usually labeled as 'FES', 'MOS', or specific residue names
    cluster_atoms = []

    for line in lines:
        if not line.startswith('HETATM'):
            continue

        # Parse PDB line
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        element = line[76:78].strip() if len(line) > 76 else atom_name[0]

        # Look for FeMo cluster components
        # Common residue names: FMO, FES, MOS, CLF (cluster with Fe)
        if res_name in ['FMO', 'FES', 'MOS', 'CLF', '4FE', 'SF4']:
            cluster_atoms.append((element, x, y, z, res_name, atom_name))
        # Also catch individual Fe, Mo, S atoms that might be labeled differently
        elif element in ['FE', 'MO', 'S'] and res_name in ['FE2', 'MO', 'SF4', 'CLF']:
            cluster_atoms.append((element, x, y, z, res_name, atom_name))

    if not cluster_atoms:
        print("No FeMo-cofactor atoms found. Searching for any Fe/Mo/S clusters...", file=sys.stderr)

        # Broader search - any Fe, Mo, S
        for line in lines:
            if not line.startswith('HETATM'):
                continue

            element = line[76:78].strip() if len(line) > 76 else line[12:16].strip()[0]
            res_name = line[17:20].strip()

            if element in ['FE', 'MO', 'S', 'C'] and res_name != 'HOH':
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_name = line[12:16].strip()
                cluster_atoms.append((element, x, y, z, res_name, atom_name))

    if not cluster_atoms:
        print("ERROR: Could not find FeMo-cofactor atoms!", file=sys.stderr)
        return False

    print(f"Found {len(cluster_atoms)} potential cluster atoms", file=sys.stderr)

    # Group by proximity - FeMo cluster should be within ~10 Angstroms
    # Find the central cluster by looking for high concentration of Fe/Mo/S

    # Count elements
    element_counts = {}
    for atom in cluster_atoms:
        elem = atom[0]
        element_counts[elem] = element_counts.get(elem, 0) + 1

    print(f"Element counts: {element_counts}", file=sys.stderr)

    # Filter to get reasonable FeMo-cofactor size
    # Should have ~7 Fe, 1 Mo, 9 S, plus C and homocitrate (C, O, H)
    femo_atoms = []
    for atom in cluster_atoms:
        elem, x, y, z, res, name = atom
        # Include metal and sulfur atoms
        if elem in ['FE', 'MO', 'S']:
            femo_atoms.append(atom)
        # Include carbon if it might be central carbide
        elif elem == 'C' and res in ['FMO', 'FES', 'CLF', 'C']:
            femo_atoms.append(atom)

    # If we have atoms, write XYZ file
    if femo_atoms:
        with open(output_xyz, 'w') as f:
            f.write(f"{len(femo_atoms)}\n")
            f.write(f"FeMo-cofactor extracted from {pdb_file}\n")
            for elem, x, y, z, res, name in femo_atoms:
                # Clean up element symbol
                elem = elem.title() if len(elem) == 2 else elem.upper()
                f.write(f"{elem:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")

        print(f"Wrote {len(femo_atoms)} atoms to {output_xyz}", file=sys.stderr)
        print(f"Composition: {dict((e, sum(1 for a in femo_atoms if a[0] == e)) for e in set(a[0] for a in femo_atoms))}", file=sys.stderr)
        return True

    return False


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.pdb> <output.xyz>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    xyz_file = sys.argv[2]

    success = extract_femo_cofactor(pdb_file, xyz_file)
    sys.exit(0 if success else 1)
