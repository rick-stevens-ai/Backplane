#!/usr/bin/env python3
"""
Extract a single FeMo-cofactor cluster from PDB
"""
import sys

def extract_single_femo(pdb_file, output_xyz, chain='A', resid='6496'):
    """Extract CFN residue (FeMo-cofactor) from specific chain"""

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    atoms = []

    for line in lines:
        if not line.startswith('HETATM'):
            continue

        # Parse PDB line
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21:22]
        res_id = line[22:26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        element = line[76:78].strip() if len(line) > 76 else atom_name[0:2]

        # Extract CFN residue from specified chain
        if res_name == 'CFN' and chain_id == chain and res_id == resid:
            element = element.strip()
            if not element or element == '':
                # Try to guess from atom name
                if atom_name.startswith('FE'):
                    element = 'Fe'
                elif atom_name.startswith('MO'):
                    element = 'Mo'
                elif atom_name.startswith('S'):
                    element = 'S'
                elif atom_name.startswith('C'):
                    element = 'C'
                elif atom_name.startswith('O'):
                    element = 'O'
                elif atom_name.startswith('H'):
                    element = 'H'

            atoms.append((element, x, y, z, atom_name))

    if not atoms:
        print(f"ERROR: No CFN atoms found in chain {chain} residue {resid}", file=sys.stderr)
        return False

    # Count elements
    elem_count = {}
    for elem, _, _, _, _ in atoms:
        elem = elem.title() if len(elem) == 2 else elem.upper()
        elem_count[elem] = elem_count.get(elem, 0) + 1

    print(f"Extracted {len(atoms)} atoms from CFN {chain}:{resid}", file=sys.stderr)
    print(f"Composition: {elem_count}", file=sys.stderr)

    # Write XYZ file
    with open(output_xyz, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"FeMo-cofactor from PDB {pdb_file} chain {chain} residue {resid}\n")
        for elem, x, y, z, name in atoms:
            elem = elem.title() if len(elem) == 2 else elem.upper()
            f.write(f"{elem:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")

    print(f"Wrote to {output_xyz}", file=sys.stderr)
    return True


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.pdb> <output.xyz> [chain] [resid]")
        sys.exit(1)

    pdb_file = sys.argv[1]
    xyz_file = sys.argv[2]
    chain = sys.argv[3] if len(sys.argv) > 3 else 'A'
    resid = sys.argv[4] if len(sys.argv) > 4 else '6496'

    success = extract_single_femo(pdb_file, xyz_file, chain, resid)
    sys.exit(0 if success else 1)
