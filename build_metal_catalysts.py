#!/usr/bin/env python3
"""
Build XYZ structures for metal-based NH3 synthesis catalysts
"""
import numpy as np

def build_fe3o4_cluster():
    """
    Build a small Fe3O4 (magnetite) cluster fragment
    Represents active site on iron oxide surface
    """
    # Simple tetrahedral Fe3O4 motif
    # 3 Fe atoms + 4 O atoms
    atoms = [
        # Iron atoms (roughly tetrahedral)
        ("Fe", 0.0, 0.0, 0.0),
        ("Fe", 2.5, 0.0, 0.0),
        ("Fe", 1.25, 2.165, 0.0),
        # Oxygen atoms (bridging)
        ("O", 1.25, 0.0, 1.5),
        ("O", 1.25, 0.0, -1.5),
        ("O", 0.625, 1.083, 0.75),
        ("O", 1.875, 1.083, 0.75),
    ]
    return format_xyz(atoms, "Fe3O4 cluster - magnetite active site fragment")


def build_fek_alox_promoted():
    """
    Build Fe-K-AlOx promoted site
    K promotes electron density for N2 activation
    """
    # Fe2O3 patch with K promoter and Al2O3
    atoms = [
        # Iron oxide core
        ("Fe", 0.0, 0.0, 0.0),
        ("Fe", 2.8, 0.0, 0.0),
        ("O", 1.4, 0.0, 1.2),
        ("O", 1.4, 0.0, -1.2),
        ("O", 0.7, 1.5, 0.0),
        # Potassium promoter
        ("K", 1.4, 3.0, 0.0),
        # Aluminum oxide support
        ("Al", -1.5, -1.5, 0.0),
        ("Al", 4.3, -1.5, 0.0),
        ("O", -1.5, -1.5, 1.8),
        ("O", 4.3, -1.5, 1.8),
    ]
    return format_xyz(atoms, "Fe-K-AlOx promoted site")


def build_ru10_cluster():
    """
    Build Ru10 cluster (surface-like decagon)
    Represents Ru(0001) terrace fragment
    """
    # Create a 10-atom Ru cluster in roughly planar arrangement
    # Mimics surface terrace
    atoms = []

    # Central Ru
    atoms.append(("Ru", 0.0, 0.0, 0.0))

    # First ring (6 atoms, hexagonal)
    r1 = 2.7  # Ru-Ru distance ~2.7 Å
    for i in range(6):
        angle = i * np.pi / 3
        x = r1 * np.cos(angle)
        y = r1 * np.sin(angle)
        atoms.append(("Ru", x, y, 0.0))

    # Second layer (3 atoms above)
    for i in range(3):
        angle = i * 2 * np.pi / 3 + np.pi/6
        x = r1 * 0.7 * np.cos(angle)
        y = r1 * 0.7 * np.sin(angle)
        atoms.append(("Ru", x, y, 2.2))

    return format_xyz(atoms, "Ru10 cluster - Ru(0001) terrace fragment")


def build_ru_ba_oxide():
    """
    Build Ru-Ba/oxide promoted surface
    BaO enhances N2 activation
    """
    # Ru8 cluster with BaO promoter
    atoms = []

    # Ru cluster core (8 atoms in two layers)
    # Bottom layer (4 atoms, square)
    spacing = 2.7
    for i in range(2):
        for j in range(2):
            x = i * spacing - spacing/2
            y = j * spacing - spacing/2
            atoms.append(("Ru", x, y, 0.0))

    # Top layer (4 atoms, offset square)
    for i in range(2):
        for j in range(2):
            x = i * spacing
            y = j * spacing
            atoms.append(("Ru", x, y, 2.2))

    # BaO promoter on top
    atoms.append(("Ba", spacing/2, spacing/2, 4.5))
    atoms.append(("O", spacing/2, spacing/2, 6.0))

    return format_xyz(atoms, "Ru-Ba/oxide promoted surface")


def format_xyz(atoms, description):
    """Format atoms as XYZ string"""
    lines = [str(len(atoms)), description]
    for element, x, y, z in atoms:
        lines.append(f"{element:2s}  {x:12.8f}  {y:12.8f}  {z:12.8f}")
    return "\n".join(lines)


def save_all_structures():
    """Build and save all 4 structures"""
    structures = {
        "fe3o4_cluster": build_fe3o4_cluster(),
        "fek_alox_promoted": build_fek_alox_promoted(),
        "ru10_cluster": build_ru10_cluster(),
        "ru_ba_oxide": build_ru_ba_oxide()
    }

    for name, xyz in structures.items():
        filename = f"{name}.xyz"
        with open(filename, 'w') as f:
            f.write(xyz)
        print(f"✓ Created {filename}")

    return structures


if __name__ == "__main__":
    print("Building metal catalyst structures for NH3 synthesis...")
    print("="*80)

    structures = save_all_structures()

    print("\n" + "="*80)
    print("Structure Summary:")
    print("="*80)

    for name, xyz in structures.items():
        lines = xyz.split('\n')
        n_atoms = lines[0]
        description = lines[1]
        print(f"\n{name}:")
        print(f"  Atoms: {n_atoms}")
        print(f"  Description: {description}")

    print("\n✓ All structures created!")
    print("\nNext: Screen with MACE-MP model for rapid energy estimates")
