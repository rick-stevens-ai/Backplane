#!/usr/bin/env python3
"""
Test the fixed SMILES-to-3D conversion in GPAW wrapper
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from wrappers.gpaw_wrapper import GPAWWrapper

def test_smiles_conversion():
    """Test SMILES conversion for NH3 catalyst molecules"""

    wrapper = GPAWWrapper()

    test_molecules = [
        ('Pyridine', 'c1ccncc1', 11),      # C5H5N
        ('Imidazole', 'c1cnc[nH]1', 9),    # C3H4N2
        ('2,2\'-Bipyridine', 'c1ccnc(c1)c2ccccn2', 20),  # C10H8N2
        ('Ammonia', 'N', 4),                # NH3
        ('Hydrazine', 'NN', 6),             # N2H4
    ]

    print("=" * 80)
    print("TESTING FIXED SMILES-TO-3D CONVERSION")
    print("=" * 80)
    print()

    all_passed = True

    for name, smiles, expected_atoms in test_molecules:
        print(f"Testing {name} (SMILES: {smiles})")
        print(f"  Expected atoms: {expected_atoms}")

        try:
            structure = wrapper._smiles_to_structure(smiles)
            actual_atoms = len(structure['atoms'])
            print(f"  Actual atoms:   {actual_atoms}")

            if actual_atoms == expected_atoms:
                print(f"  ✓ PASS")
            else:
                print(f"  ✗ FAIL - Expected {expected_atoms}, got {actual_atoms}")
                all_passed = False

            # Show first few atoms
            print(f"  First 3 atoms:")
            for i, atom in enumerate(structure['atoms'][:3]):
                print(f"    {atom['element']}: ({atom['x']:.3f}, {atom['y']:.3f}, {atom['z']:.3f})")

        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            all_passed = False

        print()

    print("=" * 80)
    if all_passed:
        print("✓ ALL TESTS PASSED - SMILES conversion working correctly!")
    else:
        print("✗ SOME TESTS FAILED - Check errors above")
    print("=" * 80)

    return all_passed


if __name__ == "__main__":
    success = test_smiles_conversion()
    sys.exit(0 if success else 1)
