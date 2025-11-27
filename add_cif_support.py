#!/usr/bin/env python3
"""
Automated script to add CIF support to all computational chemistry wrappers
Adds cif_structure parsing to CP2K, GPAW, LAMMPS, and GROMACS wrappers
"""

import re
from pathlib import Path

# Define the wrappers to update
WRAPPERS = [
    'wrappers/cp2k_wrapper.py',
    'wrappers/gpaw_wrapper.py',
    'wrappers/lammps_wrapper.py',
    'wrappers/gromacs_wrapper.py'
]

# Pattern to match the structure input section
OLD_PATTERN = r"(# Get atomic structure.*?)\n(\s+if 'atomic_structure' in job_params:.*?\n\s+structure = job_params\['atomic_structure'\]\n)(\s+elif 'xyz_structure' in job_params:)"

# Replacement pattern with CIF support added
NEW_PATTERN = r"\1\n\2\3elif 'cif_structure' in job_params:\n            # Parse CIF format string\n            structure = self._parse_cif(job_params['cif_structure'])\n        \3"

# CIF parser method to add
CIF_PARSE_METHOD = '''    def _parse_cif(self, cif_string: str) -> Dict[str, Any]:
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

'''

def add_cif_to_wrapper(wrapper_path: str):
    """Add CIF support to a wrapper file"""
    path = Path(wrapper_path)

    if not path.exists():
        print(f"❌ {wrapper_path}: File not found")
        return False

    # Read file content
    with open(path, 'r') as f:
        content = f.read()

    # Check if already has CIF support
    if 'cif_structure' in content:
        print(f"✓ {wrapper_path}: Already has CIF support")
        return True

    # Step 1: Add cif_structure to input chain
    pattern1 = re.compile(
        r"(elif 'xyz_structure' in job_params:)",
        re.MULTILINE
    )

    replacement1 = r"elif 'cif_structure' in job_params:\n            # Parse CIF format string\n            structure = self._parse_cif(job_params['cif_structure'])\n        \1"

    if pattern1.search(content):
        content = pattern1.sub(replacement1, content)
        print(f"  → Added cif_structure to input chain")
    else:
        print(f"⚠ {wrapper_path}: Could not find xyz_structure pattern")
        return False

    # Step 2: Add _parse_cif method after _parse_xyz method
    pattern2 = re.compile(
        r"(def _parse_xyz\(self.*?\n(?:.*?\n)*?        return \{[^}]+\})\n",
        re.MULTILINE | re.DOTALL
    )

    if pattern2.search(content):
        # Find the end of _parse_xyz method
        match = pattern2.search(content)
        if match:
            insert_pos = match.end()
            content = content[:insert_pos] + "\n" + CIF_PARSE_METHOD + content[insert_pos:]
            print(f"  → Added _parse_cif() method")
    else:
        print(f"⚠ {wrapper_path}: Could not find _parse_xyz method")
        return False

    # Write updated content
    with open(path, 'w') as f:
        f.write(content)

    print(f"✅ {wrapper_path}: CIF support added successfully\n")
    return True

def main():
    print("="*80)
    print("Adding CIF Support to Computational Chemistry Wrappers")
    print("="*80)
    print()

    success_count = 0
    for wrapper in WRAPPERS:
        if add_cif_to_wrapper(wrapper):
            success_count += 1

    print("="*80)
    print(f"Summary: {success_count}/{len(WRAPPERS)} wrappers updated successfully")
    print("="*80)

    if success_count == len(WRAPPERS):
        print("\n✅ All wrappers now support CIF format!")
        print("\nUsage:")
        print("  job_params = {'cif_structure': cif_string, ...}")
        print("\nPriority chain: atomic_structure → cif_structure → xyz_structure → smiles → default")
    else:
        print(f"\n⚠ {len(WRAPPERS) - success_count} wrappers need manual update")

if __name__ == "__main__":
    main()
