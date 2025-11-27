#!/usr/bin/env python3
"""
Minimal validation test for computational chemistry applications
Quick check to verify all wrappers can initialize and find executables
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from wrappers.quantum_espresso import QuantumEspressoWrapper
from wrappers.cp2k_wrapper import CP2KWrapper
from wrappers.gpaw_wrapper import GPAWWrapper
from wrappers.lammps_wrapper import LAMMPSWrapper
from wrappers.gromacs_wrapper import GROMACSWrapper


def test_wrapper(wrapper_class, app_name):
    """Test if wrapper can initialize"""
    print(f"\nTesting {app_name}...")
    print("─" * 70)

    try:
        wrapper = wrapper_class()
        print(f"✓ Wrapper initialized")
        print(f"  App path: {wrapper.app_path}")
        print(f"  Version: {wrapper.version}")

        if wrapper.version == "not_installed":
            print(f"  ⚠ Warning: Executable not found or not working")
            return False

        # Try to prepare a test job (doesn't execute)
        print(f"  Testing input generation...")
        job_dir = wrapper.prepare_job(
            job_id="test_validation",
            job_params={
                "molecule_smiles": "O",  # Simple water molecule
                "system_name": "validation_test"
            }
        )
        print(f"✓ Input generation successful")
        print(f"  Test directory: {job_dir}")

        # Check if input file was created
        input_file = job_dir / wrapper._get_input_filename()
        if input_file.exists():
            print(f"✓ Input file created: {input_file.name}")
            size = input_file.stat().st_size
            print(f"  Size: {size} bytes")

        return True

    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run minimal validation tests"""
    print("="*70)
    print("MINIMAL APPLICATION VALIDATION TEST")
    print("="*70)
    print("\nThis test checks if wrappers can initialize and generate input files")
    print("without actually running calculations.\n")

    results = {}

    # Test each application
    results['Quantum ESPRESSO'] = test_wrapper(QuantumEspressoWrapper, "Quantum ESPRESSO")
    results['CP2K'] = test_wrapper(CP2KWrapper, "CP2K")
    results['GPAW'] = test_wrapper(GPAWWrapper, "GPAW")
    results['LAMMPS'] = test_wrapper(LAMMPSWrapper, "LAMMPS")
    results['GROMACS'] = test_wrapper(GROMACSWrapper, "GROMACS")

    # Summary
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)

    success_count = sum(1 for r in results.values() if r)
    total = len(results)

    for app, success in results.items():
        status = "✓ READY" if success else "✗ NOT READY"
        print(f"{app:20s}: {status}")

    print(f"\n{success_count}/{total} applications ready")

    if success_count == total:
        print("\n✓ ALL APPLICATIONS READY!")
        print("You can now run: python test_caffeine_energy.py")
        return 0
    else:
        print(f"\n⚠ {total - success_count} application(s) need attention")
        return 1


if __name__ == "__main__":
    sys.exit(main())
