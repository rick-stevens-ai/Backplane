#!/usr/bin/env python3
"""
Test script for computational chemistry applications
Tests all five applications with gpt-oss:120b model
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent


def test_quantum_espresso(agent: ComputationalChemistryAgent):
    """Test Quantum ESPRESSO DFT calculation"""
    print("\n" + "="*80)
    print("TEST 1: QUANTUM ESPRESSO - Electronic Structure Calculation")
    print("="*80)

    user_request = """Please use Quantum ESPRESSO to perform a DFT geometry optimization
    on a water molecule (O) to analyze its electronic structure. This will help understand
    its properties for catalyst design."""

    try:
        result = agent.run_agentic_workflow(user_request)
        print("\n" + "="*80)
        print("✓ QUANTUM ESPRESSO TEST COMPLETED")
        print("="*80)
        return result
    except Exception as e:
        print(f"\n✗ QUANTUM ESPRESSO TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_cp2k(agent: ComputationalChemistryAgent):
    """Test CP2K mixed DFT calculation"""
    print("\n" + "="*80)
    print("TEST 2: CP2K - Mixed Gaussian/Plane-wave DFT")
    print("="*80)

    user_request = """Use CP2K to calculate the electronic structure and energy
    of an ethanol molecule (CCO) using the PBE functional. This is important for
    understanding alcohol-based catalysis."""

    try:
        result = agent.run_agentic_workflow(user_request)
        print("\n" + "="*80)
        print("✓ CP2K TEST COMPLETED")
        print("="*80)
        return result
    except Exception as e:
        print(f"\n✗ CP2K TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_gpaw(agent: ComputationalChemistryAgent):
    """Test GPAW real-space DFT calculation"""
    print("\n" + "="*80)
    print("TEST 3: GPAW - Real-space Grid DFT in Python")
    print("="*80)

    user_request = """Run a GPAW calculation to determine the electronic properties
    of a water molecule (O) using real-space grid methods. GPAW's Python integration
    makes it ideal for rapid prototyping."""

    try:
        result = agent.run_agentic_workflow(user_request)
        print("\n" + "="*80)
        print("✓ GPAW TEST COMPLETED")
        print("="*80)
        return result
    except Exception as e:
        print(f"\n✗ GPAW TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_lammps(agent: ComputationalChemistryAgent):
    """Test LAMMPS molecular dynamics"""
    print("\n" + "="*80)
    print("TEST 4: LAMMPS - Classical Molecular Dynamics")
    print("="*80)

    user_request = """Use LAMMPS to perform an energy minimization on an ethanol
    molecule (CCO) using classical molecular dynamics. This will help understand
    the molecular structure for catalyst screening."""

    try:
        result = agent.run_agentic_workflow(user_request)
        print("\n" + "="*80)
        print("✓ LAMMPS TEST COMPLETED")
        print("="*80)
        return result
    except Exception as e:
        print(f"\n✗ LAMMPS TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_gromacs(agent: ComputationalChemistryAgent):
    """Test GROMACS molecular dynamics"""
    print("\n" + "="*80)
    print("TEST 5: GROMACS - Biomolecular Molecular Dynamics")
    print("="*80)

    user_request = """Run a GROMACS energy minimization simulation on a water
    molecule (O) to understand its behavior in enzyme active sites for biochemical
    catalysis studies."""

    try:
        result = agent.run_agentic_workflow(user_request)
        print("\n" + "="*80)
        print("✓ GROMACS TEST COMPLETED")
        print("="*80)
        return result
    except Exception as e:
        print(f"\n✗ GROMACS TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Run all application tests with gpt-oss:120b"""
    print("\n" + "="*80)
    print("COMPUTATIONAL CHEMISTRY APPLICATIONS TEST SUITE")
    print("Model: gpt-oss:120b (120B parameters)")
    print("Server: spark-container-03")
    print("="*80)

    # Initialize agent with gpt-oss:120b
    try:
        agent = ComputationalChemistryAgent(
            server_config_path="spark_servers.yaml",
            server_name="spark-container-03"
        )
        print("✓ Agent initialized successfully")
        print(f"  Model: {agent.model}")
        print(f"  API Base: {agent.server_config['openai_api_base']}")
    except Exception as e:
        print(f"✗ Failed to initialize agent: {e}")
        return 1

    # Track results
    results = {}

    # Test each application
    print("\n" + "="*80)
    print("RUNNING APPLICATION TESTS")
    print("="*80)

    results['quantum_espresso'] = test_quantum_espresso(agent)
    results['cp2k'] = test_cp2k(agent)
    results['gpaw'] = test_gpaw(agent)
    results['lammps'] = test_lammps(agent)
    results['gromacs'] = test_gromacs(agent)

    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)

    success_count = sum(1 for r in results.values() if r is not None)
    total_count = len(results)

    for app, result in results.items():
        status = "✓ PASSED" if result is not None else "✗ FAILED"
        print(f"{app:20s}: {status}")

    print(f"\n{success_count}/{total_count} tests passed")

    if success_count == total_count:
        print("\n✓ ALL TESTS PASSED - Catalyst design applications ready!")
        return 0
    else:
        print(f"\n✗ {total_count - success_count} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
