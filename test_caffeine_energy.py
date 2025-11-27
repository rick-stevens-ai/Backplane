#!/usr/bin/env python3
"""
Caffeine Molecule Energy Calculations Across Five Applications

Tests each computational chemistry application by calculating the energy
of caffeine molecule (C8H10N4O2) - a common stimulant and drug molecule.

Caffeine SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C

This test demonstrates:
- Real molecular system (24 atoms)
- Multiple atom types (C, H, N, O)
- Energy calculations with different methods
- Cross-application comparison
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import json
import time


CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


def print_header(title):
    """Print formatted section header"""
    print("\n" + "="*80)
    print(f"{title:^80}")
    print("="*80)


def print_results(app_name, result):
    """Print formatted results"""
    print(f"\n{'─'*80}")
    print(f"Results for {app_name}:")
    print(f"{'─'*80}")

    if result and result.get('status') == 'SUCCESS':
        res = result.get('result', {})

        # Print application info
        print(f"✓ Calculation completed successfully")
        print(f"  Application: {res.get('application', 'N/A')}")
        print(f"  Experiment: {res.get('experiment_name', 'N/A')}")
        print(f"  Execution time: {res.get('execution_time', 'N/A')} seconds")

        # Print energy results
        print(f"\nEnergy Results:")
        results_data = res.get('results', {})

        # Different applications return energy in different formats
        if 'energy_hartree' in results_data:
            print(f"  Energy: {results_data['energy_hartree']:.6f} Hartree")
            print(f"        = {results_data.get('energy_ev', 'N/A')} eV")
        elif 'energy_ev' in results_data:
            print(f"  Energy: {results_data['energy_ev']:.6f} eV")
        elif 'energy_ry' in results_data:
            print(f"  Energy: {results_data['energy_ry']:.6f} Ry")
            if 'energy_ev' in results_data:
                print(f"        = {results_data['energy_ev']:.6f} eV")
        elif 'potential_energy_kj_mol' in results_data:
            print(f"  Energy: {results_data['potential_energy_kj_mol']:.6f} kJ/mol")
            if 'potential_energy_kcal_mol' in results_data:
                print(f"        = {results_data['potential_energy_kcal_mol']:.6f} kcal/mol")
        elif 'energy_kcal_mol' in results_data:
            print(f"  Energy: {results_data['energy_kcal_mol']:.6f} kcal/mol")
            if 'energy_ev' in results_data:
                print(f"        = {results_data['energy_ev']:.6f} eV")

        # Print convergence info
        if 'converged' in results_data:
            print(f"  Converged: {results_data['converged']}")
        if 'scf_converged' in results_data:
            print(f"  SCF Converged: {results_data['scf_converged']}")
        if 'scf_iterations' in results_data:
            print(f"  SCF Iterations: {results_data['scf_iterations']}")

        # Print additional properties
        if 'forces' in results_data and results_data['forces']:
            print(f"  Forces computed: Yes ({len(results_data['forces'])} atoms)")

        if 'dipole_moment' in results_data:
            print(f"  Dipole moment: {results_data['dipole_moment']}")

        # Print output files
        if 'output_files' in res:
            print(f"\nOutput Files:")
            for key, path in res['output_files'].items():
                print(f"  {key}: {path}")

        print(f"{'─'*80}\n")
        return True

    else:
        print(f"✗ Calculation failed or did not complete")
        if result:
            print(f"  Status: {result.get('status', 'UNKNOWN')}")
            if 'result' in result and 'error' in result['result']:
                print(f"  Error: {result['result']['error']}")
        print(f"{'─'*80}\n")
        return False


def test_quantum_espresso(agent):
    """Test Quantum ESPRESSO: DFT energy calculation"""
    print_header("TEST 1: QUANTUM ESPRESSO - DFT Energy Calculation")

    print("""
Quantum ESPRESSO is a plane-wave DFT code excellent for accurate electronic
structure calculations. We'll perform an SCF (self-consistent field) calculation
to determine the total energy of caffeine.

Method: Plane-wave DFT
Functional: (default, typically PBE)
Cutoff: 50 Ry (wavefunction), 400 Ry (charge density)
""")

    request = f"""Please use Quantum ESPRESSO to calculate the total electronic energy
of caffeine molecule (SMILES: {CAFFEINE_SMILES}).

Perform an SCF (self-consistent field) calculation to get the ground state energy.
Use standard cutoffs for accuracy. This will give us the DFT energy of the system."""

    try:
        print("Submitting job to gpt-oss:120b agent...")
        start = time.time()
        result = agent.run_agentic_workflow(request)
        elapsed = time.time() - start

        success = print_results("Quantum ESPRESSO", result)
        print(f"Total time: {elapsed:.1f} seconds")
        return result if success else None

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_cp2k(agent):
    """Test CP2K: Mixed Gaussian/plane-wave DFT"""
    print_header("TEST 2: CP2K - Mixed Gaussian/Plane-wave DFT")

    print("""
CP2K uses a mixed Gaussian/plane-wave approach, making it efficient for
molecular systems. We'll calculate the total energy using the PBE functional.

Method: Mixed Gaussian/plane-wave DFT
Functional: PBE
Basis: DZVP-MOLOPT-SR-GTH
""")

    request = f"""Use CP2K to calculate the total energy of caffeine molecule
(SMILES: {CAFFEINE_SMILES}).

Run an ENERGY calculation (single point) using the PBE functional.
Use the DZVP-MOLOPT basis set for accurate energetics."""

    try:
        print("Submitting job to gpt-oss:120b agent...")
        start = time.time()
        result = agent.run_agentic_workflow(request)
        elapsed = time.time() - start

        success = print_results("CP2K", result)
        print(f"Total time: {elapsed:.1f} seconds")
        return result if success else None

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_gpaw(agent):
    """Test GPAW: Real-space grid DFT"""
    print_header("TEST 3: GPAW - Real-space Grid DFT")

    print("""
GPAW uses real-space grids and is implemented in Python, making it excellent
for prototyping. We'll calculate the energy using finite-difference method.

Method: Real-space finite-difference DFT
Functional: PBE
Grid spacing: 0.2 Angstrom
""")

    request = f"""Run a GPAW calculation to determine the total energy of caffeine
(SMILES: {CAFFEINE_SMILES}).

Use the finite-difference mode with PBE functional and a grid spacing of 0.2 Angstrom.
Perform a single-point energy calculation."""

    try:
        print("Submitting job to gpt-oss:120b agent...")
        start = time.time()
        result = agent.run_agentic_workflow(request)
        elapsed = time.time() - start

        success = print_results("GPAW", result)
        print(f"Total time: {elapsed:.1f} seconds")
        return result if success else None

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_lammps(agent):
    """Test LAMMPS: Classical force field energy minimization"""
    print_header("TEST 4: LAMMPS - Classical Force Field Energy")

    print("""
LAMMPS uses classical force fields for rapid energy calculations on larger
systems. We'll perform energy minimization to find the optimal structure
and its associated potential energy.

Method: Classical molecular mechanics
Force field: Lennard-Jones (for quick test)
Temperature: 300 K
""")

    request = f"""Use LAMMPS to calculate the minimized potential energy of caffeine
molecule (SMILES: {CAFFEINE_SMILES}).

Perform an energy minimization using a classical force field at 300 K.
This will give us the classical potential energy at the optimized geometry."""

    try:
        print("Submitting job to gpt-oss:120b agent...")
        start = time.time()
        result = agent.run_agentic_workflow(request)
        elapsed = time.time() - start

        success = print_results("LAMMPS", result)
        print(f"Total time: {elapsed:.1f} seconds")
        return result if success else None

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_gromacs(agent):
    """Test GROMACS: Biomolecular force field energy minimization"""
    print_header("TEST 5: GROMACS - Biomolecular Force Field Energy")

    print("""
GROMACS is optimized for biomolecular systems. We'll perform energy minimization
using steepest descents to find the minimum energy structure.

Method: Classical molecular mechanics
Force field: OPLS-AA
Algorithm: Steepest descents
""")

    request = f"""Run a GROMACS energy minimization on caffeine molecule
(SMILES: {CAFFEINE_SMILES}).

Use energy minimization (EM) with the OPLS-AA force field to find the
minimum energy structure. Report the final potential energy."""

    try:
        print("Submitting job to gpt-oss:120b agent...")
        start = time.time()
        result = agent.run_agentic_workflow(request)
        elapsed = time.time() - start

        success = print_results("GROMACS", result)
        print(f"Total time: {elapsed:.1f} seconds")
        return result if success else None

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None


def compare_results(results):
    """Compare energies across all methods"""
    print_header("COMPARATIVE ANALYSIS")

    print("""
Caffeine Molecule: C8H10N4O2 (24 atoms)
SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C

Energy Comparison Across Methods:
""")

    # Collect energies
    energy_table = []

    for app_name, result in results.items():
        if result and result.get('status') == 'SUCCESS':
            res = result.get('result', {})
            results_data = res.get('results', {})

            # Extract energy in eV if possible
            energy_ev = None
            energy_str = "N/A"
            method = ""

            if 'energy_ev' in results_data:
                energy_ev = results_data['energy_ev']
                energy_str = f"{energy_ev:.2f} eV"
                method = "DFT" if app_name in ['Quantum ESPRESSO', 'CP2K', 'GPAW'] else "Classical"
            elif 'energy_hartree' in results_data:
                energy_ev = results_data['energy_hartree'] * 27.2114
                energy_str = f"{results_data['energy_hartree']:.6f} Ha ({energy_ev:.2f} eV)"
                method = "DFT"
            elif 'potential_energy_kcal_mol' in results_data:
                energy_kcal = results_data['potential_energy_kcal_mol']
                energy_str = f"{energy_kcal:.2f} kcal/mol"
                method = "Classical"
            elif 'energy_kcal_mol' in results_data:
                energy_kcal = results_data['energy_kcal_mol']
                energy_str = f"{energy_kcal:.2f} kcal/mol"
                method = "Classical"

            converged = results_data.get('converged', results_data.get('scf_converged', 'N/A'))
            exec_time = res.get('execution_time', 'N/A')

            energy_table.append({
                'application': app_name,
                'method': method,
                'energy': energy_str,
                'converged': converged,
                'time': exec_time
            })

    # Print table
    print(f"{'Application':<20} {'Method':<10} {'Energy':<30} {'Converged':<12} {'Time (s)':<10}")
    print("─" * 90)

    for row in energy_table:
        print(f"{row['application']:<20} {row['method']:<10} {row['energy']:<30} "
              f"{str(row['converged']):<12} {str(row['time']):<10}")

    print("\n" + "─" * 90)

    # Analysis
    print("\nKey Observations:")
    print("─" * 90)

    dft_apps = [r for r in energy_table if r['method'] == 'DFT']
    classical_apps = [r for r in energy_table if r['method'] == 'Classical']

    if dft_apps:
        print("\n1. DFT Calculations (Quantum Mechanical):")
        print("   - Provide electronic structure energies")
        print("   - Include electron-electron interactions")
        print("   - More computationally expensive but more accurate")
        print("   - Energies are typically large negative numbers")

    if classical_apps:
        print("\n2. Classical Force Field Calculations:")
        print("   - Use empirical potentials")
        print("   - Much faster but less accurate for electronic properties")
        print("   - Good for conformational analysis and large systems")
        print("   - Energies represent potential energy of the configuration")

    print("\n3. Method Comparison:")
    print("   - DFT methods (QE, CP2K, GPAW) should give similar energies")
    print("   - Classical methods (LAMMPS, GROMACS) will differ in scale")
    print("   - Energy differences are more meaningful than absolute values")

    success_rate = len(energy_table) / 5 * 100
    print(f"\n✓ Successfully completed {len(energy_table)}/5 calculations ({success_rate:.0f}%)")


def main():
    """Run all caffeine energy calculations"""
    print_header("CAFFEINE ENERGY CALCULATIONS TEST SUITE")

    print("""
This test suite calculates the energy of caffeine molecule (C8H10N4O2)
using five different computational chemistry applications.

Caffeine: A stimulant drug and common biomolecule
SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Atoms: 24 (8 C, 10 H, 4 N, 2 O)

Applications tested:
1. Quantum ESPRESSO - Plane-wave DFT
2. CP2K - Mixed Gaussian/plane-wave DFT
3. GPAW - Real-space grid DFT
4. LAMMPS - Classical MD
5. GROMACS - Biomolecular MD

Model: gpt-oss:120b (120B parameters)
Server: spark-container-03
""")

    # Initialize agent
    try:
        print("\nInitializing agent...")
        agent = ComputationalChemistryAgent(
            server_config_path="spark_servers.yaml",
            server_name="spark-container-03"
        )
        print(f"✓ Agent initialized")
        print(f"  Model: {agent.model}")
        print(f"  API: {agent.server_config['openai_api_base']}")
    except Exception as e:
        print(f"\n✗ Failed to initialize agent: {e}")
        return 1

    # Run all tests
    results = {}

    print("\n" + "="*80)
    print("STARTING CALCULATIONS")
    print("="*80)
    print("\nNote: Each calculation may take several minutes.")
    print("The agent will autonomously set appropriate parameters.\n")

    results['Quantum ESPRESSO'] = test_quantum_espresso(agent)
    results['CP2K'] = test_cp2k(agent)
    results['GPAW'] = test_gpaw(agent)
    results['LAMMPS'] = test_lammps(agent)
    results['GROMACS'] = test_gromacs(agent)

    # Compare results
    compare_results(results)

    # Summary
    print_header("TEST SUMMARY")

    success_count = sum(1 for r in results.values() if r is not None)

    print(f"\nCompleted {success_count}/5 calculations successfully")

    for app, result in results.items():
        status = "✓ PASSED" if result is not None else "✗ FAILED"
        print(f"  {app:<20}: {status}")

    if success_count == 5:
        print("\n✓ ALL CALCULATIONS COMPLETED!")
        print("All five computational chemistry applications are working correctly.")
        return 0
    elif success_count > 0:
        print(f"\n⚠ {success_count}/5 calculations completed")
        print("Some applications may need configuration or installation.")
        return 1
    else:
        print("\n✗ NO CALCULATIONS COMPLETED")
        print("Please check application installations and configurations.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
