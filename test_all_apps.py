#!/usr/bin/env python3
"""
Comprehensive test of all 5 computational chemistry applications with caffeine molecule
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time
import json

CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


def test_app(agent, app_name, request):
    """Test a single application"""
    print(f"\n{'='*80}", flush=True)
    print(f"Testing {app_name}", flush=True)
    print(f"{'='*80}", flush=True)

    print(f"\nRequest: {request[:100]}...", flush=True)

    start = time.time()
    try:
        result = agent.run_agentic_workflow(request)
        elapsed = time.time() - start

        print(f"\n{'-'*80}", flush=True)
        print(f"RESULT for {app_name}:", flush=True)
        print(f"{'-'*80}", flush=True)

        if result and result.get('status') == 'SUCCESS':
            res = result.get('result', {})
            print(f"✓ Success!", flush=True)
            print(f"  Time: {elapsed:.1f} seconds", flush=True)
            print(f"  Application: {res.get('application', 'N/A')}", flush=True)
            print(f"  Experiment: {res.get('experiment_name', 'N/A')}", flush=True)

            results_data = res.get('results', {})
            if results_data:
                print(f"\n  Key Results:", flush=True)
                for key, value in list(results_data.items())[:10]:  # First 10 keys
                    if not isinstance(value, (list, dict)):
                        print(f"    {key}: {value}", flush=True)

            return True, elapsed, res
        else:
            print(f"✗ Failed", flush=True)
            print(f"  Status: {result.get('status', 'UNKNOWN')}", flush=True)
            if result and 'result' in result and 'error' in result['result']:
                print(f"  Error: {result['result']['error']}", flush=True)
            return False, elapsed, None

    except Exception as e:
        elapsed = time.time() - start
        print(f"✗ Exception: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return False, elapsed, None


def main():
    """Test all 5 applications"""
    print(f"{'='*80}", flush=True)
    print(f"COMPREHENSIVE TEST: All 5 Applications with Caffeine", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"\nMolecule: Caffeine (C8H10N4O2)", flush=True)
    print(f"SMILES: {CAFFEINE_SMILES}", flush=True)
    print(f"Model: gpt-oss:120b", flush=True)

    # Initialize agent
    print(f"\nInitializing agent...", flush=True)
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )
    print(f"✓ Agent ready", flush=True)

    # Define tests for each application
    tests = [
        ("Quantum ESPRESSO", f"""Use Quantum ESPRESSO to calculate the total electronic energy
of caffeine molecule (SMILES: {CAFFEINE_SMILES}).
Perform an SCF calculation with standard parameters."""),

        ("CP2K", f"""Use CP2K to calculate the total energy of caffeine molecule
(SMILES: {CAFFEINE_SMILES}).
Run an ENERGY calculation using PBE functional."""),

        ("GPAW", f"""Run a GPAW calculation to determine the total energy of caffeine
(SMILES: {CAFFEINE_SMILES}).
Use finite-difference mode with PBE functional."""),

        ("LAMMPS", f"""Use LAMMPS to calculate the potential energy of caffeine molecule
(SMILES: {CAFFEINE_SMILES}).
Perform energy minimization using a classical force field."""),

        ("GROMACS", f"""Run a GROMACS energy minimization on caffeine molecule
(SMILES: {CAFFEINE_SMILES}).
Use steepest descents to find the minimum energy structure.""")
    ]

    # Run tests
    results = {}
    total_start = time.time()

    for app_name, request in tests:
        success, elapsed, res = test_app(agent, app_name, request)
        results[app_name] = {
            'success': success,
            'time': elapsed,
            'result': res
        }

    total_time = time.time() - total_start

    # Summary
    print(f"\n{'='*80}", flush=True)
    print(f"SUMMARY", flush=True)
    print(f"{'='*80}", flush=True)

    success_count = sum(1 for r in results.values() if r['success'])

    print(f"\nResults:", flush=True)
    for app_name, data in results.items():
        status = "✓ SUCCESS" if data['success'] else "✗ FAILED"
        time_str = f"{data['time']:.1f}s" if data['time'] else "N/A"
        print(f"  {app_name:20s}: {status:15s}  Time: {time_str}", flush=True)

    print(f"\nOverall: {success_count}/{len(tests)} applications succeeded", flush=True)
    print(f"Total time: {total_time:.1f} seconds", flush=True)

    if success_count == len(tests):
        print(f"\n🎉 ALL APPLICATIONS WORKING!", flush=True)
        return 0
    else:
        print(f"\n⚠ {len(tests) - success_count} application(s) need attention", flush=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
