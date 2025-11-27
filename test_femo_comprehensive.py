#!/usr/bin/env python3
"""
Comprehensive test of FeMo-cofactor across multiple computational chemistry tools
Tests with a simplified Fe-S cluster representation
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time
import json

def test_tool(agent, tool_name, request):
    """Test a single computational tool"""
    print(f"\n{'='*80}", flush=True)
    print(f"TESTING: {tool_name}", flush=True)
    print(f"{'='*80}\n", flush=True)

    start = time.time()

    try:
        result = agent.run_agentic_workflow(request)
        elapsed = time.time() - start

        if result and result.get('status') == 'SUCCESS':
            res = result.get('result', {})
            print(f"\n✓ {tool_name} SUCCESS", flush=True)
            print(f"  Time: {elapsed:.1f}s", flush=True)
            print(f"  Application: {res.get('application', 'N/A')}", flush=True)

            results_data = res.get('results', {})
            if results_data:
                print(f"  Key Results:", flush=True)
                for key, value in list(results_data.items())[:10]:
                    if not isinstance(value, (list, dict)):
                        print(f"    {key}: {value}", flush=True)

            return {
                'success': True,
                'time': elapsed,
                'application': res.get('application'),
                'results': results_data
            }
        else:
            print(f"\n✗ {tool_name} did not complete", flush=True)
            print(f"  Time: {elapsed:.1f}s", flush=True)
            return {
                'success': False,
                'time': elapsed,
                'status': result.get('status') if result else 'NO_RESULT'
            }

    except Exception as e:
        elapsed = time.time() - start
        print(f"\n✗ {tool_name} EXCEPTION: {e}", flush=True)
        return {
            'success': False,
            'time': elapsed,
            'error': str(e)
        }


def main():
    """Run comprehensive FeMo-cofactor tests"""

    print(f"{'='*80}", flush=True)
    print(f"COMPREHENSIVE FeMo-COFACTOR SIMULATION TEST", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"\nReal structure: [Fe7MoS9N] cluster from nitrogenase", flush=True)
    print(f"Simplified model: Fe-S cluster for computational feasibility", flush=True)
    print(f"Model: gpt-oss:120b\n", flush=True)

    # Initialize agent
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    # Test suite - using simplified iron-sulfur clusters
    # Full FeMo-cofactor is too large for simple SMILES, so we test with smaller models
    tests = {
        "CP2K (DFT - Fe4S4 model)": """Use CP2K to calculate the energy of an Fe4S4 cubane cluster
as a simplified model of the FeMo-cofactor active site. This is a common iron-sulfur cluster motif.
Use PBE functional with spin-polarized DFT. The cluster should be treated as a neutral species.
Focus on the electronic structure of this transition metal cluster.""",

        "GPAW (DFT - Fe2S2 model)": """Use GPAW to study a simplified Fe2S2 cluster as a model for
understanding iron-sulfur chemistry relevant to the FeMo-cofactor. Use PBE functional with
finite-difference mode. Calculate the total energy and electronic properties.""",

        "CP2K (Geometry optimization)": """Use CP2K to perform geometry optimization of a small
iron-sulfur cluster (Fe2S2) as a model system. Use PBE functional and optimize the structure to find
the minimum energy configuration. This will help understand metal-sulfur bonding.""",

        "GPAW (Electronic structure)": """Use GPAW to analyze the electronic structure of a metal-sulfide
cluster. Calculate the density of states and orbital energies for a simple Fe-S system. This provides
insight into the electronic properties of metallo-clusters like FeMo-cofactor.""",

        "Classical MD (Fe-S network)": """Use GROMACS or LAMMPS to perform molecular dynamics on a simplified
representation of metal-sulfur interactions. Model the structural dynamics and energy landscape of a
sulfur-containing metal complex using classical force fields."""
    }

    results = {}
    total_start = time.time()

    for test_name, request in tests.items():
        result = test_tool(agent, test_name, request)
        results[test_name] = result
        time.sleep(2)  # Brief pause between tests

    total_time = time.time() - total_start

    # Summary
    print(f"\n\n{'='*80}", flush=True)
    print(f"FINAL SUMMARY", flush=True)
    print(f"{'='*80}\n", flush=True)

    success_count = sum(1 for r in results.values() if r.get('success', False))

    for test_name, data in results.items():
        status = "✓ SUCCESS" if data.get('success') else "✗ FAILED"
        time_str = f"{data.get('time', 0):.1f}s"
        print(f"{test_name:40s}: {status:15s}  Time: {time_str}", flush=True)

    print(f"\nTotal: {success_count}/{len(tests)} tests succeeded", flush=True)
    print(f"Total time: {total_time/60:.1f} minutes", flush=True)

    # Save detailed results
    with open('femo_results.json', 'w') as f:
        json.dump(results, indent=2, fp=f)
    print(f"\nDetailed results saved to: femo_results.json", flush=True)

    return success_count == len(tests)


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
