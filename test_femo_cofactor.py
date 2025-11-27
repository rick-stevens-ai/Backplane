#!/usr/bin/env python3
"""
Test FeMo-cofactor simulation across multiple applications
The FeMo-cofactor is the iron-molybdenum cofactor [Fe7MoS9C] in nitrogenase
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time
import json

# FeMo-cofactor is too complex for SMILES, so we'll use a simplified model
# or let the agent determine how to represent it

def test_femo_cofactor():
    """Test FeMo-cofactor with available applications"""

    print(f"{'='*80}", flush=True)
    print(f"FeMo-COFACTOR SIMULATION TEST", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"\nMolecule: FeMo-cofactor (Iron-Molybdenum cofactor)", flush=True)
    print(f"Composition: [Fe7MoS9C] with homocitrate ligand", flush=True)
    print(f"Function: Active site of nitrogenase enzyme (N2 fixation)", flush=True)
    print(f"Model: gpt-oss:120b", flush=True)
    print(f"\n{'='*80}\n", flush=True)

    # Initialize agent
    print(f"Initializing agent...", flush=True)
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )
    print(f"✓ Agent ready\n", flush=True)

    # Test with a comprehensive request that lets the agent decide approach
    user_request = """I need to study the FeMo-cofactor, which is the iron-molybdenum
cofactor [Fe7MoS9C] found in nitrogenase enzymes that catalyze nitrogen fixation.

This is a complex metal cluster with:
- 7 iron atoms
- 1 molybdenum atom
- 9 sulfur atoms
- 1 central carbide
- Homocitrate ligand

Please analyze this system using appropriate computational chemistry methods.
Since this is a large metal cluster, consider:
1. Starting with a simplified model if needed (e.g., Fe-S cluster or smaller fragment)
2. Using DFT methods suitable for transition metals
3. Classical MD if quantum methods are too computationally expensive

Please choose the most appropriate approach and run calculations to study
the electronic structure and energetics of this important catalytic center."""

    print(f"REQUEST:", flush=True)
    print(f"{user_request[:200]}...\n", flush=True)

    start = time.time()

    try:
        print(f"Running agentic workflow...\n", flush=True)
        result = agent.run_agentic_workflow(user_request)
        elapsed = time.time() - start

        print(f"\n{'='*80}", flush=True)
        print(f"RESULTS", flush=True)
        print(f"{'='*80}", flush=True)

        if result and result.get('status') == 'SUCCESS':
            res = result.get('result', {})
            print(f"✓ Success!", flush=True)
            print(f"  Total time: {elapsed:.1f} seconds", flush=True)
            print(f"  Application: {res.get('application', 'N/A')}", flush=True)
            print(f"  Experiment: {res.get('experiment_name', 'N/A')}", flush=True)

            results_data = res.get('results', {})
            if results_data:
                print(f"\n  Results:", flush=True)
                print(json.dumps(results_data, indent=4))

            return True
        else:
            print(f"✗ Did not complete with SUCCESS status", flush=True)
            print(f"  Status: {result.get('status', 'UNKNOWN')}", flush=True)
            print(f"  Time: {elapsed:.1f} seconds", flush=True)
            if result and 'result' in result:
                print(f"\n  Details:", flush=True)
                print(json.dumps(result.get('result'), indent=4))
            return False

    except Exception as e:
        elapsed = time.time() - start
        print(f"\n{'='*80}", flush=True)
        print(f"ERROR", flush=True)
        print(f"{'='*80}", flush=True)
        print(f"✗ Exception: {e}", flush=True)
        print(f"  Time: {elapsed:.1f} seconds", flush=True)
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_femo_cofactor()
    sys.exit(0 if success else 1)
