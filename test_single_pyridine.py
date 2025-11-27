#!/usr/bin/env python3
"""
Quick test: Single Pyridine calculation through full stack to verify RDKit fix
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time

def main():
    print("Testing single Pyridine calculation to verify RDKit fix...")
    print("Expected: 11 atoms (C5H5N)")
    print("="*80)

    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    request = """Calculate the energy of Pyridine (SMILES: c1ccncc1) using GPAW.
Use finite-difference mode, PBE functional, h=0.2."""

    start = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed = time.time() - start

    if result and result.get('status') == 'SUCCESS':
        res_data = result.get('result', {}).get('results', {})
        positions = res_data.get('positions', [])
        actual_atoms = len(positions)

        print(f"\n{'='*80}")
        print(f"RESULT:")
        print(f"{'='*80}")
        print(f"Status: SUCCESS ({elapsed:.1f}s)")
        print(f"Atoms: {actual_atoms} (expected: 11)")
        print(f"Energy: {res_data.get('energy_ev', 'N/A')} eV")

        if actual_atoms == 11:
            print(f"\n✓✓✓ RDKIT FIX WORKING! Pyridine has correct 11 atoms! ✓✓✓")

            # Check input file
            job_id = result.get('result', {}).get('job_id')
            print(f"\nJob ID: {job_id}")
            return True
        else:
            print(f"\n✗✗✗ STILL BROKEN! Got {actual_atoms} atoms instead of 11 ✗✗✗")
            return False
    else:
        print(f"\n✗ Calculation FAILED")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
