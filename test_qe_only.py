#!/usr/bin/env python3
"""Quick test: Quantum ESPRESSO energy calculation for caffeine"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time

CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

print("="*80, flush=True)
print("QUANTUM ESPRESSO - Caffeine Energy Calculation Test", flush=True)
print("="*80, flush=True)

print("\nInitializing agent...", flush=True)
agent = ComputationalChemistryAgent(
    server_config_path="spark_servers.yaml",
    server_name="spark-container-03"
)
print(f"✓ Agent ready: {agent.model}", flush=True)

print("\nSubmitting Quantum ESPRESSO job for caffeine molecule...", flush=True)
print(f"SMILES: {CAFFEINE_SMILES}", flush=True)

request = f"""Please use Quantum ESPRESSO to calculate the total electronic energy
of caffeine molecule (SMILES: {CAFFEINE_SMILES}).

Perform an SCF (self-consistent field) calculation to get the ground state energy.
Use standard cutoffs for accuracy. This will give us the DFT energy of the system."""

start = time.time()
result = agent.run_agentic_workflow(request)
elapsed = time.time() - start

print("\n" + "="*80, flush=True)
print("RESULTS", flush=True)
print("="*80, flush=True)

if result and result.get('status') == 'SUCCESS':
    res = result.get('result', {})
    print(f"✓ Calculation completed successfully!", flush=True)
    print(f"  Time: {elapsed:.1f} seconds", flush=True)
    print(f"  Application: {res.get('application', 'N/A')}", flush=True)
    print(f"  Experiment: {res.get('experiment_name', 'N/A')}", flush=True)

    results_data = res.get('results', {})
    if results_data:
        print("\nEnergy Results:", flush=True)
        for key, value in results_data.items():
            print(f"  {key}: {value}", flush=True)
else:
    print(f"✗ Calculation failed", flush=True)
    print(f"  Status: {result.get('status', 'UNKNOWN')}", flush=True)

print("\n" + "="*80, flush=True)
print("Test complete!", flush=True)
