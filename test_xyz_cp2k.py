#!/usr/bin/env python3
"""
Test CP2K wrapper with XYZ input (Ru10 metal cluster)
"""
import sys
import time
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent

def load_xyz(filename):
    """Load XYZ file"""
    with open(filename, 'r') as f:
        return f.read()

def main():
    print("="*80)
    print("TESTING XYZ INPUT WITH CP2K WRAPPER")
    print("="*80)
    print("Test: Ru₁₀ metal cluster DFT calculation")
    print("="*80)

    # Load Ru10 XYZ structure
    xyz_content = load_xyz("ru10_cluster.xyz")

    print("\nXYZ Structure:")
    print("-"*80)
    print(xyz_content[:200] + "...")
    print("-"*80)

    # Initialize agent
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    # Create request with XYZ structure
    request = f"""Calculate the energy of a Ru₁₀ metal cluster using CP2K with XYZ coordinates.

This is a ruthenium cluster representing a Ru(0001) terrace fragment.

XYZ Structure:
{xyz_content}

Settings:
- Run type: ENERGY (single point energy calculation)
- Functional: PBE
- Basis: DZVP (appropriate for transition metals)
- This is a metal cluster system (non-periodic)

Please calculate the total energy."""

    print("\nSubmitting to Agent...")
    print("-"*80)

    start = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed = time.time() - start

    print("\n" + "="*80)
    print("RESULTS")
    print("="*80)

    if result and result.get('status') in ['SUCCESS', 'COMPLETED']:
        print("✓ SUCCESS!")
        print(f"\nCalculation completed in {elapsed:.1f} seconds")

        # Extract energy from results
        res_data = result.get('result', {})
        if isinstance(res_data, dict):
            results_dict = res_data.get('results', {})

            energy_ev = results_dict.get('energy_ev')
            energy_hartree = results_dict.get('energy_hartree')

            if energy_ev:
                print(f"\nEnergy Results:")
                print(f"  CP2K Energy:     {energy_ev:.3f} eV")
                if energy_hartree:
                    print(f"  CP2K Energy:     {energy_hartree:.6f} Hartree")

                # Compare with MACE-MP prediction
                mace_energy = -57.217  # From MACE-MP screening
                difference = abs(energy_ev - mace_energy)
                percent_error = (difference / abs(energy_ev)) * 100 if energy_ev != 0 else 0

                print(f"\nMAC-MP vs DFT Comparison:")
                print(f"  MACE-MP:         {mace_energy:.3f} eV (2.2s)")
                print(f"  CP2K DFT:        {energy_ev:.3f} eV ({elapsed:.1f}s)")
                print(f"  Difference:      {difference:.3f} eV")
                print(f"  Error:           {percent_error:.2f}%")
                print(f"  Speedup:         {elapsed/2.2:.0f}x faster with MACE")

                print("\n" + "="*80)
                print("✓ XYZ INPUT EXTENSION WORKS!")
                print("="*80)
                print("✓ CP2K successfully processed XYZ coordinates for metal cluster")
                print("✓ Wrappers can now handle non-organic systems")
                print("✓ Full metal catalyst workflow operational")
            else:
                print("⚠ Could not extract energy from results")
                print(f"Results: {results_dict}")
        else:
            print(f"⚠ Unexpected result format: {result}")
    else:
        print(f"✗ Calculation failed or incomplete")
        print(f"Status: {result.get('status') if result else 'No result'}")
        if result:
            print(f"Result: {result}")

    print("\n" + "="*80)

if __name__ == "__main__":
    main()
