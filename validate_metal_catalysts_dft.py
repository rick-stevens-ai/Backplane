#!/usr/bin/env python3
"""
DFT Validation of Top Metal Catalysts
Compare MACE-MP predictions with high-accuracy CP2K DFT
"""
import sys
import time
import json
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent

def load_xyz(filename):
    """Load XYZ file"""
    with open(filename, 'r') as f:
        content = f.read()
    # Extract just the structure part for conversion
    lines = content.strip().split('\n')
    n_atoms = int(lines[0])
    return content, n_atoms

def xyz_to_simple_format(xyz_content):
    """Convert XYZ to simple format the agent can understand"""
    lines = xyz_content.strip().split('\n')
    atoms = []
    for line in lines[2:]:  # Skip count and comment
        if line.strip():
            parts = line.split()
            if len(parts) >= 4:
                atoms.append({
                    'element': parts[0],
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'z': float(parts[3])
                })
    return atoms

def main():
    print("="*80)
    print("DFT VALIDATION OF TOP METAL CATALYSTS")
    print("="*80)
    print("Validating MACE-MP predictions with CP2K DFT")
    print("Top candidate: Ru₁₀ Cluster")
    print("="*80)

    # Initialize agent
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    # MACE-MP predictions (from screening)
    mace_results = {
        "Ru₁₀ Cluster": -57.217,
        "Ru-Ba/oxide": -42.285,
        "Fe-K-AlOx": -38.605,
        "Fe₃O₄": -8.429
    }

    print("\nMACE-MP Predictions (Materials Project model):")
    print("-"*80)
    for name, energy in mace_results.items():
        print(f"  {name:20s} {energy:10.3f} eV")

    # Validate Ru₁₀ with CP2K
    print("\n" + "="*80)
    print("VALIDATING Ru₁₀ CLUSTER WITH CP2K DFT")
    print("="*80)

    xyz_content, n_atoms = load_xyz("ru10_cluster.xyz")
    print(f"\nStructure: Ru₁₀ cluster")
    print(f"Atoms: {n_atoms} Ru atoms")
    print(f"MACE-MP prediction: {mace_results['Ru₁₀ Cluster']:.3f} eV")
    print("\nSubmitting to CP2K for high-accuracy DFT calculation...")
    print("Method: PBE functional, DZVP basis")

    # Create request for CP2K
    # Note: CP2K needs XYZ format, we'll provide it directly
    request = f"""Calculate the energy of a Ru₁₀ metal cluster using CP2K.
This is a ruthenium cluster representing a Ru(0001) terrace fragment for NH₃ synthesis.

Structure (XYZ format):
{xyz_content}

Use:
- Run type: ENERGY (single point energy calculation)
- Functional: PBE
- Basis: DZVP (good for transition metals)
- This is a metal cluster, ensure appropriate settings for metallic system

This is for validating MACE-MP predictions with high-accuracy DFT."""

    print("\n" + "-"*80)
    print("CP2K Calculation Starting...")
    print("-"*80)

    start = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed = time.time() - start

    print("\n" + "="*80)
    print("CP2K VALIDATION RESULTS")
    print("="*80)

    if result and result.get('status') in ['SUCCESS', 'COMPLETED']:
        # Extract energy from results
        res_data = result.get('result', {})
        if isinstance(res_data, dict):
            results_dict = res_data.get('results', {})
            cp2k_energy = results_dict.get('energy_ev', results_dict.get('total_energy'))

            if cp2k_energy:
                mace_energy = mace_results['Ru₁₀ Cluster']
                difference = abs(cp2k_energy - mace_energy)
                percent_error = (difference / abs(cp2k_energy)) * 100 if cp2k_energy != 0 else 0

                print(f"\n✓ SUCCESS!")
                print(f"\nComparison:")
                print(f"  MACE-MP:     {mace_energy:12.3f} eV")
                print(f"  CP2K DFT:    {cp2k_energy:12.3f} eV")
                print(f"  Difference:  {difference:12.3f} eV")
                print(f"  Error:       {percent_error:12.2f}%")
                print(f"\nCalculation time:")
                print(f"  MACE-MP:     2.2 seconds")
                print(f"  CP2K:        {elapsed:.1f} seconds")
                print(f"  Speedup:     {elapsed/2.2:.0f}x faster with MACE")

                # Save results
                validation_results = {
                    "catalyst": "Ru₁₀ Cluster",
                    "mace_mp_energy_ev": mace_energy,
                    "cp2k_dft_energy_ev": cp2k_energy,
                    "difference_ev": difference,
                    "percent_error": percent_error,
                    "mace_time_s": 2.2,
                    "cp2k_time_s": elapsed,
                    "speedup": elapsed/2.2,
                    "method": {
                        "mace": "MACE-MP medium, Materials Project model",
                        "dft": "CP2K PBE/DZVP"
                    }
                }

                with open('metal_catalyst_validation.json', 'w') as f:
                    json.dump(validation_results, f, indent=2)
                print(f"\n✓ Results saved to: metal_catalyst_validation.json")

                # Conclusion
                print("\n" + "="*80)
                print("VALIDATION CONCLUSION")
                print("="*80)
                if percent_error < 10:
                    print(f"✓ EXCELLENT: MACE-MP within {percent_error:.1f}% of DFT")
                    print("  MACE-MP provides reliable screening for metal catalysts")
                elif percent_error < 20:
                    print(f"✓ GOOD: MACE-MP within {percent_error:.1f}% of DFT")
                    print("  Suitable for rapid pre-screening, validate top candidates")
                else:
                    print(f"⚠ MODERATE: MACE-MP differs by {percent_error:.1f}%")
                    print("  Use for qualitative ranking, always validate with DFT")

                print(f"\nWorkflow validated:")
                print(f"  1. MACE-MP rapid screening: ~2s per catalyst")
                print(f"  2. Rank candidates by predicted energy")
                print(f"  3. DFT validation of top picks: ~{elapsed/60:.0f} min each")
                print(f"  4. Total time: seconds (MACE) + hours (DFT top candidates)")
                print(f"\n  For 100 catalysts: 200s (MACE) + {elapsed*5/60:.0f} min (DFT top 5)")
                print(f"  vs all-DFT: {elapsed*100/3600:.1f} hours!")

            else:
                print("⚠ Could not extract energy from CP2K results")
                print(f"Result structure: {res_data}")
        else:
            print(f"⚠ Unexpected result format: {result}")
    else:
        print(f"✗ CP2K calculation failed or incomplete")
        print(f"Status: {result.get('status') if result else 'No result'}")
        print("\nNote: This may be expected if CP2K needs specific input format")
        print("      or if the calculation is still running")

    print("\n" + "="*80)
    print("✓ Validation workflow complete!")
    print("="*80)

if __name__ == "__main__":
    main()
