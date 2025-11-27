#!/usr/bin/env python3
"""
CORRECTED NH3 catalyst screening - 5 diverse candidates with proper SMILES conversion
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time
import json
from datetime import datetime

def main():
    """Corrected NH3 catalyst screening with RDKit-based SMILES conversion"""

    print(f"{'='*80}")
    print("CORRECTED NH3 CATALYST SCREENING (5 Candidates)")
    print(f"{'='*80}")
    print(f"Objective: Rapid screening of diverse NH3 synthesis catalyst candidates")
    print(f"Method: GPAW DFT with FIXED SMILES-to-3D conversion (RDKit)")
    print(f"Model: gpt-oss:120b")
    print(f"Fix: Replaced hardcoded H2O structures with proper molecular geometries\\n")

    # Initialize agent
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    # 5 diverse, fast-to-compute catalyst candidates
    catalysts = [
        {
            "name": "Pyridine",
            "smiles": "c1ccncc1",
            "formula": "C5H5N",
            "expected_atoms": 11,
            "rationale": "Simple N-heterocycle, Lewis base, coordinates to metals"
        },
        {
            "name": "Imidazole",
            "smiles": "c1cnc[nH]1",
            "formula": "C3H4N2",
            "expected_atoms": 9,
            "rationale": "Enzyme-inspired, multiple coordination sites, proton shuttle"
        },
        {
            "name": "2,2'-Bipyridine",
            "smiles": "c1ccnc(c1)c2ccccn2",
            "formula": "C10H8N2",
            "expected_atoms": 20,
            "rationale": "Bidentate chelator, stabilizes metal complexes, common ligand"
        },
        {
            "name": "Ammonia",
            "smiles": "N",
            "formula": "NH3",
            "expected_atoms": 4,
            "rationale": "Product molecule baseline, understand NH3 energetics"
        },
        {
            "name": "Hydrazine",
            "smiles": "NN",
            "formula": "N2H4",
            "expected_atoms": 6,
            "rationale": "N-N bond model, intermediate in N2 reduction pathway"
        }
    ]

    print(f"Testing {len(catalysts)} catalyst candidates:\\n")
    for i, cat in enumerate(catalysts, 1):
        print(f"{i}. {cat['name']:20s} ({cat['formula']:8s}, {cat['expected_atoms']:2d} atoms) - {cat['rationale']}")

    print(f"\\n{'='*80}\\n")

    # Run screening
    results = []
    total_start = time.time()

    for i, catalyst in enumerate(catalysts, 1):
        name = catalyst['name']
        smiles = catalyst['smiles']
        formula = catalyst['formula']
        expected_atoms = catalyst['expected_atoms']
        rationale = catalyst['rationale']

        print(f"\\n{'='*80}")
        print(f"SCREENING {i}/{len(catalysts)}: {name}")
        print(f"{'='*80}")
        print(f"Formula: {formula}")
        print(f"SMILES: {smiles}")
        print(f"Expected atoms: {expected_atoms}")
        print(f"Rationale: {rationale}\\n")

        request = f"""Calculate the energy and electronic properties of {name} (SMILES: {smiles})
using GPAW with finite-difference mode and PBE functional. This is for NH3 catalyst screening.
Calculate:
- Total energy
- Electronic structure
- Forces
- Dipole moment
Use grid spacing h=0.2 for speed."""

        start = time.time()

        try:
            result = agent.run_agentic_workflow(request)
            elapsed = time.time() - start

            if result and result.get('status') == 'SUCCESS':
                res_data = result.get('result', {}).get('results', {})

                print(f"✓ SUCCESS ({elapsed:.1f}s)")
                print(f"  Energy: {res_data.get('energy_ev', 'N/A')} eV")
                print(f"  Energy: {res_data.get('energy_hartree', 'N/A')} Hartree")

                # Verify atom count
                positions = res_data.get('positions', [])
                actual_atoms = len(positions)
                if actual_atoms != expected_atoms:
                    print(f"  ⚠ WARNING: Expected {expected_atoms} atoms, got {actual_atoms}")

                results.append({
                    'name': name,
                    'formula': formula,
                    'smiles': smiles,
                    'expected_atoms': expected_atoms,
                    'actual_atoms': actual_atoms,
                    'rationale': rationale,
                    'success': True,
                    'time': elapsed,
                    'energy_ev': res_data.get('energy_ev'),
                    'energy_hartree': res_data.get('energy_hartree'),
                    'dipole_moment': res_data.get('dipole_moment'),
                    'results': res_data
                })
            else:
                print(f"✗ FAILED ({elapsed:.1f}s)")
                results.append({
                    'name': name,
                    'formula': formula,
                    'smiles': smiles,
                    'expected_atoms': expected_atoms,
                    'rationale': rationale,
                    'success': False,
                    'time': elapsed,
                    'status': result.get('status', 'UNKNOWN') if result else 'NO_RESULT'
                })

        except Exception as e:
            elapsed = time.time() - start
            print(f"✗ EXCEPTION: {e} ({elapsed:.1f}s)")
            results.append({
                'name': name,
                'formula': formula,
                'smiles': smiles,
                'expected_atoms': expected_atoms,
                'rationale': rationale,
                'success': False,
                'time': elapsed,
                'error': str(e)
            })

        time.sleep(1)

    total_time = time.time() - total_start

    # Generate report
    print(f"\\n\\n{'='*80}")
    print("SCREENING SUMMARY")
    print(f"{'='*80}\\n")

    success_count = sum(1 for r in results if r.get('success', False))

    print(f"Total candidates: {len(results)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {len(results) - success_count}")
    print(f"Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    print(f"Avg time/molecule: {total_time/len(results):.1f} seconds\\n")

    # Rank by energy
    successful = [r for r in results if r.get('success') and r.get('energy_ev')]
    successful.sort(key=lambda x: x['energy_ev'])

    if successful:
        print("\\nRANKED RESULTS (by total energy):")
        print(f"{'Rank':<6}{'Name':<20}{'Formula':<10}{'Atoms':<7}{'Energy (eV)':<15}{'Energy (Ha)':<15}{'Time (s)':<10}")
        print("-" * 95)

        for i, cat in enumerate(successful, 1):
            e_ev = cat['energy_ev']
            e_ha = cat.get('energy_hartree', e_ev/27.2114)
            actual = cat.get('actual_atoms', '?')
            expected = cat['expected_atoms']
            atom_str = f"{actual}/{expected}"
            print(f"{i:<6}{cat['name']:<20}{cat['formula']:<10}{atom_str:<7}{e_ev:<15.4f}{e_ha:<15.6f}{cat['time']:<10.1f}")

    # Comparison with original (invalid) results
    print(f"\\n{'='*80}")
    print("COMPARISON WITH ORIGINAL (INVALID) RESULTS:")
    print(f"{'='*80}")
    print("Original screening (nh3_fast_results.json):")
    print("  - All molecules: -14.614 eV (all calculated as H2O)")
    print("  - All geometries: Identical 3-atom structures")
    print("\\nCorrected screening (current run):")
    print("  - Distinct energies expected for different molecules")
    print("  - Correct atom counts (11, 9, 20, 4, 6 atoms)")
    print(f"{'='*80}\\n")

    # Save results
    output = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'total_candidates': len(results),
            'successful': success_count,
            'total_time_seconds': total_time,
            'model': 'gpt-oss:120b',
            'method': 'GPAW PBE',
            'fix_applied': 'RDKit SMILES-to-3D conversion',
            'previous_results': 'nh3_fast_results.json (INVALID - all H2O)'
        },
        'results': results
    }

    with open('nh3_corrected_results.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\\n✓ Results saved to: nh3_corrected_results.json")

    print(f"\\n{'='*80}")
    print("CORRECTED SCREENING COMPLETE!")
    print(f"{'='*80}\\n")

    return success_count > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
