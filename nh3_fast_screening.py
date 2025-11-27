#!/usr/bin/env python3
"""
Fast NH3 catalyst screening - 5 diverse candidates
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time
import json
from datetime import datetime

def main():
    """Fast NH3 catalyst screening"""

    print(f"{'='*80}")
    print("FAST NH3 CATALYST SCREENING (5 Candidates)")
    print(f"{'='*80}")
    print(f"Objective: Rapid screening of diverse NH3 synthesis catalyst candidates")
    print(f"Method: GPAW DFT for speed")
    print(f"Model: gpt-oss:120b\n")

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
            "rationale": "Simple N-heterocycle, Lewis base, coordinates to metals"
        },
        {
            "name": "Imidazole",
            "smiles": "c1cnc[nH]1",
            "rationale": "Enzyme-inspired, multiple coordination sites, proton shuttle"
        },
        {
            "name": "2,2'-Bipyridine",
            "smiles": "c1ccnc(c1)c2ccccn2",
            "rationale": "Bidentate chelator, stabilizes metal complexes, common ligand"
        },
        {
            "name": "Ammonia",
            "smiles": "N",
            "rationale": "Product molecule baseline, understand NH3 energetics"
        },
        {
            "name": "Hydrazine",
            "smiles": "NN",
            "rationale": "N-N bond model, intermediate in N2 reduction pathway"
        }
    ]

    print(f"Testing {len(catalysts)} catalyst candidates:\n")
    for i, cat in enumerate(catalysts, 1):
        print(f"{i}. {cat['name']:20s} - {cat['rationale']}")

    print(f"\n{'='*80}\n")

    # Run screening
    results = []
    total_start = time.time()

    for i, catalyst in enumerate(catalysts, 1):
        name = catalyst['name']
        smiles = catalyst['smiles']
        rationale = catalyst['rationale']

        print(f"\n{'='*80}")
        print(f"SCREENING {i}/{len(catalysts)}: {name}")
        print(f"{'='*80}")
        print(f"SMILES: {smiles}")
        print(f"Rationale: {rationale}\n")

        request = f"""Calculate the energy and electronic properties of {name} (SMILES: {smiles})
using GPAW with finite-difference mode and PBE functional. This is for NH3 catalyst screening.
Calculate:
- Total energy
- Electronic structure
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

                results.append({
                    'name': name,
                    'smiles': smiles,
                    'rationale': rationale,
                    'success': True,
                    'time': elapsed,
                    'energy_ev': res_data.get('energy_ev'),
                    'energy_hartree': res_data.get('energy_hartree'),
                    'results': res_data
                })
            else:
                print(f"✗ FAILED ({elapsed:.1f}s)")
                results.append({
                    'name': name,
                    'smiles': smiles,
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
                'smiles': smiles,
                'rationale': rationale,
                'success': False,
                'time': elapsed,
                'error': str(e)
            })

        time.sleep(1)

    total_time = time.time() - total_start

    # Generate report
    print(f"\n\n{'='*80}")
    print("SCREENING SUMMARY")
    print(f"{'='*80}\n")

    success_count = sum(1 for r in results if r.get('success', False))

    print(f"Total candidates: {len(results)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {len(results) - success_count}")
    print(f"Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    print(f"Avg time/molecule: {total_time/len(results):.1f} seconds\n")

    # Rank by energy
    successful = [r for r in results if r.get('success') and r.get('energy_ev')]
    successful.sort(key=lambda x: x['energy_ev'])

    if successful:
        print("\nRANKED RESULTS (by total energy):")
        print(f"{'Rank':<6}{'Name':<20}{'Energy (eV)':<15}{'Energy (Ha)':<15}{'Time (s)':<10}")
        print("-" * 80)

        for i, cat in enumerate(successful, 1):
            e_ev = cat['energy_ev']
            e_ha = cat.get('energy_hartree', e_ev/27.2114)
            print(f"{i:<6}{cat['name']:<20}{e_ev:<15.4f}{e_ha:<15.6f}{cat['time']:<10.1f}")

    # Save results
    output = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'total_candidates': len(results),
            'successful': success_count,
            'total_time_seconds': total_time,
            'model': 'gpt-oss:120b',
            'method': 'GPAW PBE'
        },
        'results': results
    }

    with open('nh3_fast_results.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n✓ Results saved to: nh3_fast_results.json")

    print(f"\n{'='*80}")
    print("FAST SCREENING COMPLETE!")
    print(f"{'='*80}\n")

    return success_count > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
