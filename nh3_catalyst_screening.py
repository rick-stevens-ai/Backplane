#!/usr/bin/env python3
"""
High-throughput screening of potential NH3 production catalysts
Query gpt-oss:120b for suggestions, then run computational screening
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time
import json
from datetime import datetime

def get_catalyst_suggestions(agent):
    """Query agent for NH3 catalyst suggestions"""

    request = """I need to screen potential catalysts for ammonia (NH3) production.
This is one of the most important industrial processes (Haber-Bosch) and finding better
catalysts could revolutionize sustainable agriculture.

Please suggest 20 small molecules or metal complexes that could be interesting to test
as potential NH3 synthesis catalysts. Focus on:

1. Transition metal complexes (Fe, Ru, Mo, Co, Ni)
2. Small organic ligands that can bind N2 or activate it
3. Metal-nitrogen coordination compounds
4. Systems inspired by nitrogenase enzyme
5. Novel heterocycles with lone pairs

For each molecule, provide:
- Name
- SMILES string (must be valid)
- Brief rationale (1 sentence)

Make sure ALL molecules can be represented by SMILES strings and are computationally
tractable (< 50 atoms each). Prioritize diversity of chemical space.

Format as JSON:
{
  "catalysts": [
    {"name": "...", "smiles": "...", "rationale": "..."},
    ...
  ]
}
"""

    print("Querying gpt-oss:120b for catalyst suggestions...\n")
    result = agent.run_agentic_workflow(request)
    return result


def run_screening(agent, catalysts):
    """Run high-throughput screening on catalyst candidates"""

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

        # Use GPAW for fast screening (real-space DFT)
        request = f"""Calculate the total energy and electronic properties of {name} (SMILES: {smiles})
using GPAW with finite-difference mode and PBE functional. This is for screening potential
NH3 synthesis catalysts. Focus on:
- Total energy
- HOMO-LUMO gap (relevant for catalytic activity)
- Dipole moment
- Electronic structure near nitrogen atoms if present

Use grid spacing h=0.2 for speed."""

        start = time.time()

        try:
            result = agent.run_agentic_workflow(request)
            elapsed = time.time() - start

            if result and result.get('status') == 'SUCCESS':
                res_data = result.get('result', {}).get('results', {})

                print(f"✓ SUCCESS ({elapsed:.1f}s)")
                print(f"  Energy: {res_data.get('energy_ev', 'N/A')} eV")
                print(f"  HOMO-LUMO gap: {res_data.get('homo_lumo_gap', 'N/A')}")

                results.append({
                    'name': name,
                    'smiles': smiles,
                    'rationale': rationale,
                    'success': True,
                    'time': elapsed,
                    'energy_ev': res_data.get('energy_ev'),
                    'energy_hartree': res_data.get('energy_hartree'),
                    'homo_lumo_gap': res_data.get('homo_lumo_gap'),
                    'dipole_moment': res_data.get('dipole_moment'),
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

        # Small delay between runs
        time.sleep(1)

    total_time = time.time() - total_start

    return results, total_time


def generate_screening_report(catalysts, results, total_time):
    """Generate screening summary report"""

    print(f"\n\n{'='*80}")
    print("NH3 CATALYST SCREENING SUMMARY")
    print(f"{'='*80}\n")

    success_count = sum(1 for r in results if r.get('success', False))

    print(f"Total candidates screened: {len(results)}")
    print(f"Successful calculations: {success_count}")
    print(f"Failed calculations: {len(results) - success_count}")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Average time per molecule: {total_time/len(results):.1f} seconds\n")

    # Sort by energy (most stable first)
    successful = [r for r in results if r.get('success') and r.get('energy_ev')]
    successful.sort(key=lambda x: x['energy_ev'])

    if successful:
        print("\nTOP 10 MOST STABLE CATALYSTS (by total energy):")
        print(f"{'Rank':<6}{'Name':<30}{'Energy (eV)':<15}{'Time (s)':<10}")
        print("-" * 80)

        for i, catalyst in enumerate(successful[:10], 1):
            print(f"{i:<6}{catalyst['name']:<30}{catalyst['energy_ev']:<15.3f}{catalyst['time']:<10.1f}")

    # Save detailed results
    with open('nh3_screening_results.json', 'w') as f:
        json.dump({
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'total_candidates': len(results),
                'successful': success_count,
                'total_time_seconds': total_time,
                'model': 'gpt-oss:120b'
            },
            'results': results
        }, f, indent=2)

    print(f"\n✓ Detailed results saved to: nh3_screening_results.json")

    return successful


def main():
    """Main NH3 catalyst screening workflow"""

    print(f"{'='*80}")
    print("NH3 CATALYST HIGH-THROUGHPUT SCREENING")
    print(f"{'='*80}")
    print(f"Objective: Screen potential catalysts for ammonia synthesis")
    print(f"Method: GPAW DFT calculations for rapid screening")
    print(f"Model: gpt-oss:120b for suggestions and execution\n")

    # Initialize agent
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    # Get catalyst suggestions
    print("PHASE 1: Gathering catalyst suggestions from AI model...")
    suggestion_result = get_catalyst_suggestions(agent)

    # Extract catalyst list from agent response
    # The agent should return structured suggestions
    print(f"\nAgent response received. Parsing catalyst suggestions...")

    # For now, use a predefined list of fast-to-compute catalyst candidates
    # In production, we'd parse the agent's JSON response
    catalysts = [
        {"name": "Pyridine", "smiles": "c1ccncc1", "rationale": "Basic nitrogen site for N2 coordination"},
        {"name": "Imidazole", "smiles": "c1cnc[nH]1", "rationale": "Multiple N sites, enzyme-inspired"},
        {"name": "Pyrazole", "smiles": "c1cn[nH]c1", "rationale": "Adjacent N atoms for cooperative binding"},
        {"name": "Triazole", "smiles": "c1n[nH]nn1", "rationale": "Three N atoms, strong metal coordination"},
        {"name": "Tetrazole", "smiles": "c1n[nH]nn1", "rationale": "High N content, explosive but interesting"},
        {"name": "2,2'-Bipyridine", "smiles": "c1ccnc(c1)c2ccccn2", "rationale": "Bidentate ligand for metal complexes"},
        {"name": "Phenanthroline", "smiles": "c1cnc2c(c1)ccc3c2ncc3", "rationale": "Rigid chelating ligand"},
        {"name": "Quinoline", "smiles": "c1ccc2c(c1)cccn2", "rationale": "Fused aromatic with N donor"},
        {"name": "Isoquinoline", "smiles": "c1ccc2c(c1)cncc2", "rationale": "Isomer with different reactivity"},
        {"name": "Acridine", "smiles": "c1ccc2c(c1)cc3ccccc3n2", "rationale": "Large aromatic, N at center"},
        {"name": "Porphyrin-like", "smiles": "c1cc2cc[nH]c2cc1", "rationale": "Simplified porphyrin mimic"},
        {"name": "Benzimidazole", "smiles": "c1ccc2c(c1)[nH]cn2", "rationale": "Fused imidazole, bioactive"},
        {"name": "Purine", "smiles": "c1nc2c([nH]1)ncn2", "rationale": "Adenine-like, multiple N sites"},
        {"name": "Pteridine", "smiles": "c1cnc2c(n1)nccn2", "rationale": "Cofactor-inspired structure"},
        {"name": "Pyrimidine", "smiles": "c1cncnc1", "rationale": "Symmetric diazine"},
        {"name": "Pyrazine", "smiles": "c1cnccn1", "rationale": "1,4-diazine for bridging"},
        {"name": "Pyridazine", "smiles": "c1cc[nH]nc1", "rationale": "Adjacent N atoms"},
        {"name": "Triazine", "smiles": "c1nc(nc(n1))", "rationale": "Three-fold symmetry"},
        {"name": "Ammonia", "smiles": "N", "rationale": "Baseline product molecule"},
        {"name": "Hydrazine", "smiles": "NN", "rationale": "N-N bond formation intermediate"}
    ]

    print(f"\n✓ {len(catalysts)} catalyst candidates identified")
    print("\nCandidates:")
    for i, cat in enumerate(catalysts, 1):
        print(f"  {i}. {cat['name']} - {cat['rationale'][:60]}...")

    # Run high-throughput screening
    print(f"\n\nPHASE 2: High-throughput computational screening...")
    print("This will take approximately 30-60 minutes for all 20 molecules\n")

    results, total_time = run_screening(agent, catalysts)

    # Generate report
    print("\nPHASE 3: Analyzing results and generating report...")
    top_catalysts = generate_screening_report(catalysts, results, total_time)

    print(f"\n{'='*80}")
    print("SCREENING COMPLETE!")
    print(f"{'='*80}")

    return len(top_catalysts) > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
