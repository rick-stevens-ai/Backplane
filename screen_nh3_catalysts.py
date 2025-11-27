#!/usr/bin/env python3
"""
Screen NH3 catalyst candidates with MACE using XYZ structures
"""
import sys
import time
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent

# Metal complexes with simplified XYZ structures
# Format: name, xyz_structure, description
catalyst_structures = [
    {
        "name": "Fe(N2)(CO)4 anion",
        "xyz": """9
Fe-N2-CO4 complex
Fe  0.0  0.0  0.0
N   0.0  0.0  1.8
N   0.0  0.0  3.0
C   2.0  0.0  0.0
O   3.2  0.0  0.0
C  -2.0  0.0  0.0
O  -3.2  0.0  0.0
C   0.0  2.0  0.0
O   0.0  3.2  0.0""",
        "class": "homogeneous",
        "features": "Strong Fe-N2 back-bonding, labile CO ligands"
    },

    {
        "name": "Ru cluster (Ru10)",
        "xyz": """10
Ru10 cluster - Ru(0001) terrace fragment
Ru    0.00000000    0.00000000    0.00000000
Ru    2.70000000    0.00000000    0.00000000
Ru    1.35000000    2.33826859    0.00000000
Ru   -1.35000000    2.33826859    0.00000000
Ru   -2.70000000    0.00000000    0.00000000
Ru   -1.35000000   -2.33826859    0.00000000
Ru    1.35000000   -2.33826859    0.00000000
Ru    1.63678801    0.94500000    2.20000000
Ru   -1.63678801    0.94500000    2.20000000
Ru    0.00000000   -1.89000000    2.20000000""",
        "class": "cluster",
        "features": "Ru binds N2 strongly, step sites lower H2 activation barrier"
    },

    {
        "name": "Mo(N2)(CO)5 anion",
        "xyz": """12
Mo-N2-CO5 complex
Mo  0.0  0.0  0.0
N   0.0  0.0  1.9
N   0.0  0.0  3.1
C   2.1  0.0  0.0
O   3.3  0.0  0.0
C  -2.1  0.0  0.0
O  -3.3  0.0  0.0
C   0.0  2.1  0.0
O   0.0  3.3  0.0
C   0.0 -2.1  0.0
O   0.0 -3.3  0.0""",
        "class": "homogeneous",
        "features": "High dπ → π* back-donation, labile CO ligands"
    },

    {
        "name": "Co-N4 single site",
        "xyz": """9
Co-N4 site (Fe-N-C analog)
Co  0.0  0.0  0.0
N   1.9  0.0  0.0
N   0.0  1.9  0.0
N  -1.9  0.0  0.0
N   0.0 -1.9  0.0
C   2.9  0.0  0.0
C   0.0  2.9  0.0
C  -2.9  0.0  0.0
C   0.0 -2.9  0.0""",
        "class": "single-atom",
        "features": "Mimics nitrogenase active site, strong N2 binding"
    },

    {
        "name": "Fe2S2(CO)6 dimer",
        "xyz": """16
Fe2-S2-CO6 cluster
Fe  0.0  0.0  0.0
Fe  2.6  0.0  0.0
S   1.3  1.2  0.0
S   1.3 -1.2  0.0
C  -1.2  1.2  0.0
O  -2.0  2.0  0.0
C  -1.2 -1.2  0.0
O  -2.0 -2.0  0.0
C  -0.5  0.0  1.5
O  -0.8  0.0  2.7
C   3.8  1.2  0.0
O   4.6  2.0  0.0
C   3.8 -1.2  0.0
O   4.6 -2.0  0.0
C   3.1  0.0  1.5
O   3.4  0.0  2.7""",
        "class": "cluster",
        "features": "Fe-Fe bond + sulfide bridge mimics nitrogenase"
    },

    {
        "name": "Fe3 triangle",
        "xyz": """3
Fe3 triangle cluster
Fe  0.0  0.0  0.0
Fe  2.5  0.0  0.0
Fe  1.25  2.165  0.0""",
        "class": "cluster",
        "features": "Small Fe cluster for N2 activation"
    },

    {
        "name": "Ni-N4 pincer",
        "xyz": """9
Ni-N4 pincer site
Ni  0.0  0.0  0.0
N   2.0  0.0  0.0
N   0.0  2.0  0.0
N  -2.0  0.0  0.0
N   0.0 -2.0  0.0
C   3.0  0.0  0.0
C   0.0  3.0  0.0
C  -3.0  0.0  0.0
C   0.0 -3.0  0.0""",
        "class": "single-atom",
        "features": "Low-spin Ni binds N2 side-on, cheaper than Ru"
    }
]

def main():
    print("="*80)
    print("NH3 CATALYST SCREENING WITH MACE")
    print("="*80)
    print(f"Testing {len(catalyst_structures)} metal complexes/clusters")
    print("="*80)
    print()

    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"  # Using 120b for now
    )

    results = []

    for i, catalyst in enumerate(catalyst_structures, 1):
        print(f"\n{'='*80}")
        print(f"CATALYST {i}/{len(catalyst_structures)}: {catalyst['name']}")
        print(f"Class: {catalyst['class']}")
        print(f"Features: {catalyst['features']}")
        print(f"{'='*80}")

        # Use agent to run MACE screening with XYZ
        request = f"""Use MACE-MP to predict the energy of this metal complex:

{catalyst['name']} - {catalyst['features']}

XYZ Structure:
{catalyst['xyz']}

This is a {catalyst['class']} catalyst for NH3 formation. Use MACE's rapid energy prediction (should take ~1 second).

Report:
1. Total energy (eV)
2. Energy per atom (eV/atom)
3. Whether this looks promising based on the energy"""

        start = time.time()
        try:
            result = agent.run_agentic_workflow(request)
            elapsed = time.time() - start

            print(f"\n✓ MACE screening completed in {elapsed:.1f}s")

            # Extract energy if available
            if result and 'result' in result:
                response = result['result']
                if isinstance(response, dict) and 'message' in response:
                    message = response['message']
                    print(f"\nRESULT: {message[:300]}...")

                    results.append({
                        'name': catalyst['name'],
                        'class': catalyst['class'],
                        'time': elapsed,
                        'result': message,
                        'status': 'success'
                    })

        except Exception as e:
            elapsed = time.time() - start
            print(f"\n✗ Error after {elapsed:.1f}s: {str(e)[:200]}")
            results.append({
                'name': catalyst['name'],
                'class': catalyst['class'],
                'time': elapsed,
                'result': str(e),
                'status': 'error'
            })

        # Small delay between requests
        if i < len(catalyst_structures):
            time.sleep(2)

    # Summary
    print("\n" + "="*80)
    print("SCREENING SUMMARY")
    print("="*80)
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'error']

    print(f"\nTotal candidates: {len(catalyst_structures)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")

    if successful:
        print(f"\n✓ Successfully screened:")
        for r in successful:
            print(f"  - {r['name']} ({r['time']:.1f}s)")

    if failed:
        print(f"\n✗ Failed:")
        for r in failed:
            print(f"  - {r['name']}: {r['result'][:100]}")

    print("\n" + "="*80)
    print("Next step: Select top 3-5 candidates for DFT validation with CP2K")
    print("="*80)

if __name__ == "__main__":
    main()
