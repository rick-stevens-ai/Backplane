#!/usr/bin/env python3
"""
EXPANDED NH3 catalyst screening - 8 diverse candidates with corrected SMILES conversion
Original 5 + 3 additional molecules
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import time
import json
from datetime import datetime

def main():
    """Expanded NH3 catalyst screening with 8 molecules"""

    print(f"{'='*80}")
    print("EXPANDED NH3 CATALYST SCREENING (8 Candidates)")
    print(f"{'='*80}")
    print(f"Objective: Screen diverse NH3 synthesis catalyst candidates")
    print(f"Method: GPAW DFT with RDKit SMILES-to-3D conversion")
    print(f"Model: gpt-oss:120b")
    print(f"\\nOriginal 5 + 3 additional diverse molecules\\n")

    # Initialize agent
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    # 8 diverse catalyst candidates
    catalysts = [
        # ORIGINAL 5 MOLECULES
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
        },
        # NEW 3 MOLECULES
        {
            "name": "Pyrazole",
            "smiles": "c1c[nH]nc1",
            "formula": "C3H4N2",
            "expected_atoms": 9,
            "rationale": "Imidazole isomer, different N positions, metal coordination"
        },
        {
            "name": "1,2,4-Triazole",
            "smiles": "c1ncn[nH]1",
            "formula": "C2H3N3",
            "expected_atoms": 8,
            "rationale": "Three nitrogen atoms, potential N2 activation, multiple binding modes"
        },
        {
            "name": "Methylamine",
            "smiles": "CN",
            "formula": "CH5N",
            "expected_atoms": 7,
            "rationale": "Simplest primary amine, baseline for alkyl substitution effects"
        }
    ]

    print(f"Testing {len(catalysts)} catalyst candidates:\\n")
    print(f"{'#':<4}{'Name':<20}{'Formula':<10}{'Atoms':<7}{'Rationale'}")
    print("-" * 95)
    for i, cat in enumerate(catalysts, 1):
        status = "ORIGINAL" if i <= 5 else "NEW"
        print(f"{i:<4}{cat['name']:<20}{cat['formula']:<10}{cat['expected_atoms']:<7}{cat['rationale']}")

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
        print(f"SCREENING {i}/{len(catalysts)}: {name} {'[ORIGINAL]' if i <= 5 else '[NEW]'}")
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
                print(f"  Energy: {res_data.get('energy_ev', 'N/A'):.6f} eV")
                print(f"  Energy: {res_data.get('energy_hartree', 'N/A'):.8f} Hartree")

                # Verify atom count
                positions = res_data.get('positions', [])
                actual_atoms = len(positions)
                match_symbol = "✓" if actual_atoms == expected_atoms else "⚠"
                print(f"  Atoms: {actual_atoms}/{expected_atoms} {match_symbol}")

                if actual_atoms != expected_atoms:
                    print(f"  ⚠ WARNING: Atom count mismatch!")

                # Dipole moment magnitude
                dipole = res_data.get('dipole_moment', [0, 0, 0])
                dipole_mag = (dipole[0]**2 + dipole[1]**2 + dipole[2]**2)**0.5
                print(f"  Dipole: {dipole_mag:.3f} e·Å")

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
                    'dipole_moment': dipole,
                    'dipole_magnitude': dipole_mag,
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
    print("EXPANDED SCREENING SUMMARY")
    print(f"{'='*80}\\n")

    success_count = sum(1 for r in results if r.get('success', False))

    print(f"Total candidates: {len(results)}")
    print(f"  Original set: 5")
    print(f"  Additional: 3")
    print(f"Successful: {success_count}")
    print(f"Failed: {len(results) - success_count}")
    print(f"Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    print(f"Avg time/molecule: {total_time/len(results):.1f} seconds\\n")

    # Rank by energy
    successful = [r for r in results if r.get('success') and r.get('energy_ev')]
    successful.sort(key=lambda x: x['energy_ev'])

    if successful:
        print("\\nRANKED RESULTS (by total energy):")
        print(f"{'Rank':<6}{'Name':<20}{'Formula':<10}{'Atoms':<8}{'Energy (eV)':<16}{'Dipole (e·Å)':<14}{'Time (s)':<10}")
        print("-" * 100)

        for i, cat in enumerate(successful, 1):
            e_ev = cat['energy_ev']
            actual = cat.get('actual_atoms', '?')
            expected = cat['expected_atoms']
            atom_str = f"{actual}/{expected}"
            dipole_mag = cat.get('dipole_magnitude', 0.0)
            print(f"{i:<6}{cat['name']:<20}{cat['formula']:<10}{atom_str:<8}{e_ev:<16.6f}{dipole_mag:<14.3f}{cat['time']:<10.1f}")

        # Energy range analysis
        e_min = min(r['energy_ev'] for r in successful)
        e_max = max(r['energy_ev'] for r in successful)
        e_range = e_max - e_min

        print(f"\\n{'='*80}")
        print("ENERGY ANALYSIS:")
        print(f"{'='*80}")
        print(f"Lowest energy:  {e_min:.6f} eV ({successful[0]['name']})")
        print(f"Highest energy: {e_max:.6f} eV ({successful[-1]['name']})")
        print(f"Energy range:   {e_range:.6f} eV ({e_range*27.2114:.4f} Hartree)")
        print(f"\\nNote: Lower energy indicates more stable isolated molecule.")
        print(f"For catalysis, consider HOMO-LUMO gaps and binding energies.\\n")

    # Save results
    output = {
        'metadata': {
            'timestamp': datetime.now().isoformat(),
            'total_candidates': len(results),
            'original_molecules': 5,
            'additional_molecules': 3,
            'successful': success_count,
            'total_time_seconds': total_time,
            'model': 'gpt-oss:120b',
            'method': 'GPAW PBE finite-difference',
            'grid_spacing_angstrom': 0.2,
            'functional': 'PBE',
            'smiles_conversion': 'RDKit with ETKDG + UFF optimization'
        },
        'results': results
    }

    with open('nh3_expanded_results.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\\n{'='*80}")
    print(f"✓ Results saved to: nh3_expanded_results.json")
    print(f"{'='*80}")

    print(f"\\n{'='*80}")
    print("EXPANDED SCREENING COMPLETE!")
    print(f"{'='*80}\\n")

    return success_count > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
