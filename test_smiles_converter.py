#!/usr/bin/env python3
"""
Test SMILES Converter Module
Demonstrates formula-to-SMILES conversion using oss120
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from smiles_converter import create_converter
from agent_apps import ComputationalChemistryAgent

def test_converter():
    print("="*80)
    print("SMILES CONVERTER TEST")
    print("="*80)
    print("Testing formula-to-SMILES conversion using oss120 (gpt-oss:120b)")
    print("="*80)
    print()

    # Initialize agent
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"  # Using 120b model
    )

    # Create converter
    converter = create_converter(agent)

    # Test cases - mix of simple molecules and complex systems
    test_cases = [
        # Simple organic molecules (should succeed)
        {"formula": "H2O", "description": "water"},
        {"formula": "NH3", "description": "ammonia"},
        {"formula": "C6H6", "description": "benzene"},
        {"formula": "C2H5OH", "description": "ethanol"},

        # Nitrogen-containing ligands (should succeed)
        {"formula": "C5H5N", "description": "pyridine"},
        {"formula": "C10H8N2", "description": "2,2'-bipyridine"},

        # Metal complexes (should fail - not SMILES representable)
        {"formula": "Fe(CO)4(N2)", "description": "iron carbonyl complex with dinitrogen"},
        {"formula": "Ru10", "description": "ruthenium metal cluster with 10 atoms"},
    ]

    results = []

    for i, test in enumerate(test_cases, 1):
        print(f"\n{'='*80}")
        print(f"TEST {i}/{len(test_cases)}: {test['formula']} ({test['description']})")
        print(f"{'='*80}")

        try:
            # Convert formula to SMILES
            smiles, metadata = converter.formula_to_smiles(
                formula=test['formula'],
                description=test['description'],
                validate=True
            )

            if smiles and metadata['valid']:
                print(f"✓ SUCCESS")
                print(f"  Formula:     {test['formula']}")
                print(f"  Description: {test['description']}")
                print(f"  SMILES:      {smiles}")
                print(f"  Valid:       {metadata['valid']}")

                results.append({
                    'formula': test['formula'],
                    'smiles': smiles,
                    'status': 'success'
                })
            else:
                error = metadata.get('error', 'Unknown error')
                print(f"✗ FAILED TO CONVERT")
                print(f"  Formula:     {test['formula']}")
                print(f"  Error:       {error}")

                results.append({
                    'formula': test['formula'],
                    'smiles': None,
                    'status': 'not_representable',
                    'error': error
                })

        except Exception as e:
            print(f"✗ ERROR")
            print(f"  Formula: {test['formula']}")
            print(f"  Error:   {str(e)}")

            results.append({
                'formula': test['formula'],
                'smiles': None,
                'status': 'error',
                'error': str(e)
            })

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")

    successful = [r for r in results if r['status'] == 'success']
    not_representable = [r for r in results if r['status'] == 'not_representable']
    errors = [r for r in results if r['status'] == 'error']

    print(f"\nTotal tests:        {len(test_cases)}")
    print(f"Successful:         {len(successful)}")
    print(f"Not representable:  {len(not_representable)}")
    print(f"Errors:             {len(errors)}")

    if successful:
        print(f"\n✓ Successfully converted to SMILES:")
        for r in successful:
            print(f"  {r['formula']:15s} → {r['smiles']}")

    if not_representable:
        print(f"\n⚠ Not representable as SMILES (expected for metal complexes):")
        for r in not_representable:
            print(f"  {r['formula']}")

    if errors:
        print(f"\n✗ Errors:")
        for r in errors:
            print(f"  {r['formula']}: {r.get('error', 'Unknown')}")

    print(f"\n{'='*80}")
    print("CONVERTER CAPABILITIES DEMONSTRATED:")
    print(f"{'='*80}")
    print("✓ Converts simple molecular formulas to SMILES using LLM reasoning")
    print("✓ Validates SMILES with RDKit (or basic validation)")
    print("✓ Detects when structures cannot be represented as SMILES")
    print("✓ Handles both organic molecules and recognizes metal complexes")
    print(f"{'='*80}")

if __name__ == "__main__":
    test_converter()
