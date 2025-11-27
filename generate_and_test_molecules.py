#!/usr/bin/env python3
"""
Use gpt-oss:120b to generate 100 molecules, then evaluate with MACE
"""
import sys
import time
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent

def main():
    print("="*80)
    print("MOLECULE GENERATION AND MACE EVALUATION WORKFLOW")
    print("="*80)
    print()
    
    # Initialize agent
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )
    
    # Request oss120b to generate 100 molecules
    print("STEP 1: Requesting gpt-oss:120b to generate 100 molecules...")
    print("-"*80)
    
    generation_request = """Generate 100 diverse and chemically interesting molecules for computational chemistry evaluation. 

Focus on molecules that would be useful for:
1. Catalyst development (ligands, small organic catalysts)
2. Drug discovery (small drug-like molecules)
3. Materials science (building blocks, functional molecules)
4. Energy storage (electrolytes, organic battery materials)

For each molecule, provide:
- SMILES notation
- Brief description (one line)

Format as a simple list:
1. SMILES | Description
2. SMILES | Description
...

Make sure all SMILES are valid and chemically reasonable. Include variety in:
- Molecular weight (50-400 Da)
- Functional groups (amines, alcohols, carbonyls, aromatics, heterocycles)
- Ring systems (5-membered, 6-membered, fused rings)
- Heteroatoms (N, O, S, P)

Generate the full list of 100 molecules now."""

    start = time.time()
    result = agent.run_agentic_workflow(generation_request)
    elapsed = time.time() - start
    
    print(f"\nGeneration completed in {elapsed:.1f}s")
    print()
    
    # Extract the molecule list from the result
    if result and 'result' in result:
        response = result['result']
        if isinstance(response, dict) and 'message' in response:
            molecule_text = response['message']
        else:
            molecule_text = str(response)
    else:
        molecule_text = str(result)
    
    # Save the generated molecules
    with open('generated_molecules.txt', 'w') as f:
        f.write(molecule_text)
    
    print("="*80)
    print("GENERATED MOLECULES")
    print("="*80)
    print(molecule_text[:2000] + "..." if len(molecule_text) > 2000 else molecule_text)
    print()
    
    # Parse SMILES from the response
    print("="*80)
    print("STEP 2: Parsing SMILES from generated molecules...")
    print("-"*80)
    
    smiles_list = []
    for line in molecule_text.split('\n'):
        if '|' in line and 'SMILES' not in line.upper():
            parts = line.split('|')
            if len(parts) >= 1:
                # Extract SMILES (remove numbering like "1. ")
                smiles_part = parts[0].strip()
                # Remove leading numbers and dots
                import re
                smiles = re.sub(r'^\d+\.\s*', '', smiles_part).strip()
                if smiles and len(smiles) > 3:  # Basic validation
                    smiles_list.append(smiles)
    
    print(f"Parsed {len(smiles_list)} SMILES from response")
    
    if len(smiles_list) < 10:
        print("⚠ Warning: Only found", len(smiles_list), "valid SMILES")
        print("Showing first few lines of response for debugging:")
        print(molecule_text[:500])
        return
    
    # Save parsed SMILES
    with open('parsed_smiles.txt', 'w') as f:
        for i, smi in enumerate(smiles_list, 1):
            f.write(f"{i}. {smi}\n")
    
    print(f"✓ Saved {len(smiles_list)} SMILES to parsed_smiles.txt")
    print()
    
    # Now screen with MACE
    print("="*80)
    print("STEP 3: MACE Rapid Screening")
    print("="*80)
    print(f"Screening {len(smiles_list)} molecules with MACE-MP...")
    print()
    
    # Create MACE screening request
    smiles_json = str(smiles_list[:100])  # Limit to first 100
    
    mace_request = f"""Use MACE rapid screening to evaluate the following {len(smiles_list[:100])} molecules.

SMILES list: {smiles_json}

Use mace_rapid_screening to quickly calculate energies for all molecules. This should complete in under 1 minute for 100 molecules.

After screening, please:
1. Report the top 10 lowest energy molecules
2. Report the top 10 highest energy molecules
3. Provide statistics (mean, std, range)
4. Suggest which 3-5 molecules would be most interesting for detailed DFT validation"""
    
    start = time.time()
    mace_result = agent.run_agentic_workflow(mace_request)
    elapsed = time.time() - start
    
    print(f"\nMACE screening completed in {elapsed:.1f}s")
    print()
    
    # Extract and display results
    if mace_result and 'result' in mace_result:
        mace_response = mace_result['result']
        if isinstance(mace_response, dict) and 'message' in mace_response:
            mace_text = mace_response['message']
        else:
            mace_text = str(mace_response)
    else:
        mace_text = str(mace_result)
    
    # Save MACE results
    with open('mace_screening_results.txt', 'w') as f:
        f.write(mace_text)
    
    print("="*80)
    print("MACE SCREENING RESULTS")
    print("="*80)
    print(mace_text)
    print()
    
    print("="*80)
    print("WORKFLOW COMPLETE")
    print("="*80)
    print(f"✓ Generated {len(smiles_list)} molecules")
    print(f"✓ MACE screening completed in {elapsed:.1f}s")
    print(f"✓ Results saved to mace_screening_results.txt")
    print()
    print("Next steps:")
    print("1. Review top candidates in mace_screening_results.txt")
    print("2. Run DFT validation on most promising molecules")
    print("3. Compare MACE predictions with DFT results")
    print("="*80)

if __name__ == "__main__":
    main()
