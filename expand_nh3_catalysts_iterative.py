#!/usr/bin/env python3
"""
Iteratively expand NH3 catalyst set using oss120b logical reasoning
"""
import sys
import time
import json
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent

# Existing candidates from previous work
existing_candidates = [
    {"smiles": "c1ccncc1", "name": "Pyridine", "rationale": "N-heterocycle ligand, Lewis base"},
    {"smiles": "c1cnc[nH]1", "name": "Imidazole", "rationale": "Enzyme-inspired, proton shuttle"},
    {"smiles": "c1ccnc(c1)c2ccccn2", "name": "2,2'-Bipyridine", "rationale": "Bidentate chelator"},
    {"smiles": "N", "name": "Ammonia", "rationale": "Product molecule"},
    {"smiles": "NN", "name": "Hydrazine", "rationale": "N-N bond intermediate"},
    {"smiles": "c1c[nH]nc1", "name": "Pyrazole", "rationale": "Imidazole isomer"},
    {"smiles": "c1ncn[nH]1", "name": "1,2,4-Triazole", "rationale": "Three N atoms for N2 activation"},
]

def expand_batch(agent, existing, batch_num, target_count=15):
    """Ask oss120b to suggest next batch based on existing"""
    
    # Create summary of existing candidates
    existing_summary = "\n".join([
        f"  {i+1}. {c['name']} ({c['smiles']}) - {c['rationale']}" 
        for i, c in enumerate(existing)
    ])
    
    request = f"""Based on the following NH3 formation catalyst candidates we've already identified:

{existing_summary}

Using chemical reasoning, suggest {target_count} NEW molecules that would be promising for NH3 formation catalysis. Consider:

1. **Nitrogen-containing ligands** that could coordinate to metal centers (Ru, Fe, Mo, Re)
2. **Proton shuttles** and bases that facilitate H+ transfer
3. **Molecules with multiple coordination sites** for binding N2, H2, or intermediates 
4. **Systematic variations** of molecules above (substituted pyridines, different heterocycles)
5. **Biomimetic structures** inspired by nitrogenase (FeMo-cofactor has imidazole from histidine)

For each molecule provide:
- SMILES notation (verify it's valid!)
- Name
- Brief rationale (1-2 sentences explaining why it's good for NH3 catalysis)

Format as:
1. SMILES | Name | Rationale
2. SMILES | Name | Rationale
...

Focus on diversity and chemical logic. Avoid duplicating molecules already listed above."""

    print(f"\n{'='*80}")
    print(f"BATCH {batch_num}: Requesting {target_count} new candidates")
    print(f"{'='*80}\n")
    
    start = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed = time.time() - start
    
    print(f"✓ Response received in {elapsed:.1f}s")
    
    # Extract text from result
    if result and 'result' in result:
        response = result['result']
        if isinstance(response, dict) and 'message' in response:
            text = response['message']
        else:
            text = str(response)
    else:
        text = str(result)
    
    return text, elapsed

def parse_molecules(text):
    """Parse SMILES from response"""
    molecules = []
    for line in text.split('\n'):
        if '|' in line and 'SMILES' not in line.upper():
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 3:
                # Remove numbering
                import re
                smiles = re.sub(r'^\d+\.\s*', '', parts[0]).strip()
                name = parts[1].strip()
                rationale = parts[2].strip()
                
                if smiles and len(smiles) > 2:
                    molecules.append({
                        'smiles': smiles,
                        'name': name,
                        'rationale': rationale
                    })
    return molecules

def main():
    print("="*80)
    print("NH3 CATALYST EXPANSION - ITERATIVE APPROACH")
    print("="*80)
    print(f"Starting with {len(existing_candidates)} candidates")
    print(f"Target: 100 total candidates")
    print("="*80)
    
    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-01"  # Using 20B model (faster)
    )
    
    all_candidates = existing_candidates.copy()
    batch_num = 1
    target_total = 100
    batch_size = 15
    
    while len(all_candidates) < target_total and batch_num <= 10:  # Max 10 batches
        remaining = target_total - len(all_candidates)
        batch_count = min(batch_size, remaining)
        
        # Get next batch
        response_text, elapsed = expand_batch(agent, all_candidates, batch_num, batch_count)
        
        # Save raw response
        with open(f'batch_{batch_num}_response.txt', 'w') as f:
            f.write(response_text)
        
        # Parse molecules
        new_molecules = parse_molecules(response_text)
        print(f"✓ Parsed {len(new_molecules)} molecules from batch {batch_num}")
        
        if not new_molecules:
            print("⚠ No molecules parsed, stopping")
            break
        
        # Add to collection
        all_candidates.extend(new_molecules)
        print(f"✓ Total candidates: {len(all_candidates)}/{target_total}\n")
        
        # Save progress
        with open('nh3_candidates_expanded.json', 'w') as f:
            json.dump({
                'total': len(all_candidates),
                'batches_completed': batch_num,
                'candidates': all_candidates
            }, f, indent=2)
        
        batch_num += 1
        
        # Small delay between batches
        if len(all_candidates) < target_total:
            time.sleep(2)
    
    print("="*80)
    print("EXPANSION COMPLETE")
    print("="*80)
    print(f"Final count: {len(all_candidates)} candidates")
    print(f"Batches: {batch_num - 1}")
    print()
    print("✓ Saved to: nh3_candidates_expanded.json")
    print()
    print("Next steps:")
    print("1. Screen all candidates with MACE")
    print("2. Validate top 10 with DFT")
    print("="*80)

if __name__ == "__main__":
    main()
