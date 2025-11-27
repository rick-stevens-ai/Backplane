#!/usr/bin/env python3
"""
Direct LLM calls to expand NH3 catalyst set (bypass agentic workflow)
"""
import sys
import time
import json
import openai
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

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

def expand_batch_direct(client, existing, batch_num, target_count=10):
    """Direct API call to gpt-oss:120b (no agentic workflow)"""

    # Create summary of existing candidates
    existing_summary = "\n".join([
        f"  {i+1}. {c['name']} ({c['smiles']}) - {c['rationale']}"
        for i, c in enumerate(existing)
    ])

    # Simplified, focused prompt
    request = f"""You are a computational chemistry expert. Based on these {len(existing)} NH3 catalyst candidates:

{existing_summary}

Suggest {target_count} NEW molecules for NH3 formation catalysis. Focus on:
- N-containing ligands (pyridine/imidazole derivatives)
- Proton shuttles and bases
- Metal coordination sites

Format:
SMILES | Name | Brief rationale

Be concise. Provide exactly {target_count} molecules."""

    print(f"\n{'='*80}")
    print(f"BATCH {batch_num}: Requesting {target_count} new candidates")
    print(f"{'='*80}\n")
    print(f"Prompt length: {len(request)} chars")

    start = time.time()

    try:
        # Direct API call with increased timeout
        response = client.chat.completions.create(
            model="gpt-oss:120b",
            messages=[{"role": "user", "content": request}],
            temperature=0.7,
            max_tokens=2000,
            timeout=300.0  # 5 minute timeout
        )

        elapsed = time.time() - start
        text = response.choices[0].message.content

        print(f"✓ Response received in {elapsed:.1f}s")
        return text, elapsed

    except Exception as e:
        elapsed = time.time() - start
        print(f"✗ Error after {elapsed:.1f}s: {type(e).__name__}: {str(e)[:100]}")
        return None, elapsed

def parse_molecules(text):
    """Parse SMILES from response"""
    molecules = []
    for line in text.split('\n'):
        if '|' in line and 'SMILES' not in line.upper():
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 3:
                # Remove numbering
                import re
                smiles = re.sub(r'^\\d+\\.\\s*', '', parts[0]).strip()
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
    print("NH3 CATALYST EXPANSION - DIRECT API APPROACH")
    print("="*80)
    print(f"Starting with {len(existing_candidates)} candidates")
    print(f"Target: 50 total candidates (smaller batches)")
    print("="*80)

    # Initialize OpenAI client directly
    client = openai.OpenAI(
        base_url="http://100.94.58.120:12000/v1",
        api_key="EMPTY"
    )

    all_candidates = existing_candidates.copy()
    batch_num = 1
    target_total = 50  # Reduced target
    batch_size = 10    # Smaller batches

    while len(all_candidates) < target_total and batch_num <= 6:  # Max 6 batches
        remaining = target_total - len(all_candidates)
        batch_count = min(batch_size, remaining)

        # Get next batch
        response_text, elapsed = expand_batch_direct(client, all_candidates, batch_num, batch_count)

        if response_text is None:
            print(f"⚠ Batch {batch_num} failed, stopping")
            break

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
        with open('nh3_candidates_expanded_direct.json', 'w') as f:
            json.dump({
                'total': len(all_candidates),
                'batches_completed': batch_num,
                'candidates': all_candidates
            }, f, indent=2)

        batch_num += 1

        # Small delay between batches
        if len(all_candidates) < target_total:
            print("Waiting 3s before next batch...")
            time.sleep(3)

    print("="*80)
    print("EXPANSION COMPLETE")
    print("="*80)
    print(f"Final count: {len(all_candidates)} candidates")
    print(f"Batches: {batch_num - 1}")
    print()
    print("✓ Saved to: nh3_candidates_expanded_direct.json")
    print("="*80)

if __name__ == "__main__":
    main()
