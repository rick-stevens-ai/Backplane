#!/usr/bin/env python3
"""
Query gpt-oss:120b for 3 additional NH3 catalyst candidates
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
import json

def main():
    """Query agent for additional catalyst suggestions"""

    print("="*80)
    print("QUERYING GPT-OSS:120B FOR ADDITIONAL NH3 CATALYST CANDIDATES")
    print("="*80)
    print()

    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    request = """We are screening potential catalysts for NH3 synthesis from N2 and H2.

We have already tested these 5 molecules:
1. Pyridine (c1ccncc1) - Simple N-heterocycle, Lewis base
2. Imidazole (c1cnc[nH]1) - Enzyme-inspired, proton shuttle
3. 2,2'-Bipyridine (c1ccnc(c1)c2ccccn2) - Bidentate chelator
4. Ammonia (N) - Product molecule baseline
5. Hydrazine (NN) - N-N bond intermediate

Please suggest 3 ADDITIONAL diverse molecular catalysts or ligands that could be interesting
for NH3 synthesis. For each molecule provide:
- Name
- SMILES string (must be valid)
- Chemical formula
- Expected number of atoms (including hydrogens)
- Brief rationale (how it could help NH3 synthesis)

Focus on molecules that:
- Have nitrogen coordination sites
- Are small-to-medium sized (< 30 atoms for fast DFT)
- Are chemically diverse from the 5 already tested
- Could coordinate to metal centers or activate N2/H2

Format your response as:
1. Name | SMILES | Formula | Atoms | Rationale
2. Name | SMILES | Formula | Atoms | Rationale
3. Name | SMILES | Formula | Atoms | Rationale

Be concise and use standard SMILES notation."""

    print("Request sent to agent...\n")

    result = agent.run_agentic_workflow(request)

    if result and result.get('status') == 'SUCCESS':
        response = result.get('result', {}).get('response', '')

        print("="*80)
        print("AGENT RESPONSE:")
        print("="*80)
        print(response)
        print("="*80)

        # Save response
        with open('additional_catalysts_suggestions.txt', 'w') as f:
            f.write(response)

        print("\n✓ Response saved to: additional_catalysts_suggestions.txt")
        return True
    else:
        print("✗ Failed to get response from agent")
        print(f"Status: {result.get('status') if result else 'NO_RESULT'}")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
