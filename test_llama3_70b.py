#!/usr/bin/env python3
"""
Test the agentic workflow with llama3:70b model (largest available)
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent import SimulationAgent

def main():
    """Test agent with llama3:70b model"""

    print("\n" + "="*80)
    print("TESTING AGENT WITH LLAMA3:70B MODEL")
    print("="*80)

    try:
        agent = SimulationAgent(
            server_config_path="spark_servers_local.yaml",
            server_name="local-llama3-70b"
        )
        print("✓ Agent initialized successfully")
        print(f"  Model: {agent.model}")
        print(f"  API Base: {agent.server_config['openai_api_base']}")

        # Test with protein simulation request
        user_request = """Please run a molecular dynamics simulation on an aspirin molecule (CC(=O)Oc1ccccc1C(=O)O)
        for an experiment called 'Aspirin Degradation Study'. Use temperature 300 K and study the stability."""

        print(f"\n{'-'*80}")
        print("Running agentic workflow...")
        print(f"{'-'*80}")

        result = agent.run_agentic_workflow(user_request)

        print("\n" + "="*80)
        print("✓ TEST COMPLETED SUCCESSFULLY")
        print("="*80)
        print(f"\nFinal Job Status: {result['status']}")
        if result.get('result'):
            print(f"Experiment: {result['result']['experiment_name']}")
            print(f"Molecule: {result['result']['molecule_smiles']}")
            print(f"Simulation Type: {result['result']['simulation_type']}")
            print(f"Energy Level: {result['result']['results']['energy_level_hartree']} Hartree")
            print(f"HOMO-LUMO Gap: {result['result']['results']['homo_lumo_gap_ev']} eV")

        return 0

    except Exception as e:
        print(f"\n✗ TEST FAILED")
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
