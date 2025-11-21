#!/usr/bin/env python3
"""
Test the agentic workflow with local Ollama server
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent import SimulationAgent

def main():
    """Test agent with local Ollama models"""

    print("\n" + "="*80)
    print("TESTING AGENT WITH LOCAL OLLAMA SERVER")
    print("="*80)

    # Test with gpt-oss:20b model
    print("\nInitializing agent with gpt-oss:20b model...")
    try:
        agent = SimulationAgent(
            server_config_path="spark_servers_local.yaml",
            server_name="local-gpt-oss-20b"
        )
        print("✓ Agent initialized successfully")
        print(f"  Model: {agent.model}")
        print(f"  API Base: {agent.server_config['openai_api_base']}")

        # Test user request
        user_request = """Please run a DFT optimization simulation on a water molecule (H2O)
        for an experiment called 'Water Molecule Analysis'. Set the temperature parameter to 298.15 K."""

        print(f"\n{'-'*80}")
        print("Running agentic workflow...")
        print(f"{'-'*80}")

        result = agent.run_agentic_workflow(user_request)

        print("\n" + "="*80)
        print("✓ TEST COMPLETED SUCCESSFULLY")
        print("="*80)
        print(f"\nFinal Job Status: {result['status']}")
        if result.get('result'):
            print(f"Energy Level: {result['result']['results']['energy_level_hartree']} Hartree")
            print(f"Convergence: {result['result']['results']['convergence_achieved']}")
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
