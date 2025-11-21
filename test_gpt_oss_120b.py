#!/usr/bin/env python3
"""
Test the agentic workflow with gpt-oss:120b model on containerized server
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent import SimulationAgent

def main():
    """Test agent with gpt-oss:120b model"""

    print("\n" + "="*80)
    print("TESTING AGENT WITH GPT-OSS:120B MODEL (CONTAINERIZED SERVER)")
    print("="*80)

    try:
        agent = SimulationAgent(
            server_config_path="spark_servers.yaml",
            server_name="spark-container-03"
        )
        print("✓ Agent initialized successfully")
        print(f"  Model: {agent.model}")
        print(f"  API Base: {agent.server_config['openai_api_base']}")
        print(f"  Server: {agent.server_config['server']}")

        # Test with complex protein simulation request
        user_request = """Please run a DFT optimization simulation on an ethanol molecule (CCO)
        for an experiment called 'Ethanol Molecular Structure Analysis'.
        Set the temperature parameter to 298.15 K and analyze the molecular properties."""

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
            print(f"Dipole Moment: {result['result']['results']['dipole_moment_debye']} Debye")

        return 0

    except Exception as e:
        print(f"\n✗ TEST FAILED")
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
