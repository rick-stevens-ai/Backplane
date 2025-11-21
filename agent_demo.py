#!/usr/bin/env python3
"""
Demo Agentic Workflow for Scientific Simulation
Simulates an AI agent submitting and monitoring simulations WITHOUT requiring external LLM
"""
import time
import json
import requests
from typing import Dict, Any


class SimulationAgentDemo:
    """Demo agent that submits and monitors scientific simulations"""

    def __init__(self):
        """Initialize the demo agent"""
        self.api_base = "http://127.0.0.1:8000"  # FastAPI server

    def submit_job(self, experiment_data: Dict[str, Any]) -> str:
        """Submit a simulation job to the FastAPI server"""
        response = requests.post(
            f"{self.api_base}/submit_job",
            json=experiment_data
        )
        response.raise_for_status()
        result = response.json()
        return result['job_id']

    def check_status(self, job_id: str) -> Dict[str, Any]:
        """Check the status of a submitted job"""
        response = requests.get(f"{self.api_base}/job_status/{job_id}")
        response.raise_for_status()
        return response.json()

    def poll_until_complete(self, job_id: str, poll_interval: int = 2, max_wait: int = 300) -> Dict[str, Any]:
        """Poll job status until completion or timeout"""
        start_time = time.time()

        while (time.time() - start_time) < max_wait:
            status = self.check_status(job_id)

            if status['status'] == 'SUCCESS':
                return status
            elif status['status'] == 'FAILURE':
                raise Exception(f"Job {job_id} failed: {status}")

            elapsed = time.time() - start_time
            print(f"  [{elapsed:.1f}s] Job {job_id[:8]}... status: {status['status']}")
            time.sleep(poll_interval)

        raise TimeoutError(f"Job {job_id} did not complete within {max_wait} seconds")

    def simulate_agent_thinking(self, step: str, content: str):
        """Simulate agent's thinking process"""
        print(f"\n{'='*80}")
        print(f"[AGENT REASONING] {step}")
        print(f"{'='*80}")
        print(content)
        print(f"{'='*80}\n")
        time.sleep(0.5)  # Dramatic pause

    def run_agentic_workflow_demo(self, user_request: str) -> Dict[str, Any]:
        """
        Demo agentic workflow that simulates AI agent behavior:
        1. Parse user request
        2. Decide to submit simulation
        3. Submit job
        4. Poll until complete
        5. Interpret and return results
        """
        print(f"\n{'='*80}")
        print(f"AGENTIC WORKFLOW DEMO - SIMULATED AI AGENT")
        print(f"{'='*80}")
        print(f"User Request: {user_request}\n")

        # Step 1: Agent "understands" the request
        self.simulate_agent_thinking(
            "Step 1: Understanding User Request",
            "I received a request to run a DFT optimization simulation.\n"
            "Parsing the request...\n"
            "- Molecule: Li[C6] (Lithium-ion battery compound)\n"
            "- Simulation type: DFT optimization\n"
            "- Experiment name: Battery Performance Analysis\n"
            "- Temperature: 298.15 K\n\n"
            "I need to use the submit_simulation tool to submit this job."
        )

        # Step 2: Agent decides to submit
        self.simulate_agent_thinking(
            "Step 2: Preparing Simulation Parameters",
            "Constructing simulation request with parameters:\n"
            "{\n"
            '  "experiment_name": "Battery Performance Analysis",\n'
            '  "molecule_smiles": "Li[C6]",\n'
            '  "simulation_type": "dft_optimization",\n'
            '  "parameters": {"temperature": 298.15}\n'
            "}\n\n"
            "Calling submit_simulation function..."
        )

        # Step 3: Submit the job
        experiment_data = {
            "experiment_name": "Battery Performance Analysis",
            "molecule_smiles": "Li[C6]",
            "simulation_type": "dft_optimization",
            "parameters": {"temperature": 298.15}
        }

        job_id = self.submit_job(experiment_data)
        print(f"[ACTION] Job submitted successfully!")
        print(f"         Job ID: {job_id}\n")

        # Step 4: Agent monitors the job
        self.simulate_agent_thinking(
            "Step 3: Monitoring Simulation Progress",
            f"Job submitted with ID: {job_id}\n"
            "The simulation is now running on the HPC cluster.\n"
            "I will poll the status every 2 seconds until completion.\n\n"
            "Using check_simulation_status function to monitor progress..."
        )

        # Poll for completion
        print(f"[MONITORING] Polling job status...")
        result = self.poll_until_complete(job_id)

        # Step 5: Agent interprets results
        self.simulate_agent_thinking(
            "Step 4: Interpreting Results",
            "Simulation completed successfully!\n\n"
            "Analyzing the results...\n"
            f"- Energy Level: {result['result']['results']['energy_level_hartree']} Hartree\n"
            f"- Convergence: {result['result']['results']['convergence_achieved']}\n"
            f"- HOMO-LUMO Gap: {result['result']['results']['homo_lumo_gap_ev']} eV\n"
            f"- Dipole Moment: {result['result']['results']['dipole_moment_debye']} Debye\n"
            f"- Computation Time: {result['result']['results']['computation_time_seconds']} seconds\n\n"
            "Preparing human-readable summary for the user..."
        )

        # Step 6: Generate final response
        print(f"\n{'='*80}")
        print(f"[AGENT FINAL RESPONSE]")
        print(f"{'='*80}")
        print(f"""
I've successfully completed the DFT optimization simulation for your lithium-ion
battery molecule (Li[C6]) under the experiment "Battery Performance Analysis".

Here are the key findings:

**Energy Analysis:**
- Total Energy: {result['result']['results']['energy_level_hartree']} Hartree ({result['result']['results']['energy_level_ev']} eV)
- The calculation converged successfully after {result['result']['results']['optimization_steps']} optimization steps

**Electronic Properties:**
- HOMO-LUMO Gap: {result['result']['results']['homo_lumo_gap_ev']} eV
  (This indicates good electronic conductivity for battery applications)
- Dipole Moment: {result['result']['results']['dipole_moment_debye']} Debye

**Structural Analysis:**
- Final geometry RMSD: {result['result']['results']['final_geometry_rmsd']} Å
- The structure is well-optimized with minimal structural deviation

**Output Files:**
All simulation data has been saved to:
- Checkpoint: {result['result']['output_files']['checkpoint']}
- Log file: {result['result']['output_files']['log_file']}
- Geometry: {result['result']['output_files']['geometry']}

The simulation completed in {result['result']['results']['computation_time_seconds']} seconds at the specified
temperature of {experiment_data['parameters']['temperature']} K.
""")
        print(f"{'='*80}\n")

        return result


def main():
    """Main entry point for the demo agentic workflow"""

    print("="*80)
    print("SCIENTIFIC AI AGENT DEMO")
    print("="*80)
    print("This demo simulates an AI agent that:")
    print("  1. Understands natural language requests")
    print("  2. Submits scientific simulations")
    print("  3. Monitors progress")
    print("  4. Interprets and explains results")
    print("="*80)

    # Initialize demo agent
    agent = SimulationAgentDemo()

    # Example user request
    user_request = """Please run a DFT optimization simulation on a lithium-ion battery
molecule (Li[C6]) for an experiment called 'Battery Performance Analysis'.
Set the temperature parameter to 298.15 K."""

    # Run the agentic workflow
    try:
        result = agent.run_agentic_workflow_demo(user_request)

        print("\n" + "="*80)
        print("RAW API RESPONSE (for debugging)")
        print("="*80)
        print(json.dumps(result, indent=2))
        print("="*80)

        print("\n✓ Demo completed successfully!")
        print("\nThis demonstrates how an AI agent can:")
        print("  • Parse natural language requests")
        print("  • Make decisions about which tools to use")
        print("  • Submit long-running jobs asynchronously")
        print("  • Monitor job progress")
        print("  • Interpret scientific results for users")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
