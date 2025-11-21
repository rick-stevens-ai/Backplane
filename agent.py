#!/usr/bin/env python3
"""
Agentic Workflow for Scientific Simulation
Uses OpenAI-compatible API from spark_servers.yaml to submit and monitor simulations
"""
import time
import json
import yaml
import requests
from typing import Dict, Any, Optional
from openai import OpenAI


class SimulationAgent:
    """Agent that submits and monitors scientific simulations"""

    def __init__(self, server_config_path: str = "spark_servers.yaml", server_name: str = "spark-36ae"):
        """
        Initialize the agent with a specific server configuration

        Args:
            server_config_path: Path to the YAML configuration file
            server_name: Name of the server to use from the config
        """
        self.api_base = "http://127.0.0.1:8000"  # FastAPI server
        self.server_config = self._load_server_config(server_config_path, server_name)
        self.client = OpenAI(
            api_key=self.server_config["openai_api_key"],
            base_url=self.server_config["openai_api_base"]
        )
        self.model = self.server_config["openai_model"]

    def _load_server_config(self, config_path: str, server_name: str) -> Dict[str, str]:
        """Load server configuration from YAML file"""
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        for server in config['servers']:
            if server['server'] == server_name:
                return server

        raise ValueError(f"Server {server_name} not found in {config_path}")

    def submit_job(self, experiment_data: Dict[str, Any]) -> str:
        """
        Submit a simulation job to the FastAPI server

        Args:
            experiment_data: Dictionary with experiment parameters

        Returns:
            job_id: Unique identifier for the submitted job
        """
        response = requests.post(
            f"{self.api_base}/submit_job",
            json=experiment_data
        )
        response.raise_for_status()
        result = response.json()
        return result['job_id']

    def check_status(self, job_id: str) -> Dict[str, Any]:
        """
        Check the status of a submitted job

        Args:
            job_id: Unique identifier for the job

        Returns:
            Dictionary with job status and results (if complete)
        """
        response = requests.get(f"{self.api_base}/job_status/{job_id}")
        response.raise_for_status()
        return response.json()

    def poll_until_complete(self, job_id: str, poll_interval: int = 2, max_wait: int = 300) -> Dict[str, Any]:
        """
        Poll job status until completion or timeout

        Args:
            job_id: Unique identifier for the job
            poll_interval: Seconds between status checks
            max_wait: Maximum seconds to wait

        Returns:
            Final job result dictionary
        """
        start_time = time.time()

        while (time.time() - start_time) < max_wait:
            status = self.check_status(job_id)

            if status['status'] == 'SUCCESS':
                return status
            elif status['status'] == 'FAILURE':
                raise Exception(f"Job {job_id} failed: {status}")

            print(f"[{time.time() - start_time:.1f}s] Job {job_id} status: {status['status']}")
            time.sleep(poll_interval)

        raise TimeoutError(f"Job {job_id} did not complete within {max_wait} seconds")

    def run_agentic_workflow(self, user_request: str) -> Dict[str, Any]:
        """
        Run the full agentic workflow:
        1. Agent receives user request
        2. Agent decides to submit a simulation
        3. Agent submits the job
        4. Agent polls until completion
        5. Agent returns the results

        Args:
            user_request: Natural language request from user

        Returns:
            Dictionary with simulation results
        """
        print(f"\n{'='*80}")
        print(f"AGENTIC WORKFLOW STARTED")
        print(f"{'='*80}")
        print(f"User Request: {user_request}\n")

        # Define tools/functions the agent can use
        tools = [
            {
                "type": "function",
                "function": {
                    "name": "submit_simulation",
                    "description": "Submit a scientific simulation job for molecular analysis. Use this when you need to run DFT optimization, molecular dynamics, or docking simulations.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "experiment_name": {
                                "type": "string",
                                "description": "Descriptive name for the experiment"
                            },
                            "molecule_smiles": {
                                "type": "string",
                                "description": "SMILES string representation of the molecule to analyze"
                            },
                            "simulation_type": {
                                "type": "string",
                                "enum": ["dft_optimization", "md", "docking"],
                                "description": "Type of simulation to run"
                            },
                            "parameters": {
                                "type": "object",
                                "description": "Optional physical parameters (e.g., temperature, pressure)",
                                "additionalProperties": {"type": "number"}
                            }
                        },
                        "required": ["experiment_name", "molecule_smiles", "simulation_type"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "check_simulation_status",
                    "description": "Check the status of a submitted simulation job. Returns the current status and results if complete.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "job_id": {
                                "type": "string",
                                "description": "The unique job ID returned from submit_simulation"
                            }
                        },
                        "required": ["job_id"]
                    }
                }
            }
        ]

        # Initial prompt to the agent
        messages = [
            {
                "role": "system",
                "content": """You are a scientific computing assistant. You help users submit and monitor molecular simulations.

When a user asks you to run a simulation:
1. Use submit_simulation to submit the job
2. Use check_simulation_status to monitor progress
3. Keep checking until status is SUCCESS
4. Report the final results to the user

You have access to these tools:
- submit_simulation: Submit a new simulation job
- check_simulation_status: Check if a job is complete"""
            },
            {
                "role": "user",
                "content": user_request
            }
        ]

        max_iterations = 20
        iteration = 0
        job_id = None

        while iteration < max_iterations:
            iteration += 1
            print(f"\n--- Iteration {iteration} ---")

            # Get agent's response
            response = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                tools=tools,
                tool_choice="auto"
            )

            message = response.choices[0].message
            messages.append(message)

            # Check if agent wants to use tools
            if message.tool_calls:
                print(f"Agent is calling {len(message.tool_calls)} tool(s)...")

                for tool_call in message.tool_calls:
                    function_name = tool_call.function.name
                    arguments = json.loads(tool_call.function.arguments)

                    print(f"\nTool: {function_name}")
                    print(f"Arguments: {json.dumps(arguments, indent=2)}")

                    # Execute the requested function
                    if function_name == "submit_simulation":
                        job_id = self.submit_job(arguments)
                        result = {"job_id": job_id, "status": "QUEUED"}
                        print(f"Result: Job submitted with ID: {job_id}")

                    elif function_name == "check_simulation_status":
                        result = self.check_status(arguments["job_id"])
                        print(f"Result: Status = {result['status']}")
                        if result.get('result'):
                            print(f"Simulation Results: {json.dumps(result['result'], indent=2)}")

                    else:
                        result = {"error": f"Unknown function: {function_name}"}

                    # Add function result to messages
                    messages.append({
                        "role": "tool",
                        "tool_call_id": tool_call.id,
                        "name": function_name,
                        "content": json.dumps(result)
                    })

            else:
                # Agent provided a final response
                print(f"\n{'='*80}")
                print(f"AGENT RESPONSE:")
                print(f"{'='*80}")
                print(message.content)
                print(f"{'='*80}\n")

                # If we have a job_id and agent is done, return the final result
                if job_id:
                    final_status = self.check_status(job_id)
                    return final_status

                return {"message": message.content}

        raise Exception("Agent exceeded maximum iterations without completing task")


def main():
    """Main entry point for the agentic workflow"""

    # Initialize agent with qwen3-coder model
    print("Initializing Scientific Simulation Agent...")
    agent = SimulationAgent(server_name="spark-36ae")  # Using qwen3-coder

    # Example user request
    user_request = """Please run a DFT optimization simulation on a lithium-ion battery molecule (Li[C6])
    for an experiment called 'Battery Performance Analysis'. Set the temperature parameter to 298.15 K."""

    # Run the agentic workflow
    try:
        result = agent.run_agentic_workflow(user_request)

        print("\n" + "="*80)
        print("FINAL RESULTS")
        print("="*80)
        print(json.dumps(result, indent=2))
        print("="*80)

    except Exception as e:
        print(f"\nError: {e}")
        raise


if __name__ == "__main__":
    main()
