#!/usr/bin/env python3
"""
Extended Agentic Workflow for Computational Chemistry Applications
Supports Quantum ESPRESSO, CP2K, GPAW, LAMMPS, and GROMACS
+ MACE ML foundation model for rapid molecular screening
"""
import time
import json
import yaml
import requests
from typing import Dict, Any, Optional, List
from openai import OpenAI
from mace_client_simple import MACEClient, smiles_to_xyz


class ComputationalChemistryAgent:
    """Agent that submits and monitors computational chemistry simulations across multiple applications"""

    def __init__(self, server_config_path: str = "spark_servers.yaml", server_name: str = "spark-container-03"):
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

        # Initialize MACE client for rapid ML predictions
        self.mace_client = MACEClient()

    def _load_server_config(self, config_path: str, server_name: str) -> Dict[str, str]:
        """Load server configuration from YAML file"""
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        for server in config['servers']:
            if server['server'] == server_name:
                return server

        raise ValueError(f"Server {server_name} not found in {config_path}")

    def submit_app_job(self, application: str, job_params: Dict[str, Any], experiment_name: str = "simulation") -> str:
        """
        Submit a simulation job for a specific application

        Args:
            application: Application name ('quantum_espresso', 'cp2k', 'gpaw', 'lammps', 'gromacs')
            job_params: Application-specific parameters
            experiment_name: Name of the experiment

        Returns:
            job_id: Unique identifier for the submitted job
        """
        job_data = {
            "application": application,
            "job_params": job_params,
            "experiment_name": experiment_name
        }

        # Use the new run_simulation endpoint
        response = requests.post(
            f"{self.api_base}/submit_app_job",
            json=job_data
        )
        response.raise_for_status()
        result = response.json()
        return result['job_id']

    def check_status(self, job_id: str) -> Dict[str, Any]:
        """Check the status of a submitted job"""
        response = requests.get(f"{self.api_base}/job_status/{job_id}")
        response.raise_for_status()
        return response.json()

    def run_agentic_workflow(self, user_request: str) -> Dict[str, Any]:
        """
        Run the full agentic workflow with support for multiple computational chemistry applications

        Args:
            user_request: Natural language request from user

        Returns:
            Dictionary with simulation results
        """
        print(f"\n{'='*80}")
        print(f"COMPUTATIONAL CHEMISTRY AGENTIC WORKFLOW")
        print(f"{'='*80}")
        print(f"User Request: {user_request}\n")

        # Define tools for all five applications
        tools = [
            {
                "type": "function",
                "function": {
                    "name": "run_quantum_espresso",
                    "description": "Run a Quantum ESPRESSO plane-wave DFT calculation. Use for electronic structure, band structure, or geometry optimization of molecules and materials.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "experiment_name": {
                                "type": "string",
                                "description": "Descriptive name for the experiment"
                            },
                            "calculation": {
                                "type": "string",
                                "enum": ["scf", "relax", "vc-relax", "md"],
                                "description": "Type of calculation: scf (energy), relax (geometry optimization), vc-relax (cell optimization), md (molecular dynamics)"
                            },
                            "molecule_smiles": {
                                "type": "string",
                                "description": "SMILES string of the molecule to simulate"
                            },
                            "cutoff_wfc": {
                                "type": "number",
                                "description": "Wavefunction cutoff in Rydberg (default: 50)"
                            },
                            "cutoff_rho": {
                                "type": "number",
                                "description": "Charge density cutoff in Rydberg (default: 400)"
                            }
                        },
                        "required": ["experiment_name", "calculation", "molecule_smiles"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "run_cp2k",
                    "description": "Run a CP2K mixed Gaussian/plane-wave DFT calculation. Excellent for QM/MM simulations, large systems, and metal clusters. Supports SMILES (organic molecules) or XYZ coordinates (metals, materials, any atomic system).",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "experiment_name": {
                                "type": "string",
                                "description": "Descriptive name for the experiment"
                            },
                            "run_type": {
                                "type": "string",
                                "enum": ["ENERGY", "GEO_OPT", "MD", "CELL_OPT"],
                                "description": "Type of calculation: ENERGY (single point), GEO_OPT (optimization), MD (molecular dynamics), CELL_OPT (cell optimization)"
                            },
                            "molecule_smiles": {
                                "type": "string",
                                "description": "SMILES string of the molecule to simulate (for organic molecules)"
                            },
                            "xyz_structure": {
                                "type": "string",
                                "description": "XYZ format coordinates for metal clusters, materials, or any atomic system. Format: <n_atoms>\\n<comment>\\n<element> <x> <y> <z>\\n..."
                            },
                            "functional": {
                                "type": "string",
                                "enum": ["PBE", "BLYP", "B3LYP"],
                                "description": "Exchange-correlation functional (default: PBE)"
                            },
                            "cutoff": {
                                "type": "number",
                                "description": "Plane wave cutoff in Rydberg (default: 400)"
                            }
                        },
                        "required": ["experiment_name", "run_type"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "run_gpaw",
                    "description": "Run a GPAW real-space grid DFT calculation in Python. Best for rapid prototyping and Python integration. Supports SMILES (organic molecules) or XYZ coordinates (metals, materials, any atomic system).",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "experiment_name": {
                                "type": "string",
                                "description": "Descriptive name for the experiment"
                            },
                            "calculation_type": {
                                "type": "string",
                                "enum": ["energy", "relaxation", "md"],
                                "description": "Type of calculation: energy (single point), relaxation (geometry optimization), md (molecular dynamics)"
                            },
                            "molecule_smiles": {
                                "type": "string",
                                "description": "SMILES string of the molecule to simulate (for organic molecules)"
                            },
                            "xyz_structure": {
                                "type": "string",
                                "description": "XYZ format coordinates for metal clusters, materials, or any atomic system. Format: <n_atoms>\\n<comment>\\n<element> <x> <y> <z>\\n..."
                            },
                            "mode": {
                                "type": "string",
                                "enum": ["fd", "pw", "lcao"],
                                "description": "GPAW mode: fd (finite difference), pw (plane wave), lcao (linear combination of atomic orbitals). Default: fd"
                            },
                            "xc": {
                                "type": "string",
                                "description": "Exchange-correlation functional (default: PBE)"
                            },
                            "h": {
                                "type": "number",
                                "description": "Grid spacing in Angstrom (default: 0.2)"
                            }
                        },
                        "required": ["experiment_name", "calculation_type"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "run_lammps",
                    "description": "Run a LAMMPS classical molecular dynamics simulation. Supports ReaxFF reactive force fields for catalysis studies.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "experiment_name": {
                                "type": "string",
                                "description": "Descriptive name for the experiment"
                            },
                            "simulation_type": {
                                "type": "string",
                                "enum": ["energy", "minimize", "md", "reaxff"],
                                "description": "Type of simulation: energy (single point), minimize (energy minimization), md (molecular dynamics), reaxff (reactive MD)"
                            },
                            "molecule_smiles": {
                                "type": "string",
                                "description": "SMILES string of the molecule to simulate"
                            },
                            "force_field": {
                                "type": "string",
                                "enum": ["lj", "reaxff", "eam"],
                                "description": "Force field type (default: lj)"
                            },
                            "temperature": {
                                "type": "number",
                                "description": "Temperature in Kelvin (default: 300)"
                            },
                            "timestep": {
                                "type": "number",
                                "description": "Timestep in femtoseconds (default: 1.0)"
                            },
                            "steps": {
                                "type": "integer",
                                "description": "Number of MD steps (default: 1000)"
                            }
                        },
                        "required": ["experiment_name", "simulation_type", "molecule_smiles"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "run_gromacs",
                    "description": "Run a GROMACS molecular dynamics simulation. Excellent for biomolecular systems and enzyme catalysis.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "experiment_name": {
                                "type": "string",
                                "description": "Descriptive name for the experiment"
                            },
                            "simulation_type": {
                                "type": "string",
                                "enum": ["em", "nvt", "npt", "md"],
                                "description": "Type of simulation: em (energy minimization), nvt (constant volume), npt (constant pressure), md (production MD)"
                            },
                            "molecule_smiles": {
                                "type": "string",
                                "description": "SMILES string of the molecule to simulate"
                            },
                            "temperature": {
                                "type": "number",
                                "description": "Temperature in Kelvin (default: 300)"
                            },
                            "pressure": {
                                "type": "number",
                                "description": "Pressure in bar for NPT (default: 1.0)"
                            },
                            "timestep": {
                                "type": "number",
                                "description": "Timestep in picoseconds (default: 0.002)"
                            },
                            "steps": {
                                "type": "integer",
                                "description": "Number of MD steps (default: 50000)"
                            },
                            "force_field": {
                                "type": "string",
                                "description": "Force field name (default: oplsaa)"
                            }
                        },
                        "required": ["experiment_name", "simulation_type", "molecule_smiles"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "check_simulation_status",
                    "description": "Check the status of any submitted simulation job. Returns the current status and results if complete.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "job_id": {
                                "type": "string",
                                "description": "The unique job ID returned from any run_* function"
                            }
                        },
                        "required": ["job_id"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "mace_predict_energy",
                    "description": "Fast ML energy prediction for a single molecule (~0.3s). Use for quick energy estimates or feasibility checks before running expensive DFT calculations. Provides near-quantum accuracy (98% of DFT).",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "smiles": {
                                "type": "string",
                                "description": "SMILES string of the molecule"
                            },
                            "model_type": {
                                "type": "string",
                                "enum": ["off", "mp"],
                                "description": "Model type: 'off' for organic molecules (H,C,N,O,P,S,F,Cl,Br,I), 'mp' for materials/metals. Default: off"
                            },
                            "model_size": {
                                "type": "string",
                                "enum": ["small", "medium", "large"],
                                "description": "Model size: small (fastest), medium (balanced), large (most accurate). Default: medium"
                            }
                        },
                        "required": ["smiles"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "mace_rapid_screening",
                    "description": "Fast ML screening of multiple molecules (~0.3s per molecule). Returns ranked list sorted by predicted energy. Use for large-scale screening (10-100+ molecules) before DFT validation of top candidates. Enables 10-100x speedup compared to all-DFT workflow.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "smiles_list": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "List of SMILES strings to screen"
                            },
                            "model_type": {
                                "type": "string",
                                "enum": ["off", "mp"],
                                "description": "Model type: 'off' for organic, 'mp' for materials. Default: off"
                            },
                            "top_n": {
                                "type": "integer",
                                "description": "Return only top N candidates by energy. If not specified, returns all."
                            }
                        },
                        "required": ["smiles_list"]
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "mace_optimize_geometry",
                    "description": "ML-based geometry optimization (~2s). Faster than DFT optimization but slightly less accurate. Use to prepare good starting structures for DFT calculations or for quick structure refinement.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "smiles": {
                                "type": "string",
                                "description": "SMILES string of the molecule to optimize"
                            },
                            "model_type": {
                                "type": "string",
                                "enum": ["off", "mp"],
                                "description": "Model type: 'off' for organic, 'mp' for materials. Default: off"
                            },
                            "fmax": {
                                "type": "number",
                                "description": "Force convergence criterion in eV/Å (default: 0.05)"
                            },
                            "steps": {
                                "type": "integer",
                                "description": "Maximum optimization steps (default: 200)"
                            }
                        },
                        "required": ["smiles"]
                    }
                }
            }
        ]

        # Initial prompt to the agent
        messages = [
            {
                "role": "system",
                "content": """You are an expert computational chemistry assistant with access to five DFT/MD simulation codes plus MACE ML foundation model:

**High-Accuracy DFT Simulations** (2-5 min per molecule):
1. **Quantum ESPRESSO (QE)** - Plane-wave DFT for electronic structure calculations
2. **CP2K** - Mixed Gaussian/plane-wave DFT, excellent for QM/MM and metal clusters
3. **GPAW** - Real-space grid DFT in Python, great for rapid prototyping
4. **LAMMPS** - Classical MD with ReaxFF for large-scale simulations and catalysis
5. **GROMACS** - Biomolecular MD for enzyme catalysis and protein studies

**Fast ML Predictions** (~0.3s per molecule):
6. **MACE** - Machine learning foundation model for rapid molecular screening
   - mace_predict_energy: Single molecule energy (~0.3s)
   - mace_rapid_screening: Batch screening of 10-100+ molecules (~30s total)
   - mace_optimize_geometry: ML geometry optimization (~2s)

**IMPORTANT - Molecular Input Formats**:
- **SMILES strings**: Use for organic molecules (H, C, N, O, P, S, halogens)
  - Example: "CCO" for ethanol, "N" for ammonia
- **XYZ coordinates**: Use for metal clusters, materials, surfaces, coordination complexes, or any system that cannot be represented as SMILES
  - Format: <n_atoms>\\n<comment>\\n<element> <x> <y> <z>\\n...
  - Example: Ru10 cluster, Fe3O4, Pt surfaces, FeMo-cofactor
  - CP2K and GPAW support XYZ input via the `xyz_structure` parameter
  - When user provides XYZ coordinates in their request, extract them and pass via `xyz_structure` parameter

When a user requests simulations:
1. Choose the most appropriate application based on the scientific question
2. **Determine input format**: If user provides XYZ coordinates or discusses metal clusters/materials, use `xyz_structure`; otherwise use `molecule_smiles`
3. **For screening many molecules (>10)**: Use MACE first to rank candidates, then validate top candidates with DFT
4. **For single molecules**: Use MACE for quick estimates, DFT for final accurate results
5. Set reasonable parameters for the calculation
6. Monitor DFT jobs with check_simulation_status
7. Report results with scientific interpretation

Guidelines for choosing applications:
- Use **MACE** for rapid screening, pre-filtering, or quick feasibility checks
- Use **CP2K** for metal clusters, coordination complexes, QM/MM, or systems requiring mixed basis sets
- Use **GPAW** for rapid DFT calculations, Python integration, or real-space grid methods
- Use **QE** for plane-wave DFT, periodic systems, and electronic structure
- Use **LAMMPS** for large systems, classical MD, or reactive force fields (ReaxFF)
- Use **GROMACS** for biomolecular systems, enzymes, or solvated reactions

**Optimal Workflow for Screening**: MACE rapid screening → DFT validation of top candidates (10-100x speedup!)

**Metal Catalyst Workflow**: MACE-MP screening of metal structures → CP2K/GPAW DFT validation with XYZ coordinates"""
            },
            {
                "role": "user",
                "content": user_request
            }
        ]

        max_iterations = 30
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
                    if function_name == "run_quantum_espresso":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("quantum_espresso", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "Quantum ESPRESSO"}
                        print(f"Result: Quantum ESPRESSO job submitted with ID: {job_id}")

                    elif function_name == "run_cp2k":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("cp2k", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "CP2K"}
                        print(f"Result: CP2K job submitted with ID: {job_id}")

                    elif function_name == "run_gpaw":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("gpaw", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "GPAW"}
                        print(f"Result: GPAW job submitted with ID: {job_id}")

                    elif function_name == "run_lammps":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("lammps", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "LAMMPS"}
                        print(f"Result: LAMMPS job submitted with ID: {job_id}")

                    elif function_name == "run_gromacs":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("gromacs", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "GROMACS"}
                        print(f"Result: GROMACS job submitted with ID: {job_id}")

                    elif function_name == "check_simulation_status":
                        result = self.check_status(arguments["job_id"])
                        print(f"Result: Status = {result['status']}")
                        if result.get('result'):
                            print(f"Simulation Results: {json.dumps(result['result'], indent=2)}")

                    elif function_name == "mace_predict_energy":
                        smiles = arguments["smiles"]
                        model_type = arguments.get("model_type", "off")
                        model_size = arguments.get("model_size", "medium")

                        print(f"Calling MACE for energy prediction: {smiles} (model: {model_type}, size: {model_size})")

                        try:
                            xyz = smiles_to_xyz(smiles)
                            mace_result = self.mace_client.predict_energy(
                                xyz,
                                format="xyz",
                                model_type=model_type,
                                model_size=model_size
                            )

                            if mace_result:
                                result = {
                                    "status": "SUCCESS",
                                    "smiles": smiles,
                                    "energy_ev": mace_result.get("energy"),
                                    "unit": mace_result.get("unit", "eV"),
                                    "model": mace_result.get("model"),
                                    "method": "MACE ML prediction"
                                }
                                print(f"Result: MACE predicted energy = {mace_result.get('energy')} {mace_result.get('unit', 'eV')}")
                            else:
                                result = {"status": "FAILED", "error": "MACE prediction failed"}
                                print("Result: MACE prediction FAILED")
                        except Exception as e:
                            result = {"status": "FAILED", "error": str(e)}
                            print(f"Result: Error - {e}")

                    elif function_name == "mace_rapid_screening":
                        smiles_list = arguments["smiles_list"]
                        model_type = arguments.get("model_type", "off")
                        top_n = arguments.get("top_n", None)

                        print(f"Calling MACE for rapid screening: {len(smiles_list)} molecules (model: {model_type})")

                        try:
                            from mace_client_simple import rapid_screening

                            screening_results = rapid_screening(smiles_list, model_type=model_type)

                            if screening_results:
                                # Apply top_n filter if requested
                                if top_n and top_n < len(screening_results):
                                    screening_results = screening_results[:top_n]

                                result = {
                                    "status": "SUCCESS",
                                    "total_screened": len(smiles_list),
                                    "results_returned": len(screening_results),
                                    "ranked_candidates": screening_results,
                                    "method": "MACE ML batch screening"
                                }
                                print(f"Result: Screened {len(smiles_list)} molecules, returning top {len(screening_results)}")
                                print(f"Best candidate: {screening_results[0]['smiles']} ({screening_results[0]['energy']:.3f} eV)")
                            else:
                                result = {"status": "FAILED", "error": "MACE screening failed"}
                                print("Result: MACE screening FAILED")
                        except Exception as e:
                            result = {"status": "FAILED", "error": str(e)}
                            print(f"Result: Error - {e}")

                    elif function_name == "mace_optimize_geometry":
                        smiles = arguments["smiles"]
                        model_type = arguments.get("model_type", "off")
                        fmax = arguments.get("fmax", 0.05)
                        steps = arguments.get("steps", 200)

                        print(f"Calling MACE for geometry optimization: {smiles} (model: {model_type}, fmax: {fmax})")

                        try:
                            xyz = smiles_to_xyz(smiles)
                            opt_result = self.mace_client.optimize_geometry(
                                xyz,
                                format="xyz",
                                model_type=model_type,
                                fmax=fmax,
                                steps=steps
                            )

                            if opt_result:
                                result = {
                                    "status": "SUCCESS",
                                    "smiles": smiles,
                                    "initial_energy_ev": opt_result.get("initial_energy"),
                                    "final_energy_ev": opt_result.get("final_energy"),
                                    "converged": opt_result.get("converged", False),
                                    "steps_taken": opt_result.get("steps"),
                                    "optimized_structure": opt_result.get("optimized_structure"),
                                    "method": "MACE ML optimization"
                                }
                                print(f"Result: Optimized in {opt_result.get('steps')} steps")
                                print(f"  Initial: {opt_result.get('initial_energy'):.3f} eV → Final: {opt_result.get('final_energy'):.3f} eV")
                            else:
                                result = {"status": "FAILED", "error": "MACE optimization failed"}
                                print("Result: MACE optimization FAILED")
                        except Exception as e:
                            result = {"status": "FAILED", "error": str(e)}
                            print(f"Result: Error - {e}")

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
    """Main entry point for the computational chemistry agentic workflow"""

    print("Initializing Computational Chemistry Agent...")
    agent = ComputationalChemistryAgent(server_name="spark-container-03")  # Using gpt-oss:120b

    # Example user request
    user_request = """I need to study ethanol (CCO) molecule for catalyst design.
    Please run a DFT calculation to analyze its electronic structure and molecular properties.
    Use appropriate computational chemistry methods."""

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
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
