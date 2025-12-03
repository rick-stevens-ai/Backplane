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
from backplane_logging import get_logger, log_section, log_subsection, log_dict

# Initialize logger for agent
logger = get_logger('agent')


def format_result_for_display(result: Any, indent: int = 0) -> str:
    """Format a result (dict, list, or scalar) in a human-friendly way"""
    prefix = "  " * indent

    if result is None:
        return f"{prefix}(none)"

    if isinstance(result, dict):
        lines = []
        for key, value in result.items():
            if isinstance(value, (dict, list)):
                lines.append(f"{prefix}{key}:")
                lines.append(format_result_for_display(value, indent + 1))
            elif isinstance(value, str) and len(value) > 100:
                lines.append(f"{prefix}{key}: {value[:100]}... (truncated)")
            elif isinstance(value, float):
                lines.append(f"{prefix}{key}: {value:.6f}")
            else:
                lines.append(f"{prefix}{key}: {value}")
        return "\n".join(lines)

    elif isinstance(result, list):
        if len(result) == 0:
            return f"{prefix}(empty list)"
        elif len(result) > 5:
            lines = [f"{prefix}[Showing first 5 of {len(result)} items]"]
            for i, item in enumerate(result[:5]):
                lines.append(f"{prefix}  [{i}] {item}")
            return "\n".join(lines)
        else:
            lines = []
            for i, item in enumerate(result):
                if isinstance(item, dict):
                    lines.append(f"{prefix}  [{i}]:")
                    lines.append(format_result_for_display(item, indent + 2))
                else:
                    lines.append(f"{prefix}  [{i}] {item}")
            return "\n".join(lines)

    else:
        return f"{prefix}{result}"


# Tool descriptions for logging
TOOL_DESCRIPTIONS = {
    "run_quantum_espresso": "Submit DFT calculation using Quantum ESPRESSO (plane-wave basis set)",
    "run_cp2k": "Submit DFT/QM-MM calculation using CP2K (mixed Gaussian/plane-wave basis)",
    "run_gpaw": "Submit DFT calculation using GPAW (real-space grid or plane-wave basis)",
    "run_lammps": "Submit molecular dynamics simulation using LAMMPS (classical force fields)",
    "run_gromacs": "Submit molecular dynamics simulation using GROMACS (biomolecular systems)",
    "check_simulation_status": "Check the status and retrieve results of a submitted simulation job",
    "mace_predict_energy": "Fast ML energy prediction using MACE foundation model (~0.5s)",
    "mace_optimize_geometry": "ML-based geometry optimization using MACE foundation model (~2s)",
    "mace_rapid_screening": "Rapid screening of multiple molecules using MACE ML batch prediction"
}


# Estimated execution times for logging
TOOL_TIME_ESTIMATES = {
    "run_quantum_espresso": "30-120 seconds (DFT calculation on HPC)",
    "run_cp2k": "30-120 seconds (DFT calculation on HPC)",
    "run_gpaw": "30-120 seconds (DFT calculation on HPC)",
    "run_lammps": "10-60 seconds (MD simulation, depends on steps)",
    "run_gromacs": "10-60 seconds (MD simulation, depends on steps)",
    "check_simulation_status": "~1 second (status query)",
    "mace_predict_energy": "~0.5 seconds (fast ML prediction)",
    "mace_optimize_geometry": "~2 seconds (ML optimization)",
    "mace_rapid_screening": "~0.5-2 seconds (batch ML prediction, depends on number of molecules)"
}


class ComputationalChemistryAgent:
    """Agent that submits and monitors computational chemistry simulations across multiple applications"""

    def __init__(self, server_config_path: str = "spark_servers.yaml", server_name: str = "spark-container-03"):
        """
        Initialize the agent with a specific server configuration

        Args:
            server_config_path: Path to the YAML configuration file
            server_name: Name of the server to use from the config
        """
        logger.info("Initializing Computational Chemistry Agent")
        logger.info(f"  Server config path: {server_config_path}")
        logger.info(f"  Server name: {server_name}")

        self.api_base = "http://127.0.0.1:8000"  # FastAPI server
        logger.info(f"  API base URL: {self.api_base}")

        self.server_config = self._load_server_config(server_config_path, server_name)
        logger.info(f"  LLM server: {self.server_config['openai_api_base']}")
        logger.info(f"  LLM model: {self.server_config['openai_model']}")

        self.client = OpenAI(
            api_key=self.server_config["openai_api_key"],
            base_url=self.server_config["openai_api_base"]
        )
        self.model = self.server_config["openai_model"]

        # Initialize MACE client for rapid ML predictions
        logger.info("  Initializing MACE ML client")
        self.mace_client = MACEClient()
        logger.info("Agent initialization complete")

    def _load_server_config(self, config_path: str, server_name: str) -> Dict[str, str]:
        """Load server configuration from YAML file"""
        logger.debug(f"Loading server config from {config_path}")
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        for server in config['servers']:
            if server['server'] == server_name:
                logger.debug(f"Found server config for {server_name}")
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
        logger.info(f"Submitting job to {application}")
        logger.info(f"  Experiment: {experiment_name}")
        log_dict(logger, job_params, "  Job parameters")

        job_data = {
            "application": application,
            "job_params": job_params,
            "experiment_name": experiment_name
        }

        # Use the new run_simulation endpoint
        logger.debug(f"Sending POST request to {self.api_base}/submit_app_job")
        response = requests.post(
            f"{self.api_base}/submit_app_job",
            json=job_data
        )
        response.raise_for_status()
        result = response.json()
        job_id = result['job_id']
        logger.info(f"  Job submitted successfully! Job ID: {job_id}")
        return job_id

    def check_status(self, job_id: str) -> Dict[str, Any]:
        """Check the status of a submitted job"""
        logger.debug(f"Checking status of job {job_id}")
        response = requests.get(f"{self.api_base}/job_status/{job_id}")
        response.raise_for_status()
        status = response.json()
        logger.debug(f"  Job {job_id} status: {status.get('status', 'UNKNOWN')}")
        return status

    def run_agentic_workflow(self, user_request: str) -> Dict[str, Any]:
        """
        Run the full agentic workflow with support for multiple computational chemistry applications

        Args:
            user_request: Natural language request from user

        Returns:
            Dictionary with simulation results
        """
        log_section(logger, "COMPUTATIONAL CHEMISTRY AGENTIC WORKFLOW")
        logger.info("User Request:")
        logger.info(f"  {user_request}")
        logger.info("")

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

**CRITICAL REQUIREMENT: YOU MUST THINK OUT LOUD AND EXPLAIN YOUR REASONING**

⚠️ **MANDATORY**: You are REQUIRED to provide explanatory text BEFORE every tool call.
DO NOT call tools without first explaining your reasoning in plain text.

Before calling any tools, you MUST provide a clear explanation of your computational strategy:

1. **Analyze the Problem**: What is the scientific question? What properties need to be calculated?

2. **Deliberate on Tool Selection**:
   - What tools are available for this task?
   - Which tool is most appropriate and why?
   - What are the trade-offs (speed vs. accuracy)?
   - Should I use MACE for rapid screening first, or go directly to DFT?

3. **Explain Computational Order**:
   - What sequence of calculations is needed?
   - Why this order? (e.g., "First MACE screening to narrow candidates, then DFT validation of top 3")
   - What will each step accomplish?

4. **Justify Parameter Choices**:
   - Why these specific parameters (functional, cutoff, run_type)?
   - What assumptions am I making?

**Example of Good Deliberation**:
"The user wants to screen 20 catalyst candidates for NH3 activation. Let me think about the best approach:

ANALYSIS: This is a screening task with many molecules, so computational efficiency is key.

STRATEGY: I should use a two-stage approach:
1. MACE rapid screening (all 20 molecules, ~10s total) to identify promising candidates
2. DFT validation (top 3 candidates, ~3 min each) for accurate energetics

RATIONALE:
- MACE can quickly filter out poor candidates (saves ~85 min vs. DFT for all 20)
- CP2K is best for the DFT stage because these are metal-organic complexes
- I'll use the PBE functional (good balance of accuracy/speed for catalysis)

ORDER OF EXECUTION:
1. Call mace_rapid_screening with all 20 SMILES
2. Analyze results, select top 3
3. Submit CP2K jobs for detailed analysis of top candidates
4. Monitor with check_simulation_status
5. Compare and recommend best catalyst"

Guidelines for choosing applications:
- Use **MACE** for rapid screening, pre-filtering, or quick feasibility checks
- Use **CP2K** for metal clusters, coordination complexes, QM/MM, or systems requiring mixed basis sets
- Use **GPAW** for rapid DFT calculations, Python integration, or real-space grid methods
- Use **QE** for plane-wave DFT, periodic systems, and electronic structure
- Use **LAMMPS** for large systems, classical MD, or reactive force fields (ReaxFF)
- Use **GROMACS** for biomolecular systems, enzymes, or solvated reactions

**Optimal Workflow for Screening**: MACE rapid screening → DFT validation of top candidates (10-100x speedup!)

**Metal Catalyst Workflow**: MACE-MP screening of metal structures → CP2K/GPAW DFT validation with XYZ coordinates

**RESPONSE FORMAT - MANDATORY**:
Your response MUST follow this format:
1. First provide explanatory text with your analysis and reasoning (200-500 words)
2. Then make your tool call(s)

❌ **INCORRECT** - Tool call without explanation:
<tool_call>run_gpaw(...)</tool_call>

✅ **CORRECT** - Explanation followed by tool call:
"Let me analyze this water molecule calculation request...
[200+ words of analysis, reasoning, and strategy]
Now I'll proceed with the calculation."
<tool_call>run_gpaw(...)</tool_call>

**If you skip the explanatory text, your response will be considered invalid!**"""
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
            log_subsection(logger, f"ITERATION {iteration}")
            logger.info(f"Sending request to LLM ({self.model})")
            logger.info(f"  Message history length: {len(messages)}")

            # Log the latest user/assistant message for context
            logger.info("")
            logger.info("=== PROMPT TO LLM ===")
            for msg in messages[-3:]:  # Show last 3 messages for context
                # Handle both dict messages and ChatCompletionMessage objects
                if isinstance(msg, dict):
                    role = msg.get('role', 'unknown')
                    content = msg.get('content', '')
                else:
                    # ChatCompletionMessage object
                    role = getattr(msg, 'role', 'unknown')
                    content = getattr(msg, 'content', '') or ''

                if content:
                    logger.info(f"[{role.upper()}]:")
                    for line in content.split('\n'):
                        if line.strip():
                            logger.info(f"  {line}")
                if hasattr(msg, 'tool_calls') and msg.tool_calls:
                    logger.info(f"[{role.upper()}]: <tool calls - see below>")
            logger.info("=== END PROMPT ===")
            logger.info("")

            # Get agent's response
            llm_start = time.time()
            response = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                tools=tools,
                tool_choice="auto"
            )
            llm_elapsed = time.time() - llm_start

            message = response.choices[0].message
            messages.append(message)

            logger.info(f"LLM response received in {llm_elapsed:.2f}s")
            logger.info("")

            # Check for reasoning/thinking content (for reasoning models like oss120b)
            reasoning_content = None
            if hasattr(message, 'reasoning'):
                reasoning_content = message.reasoning
            elif hasattr(message, 'thinking'):
                reasoning_content = message.thinking
            elif hasattr(message, 'reasoning_content'):
                reasoning_content = message.reasoning_content
            # Also check in the raw response object
            elif hasattr(response.choices[0], 'reasoning'):
                reasoning_content = response.choices[0].reasoning

            # Display reasoning steps if available
            if reasoning_content:
                logger.info("=== MODEL REASONING / THINKING ===")
                logger.info("[INTERNAL REASONING]:")
                for line in str(reasoning_content).split('\n'):
                    if line.strip():
                        logger.info(f"  💭 {line}")
                logger.info("=== END REASONING ===")
                logger.info("")

            logger.info("=== LLM RESPONSE ===")

            # Log the response content if present
            if message.content:
                logger.info("[ASSISTANT TEXT]:")
                for line in message.content.split('\n'):
                    logger.info(f"  {line}")
            elif message.tool_calls:
                # Show which tools are being called instead of generic message
                tool_names = [tc.function.name for tc in message.tool_calls]
                logger.warning("⚠️  LLM DID NOT PROVIDE EXPLANATORY TEXT")
                logger.warning("   The LLM is supposed to explain its reasoning before calling tools!")
                logger.warning("   This violates the deliberation requirements in the system prompt.")
                logger.info(f"[ASSISTANT]: Proceeding directly to tool execution (no explanation provided)")
                logger.info(f"  Tools to call: {', '.join(tool_names)}")
            else:
                logger.info("[ASSISTANT]: <empty response>")

            logger.info("=== END RESPONSE ===")
            logger.info("")

            # Check if agent wants to use tools
            if message.tool_calls:
                logger.info(f"Agent is calling {len(message.tool_calls)} tool(s)")

                for tool_call in message.tool_calls:
                    function_name = tool_call.function.name
                    arguments = json.loads(tool_call.function.arguments)

                    logger.info("")
                    logger.info("="*80)
                    logger.info(f"🔧 TOOL CALL: {function_name}")
                    logger.info("="*80)

                    # Log tool description and time estimate
                    description = TOOL_DESCRIPTIONS.get(function_name, "No description available")
                    time_estimate = TOOL_TIME_ESTIMATES.get(function_name, "Unknown")

                    logger.info(f"📋 Purpose: {description}")
                    logger.info(f"⏱️  Estimated time: {time_estimate}")
                    logger.info("")
                    logger.info("📥 Input Parameters:")
                    for key, value in arguments.items():
                        if isinstance(value, str) and len(value) > 200:
                            logger.info(f"  • {key}: {value[:200]}... (truncated)")
                        elif isinstance(value, list) and len(value) > 5:
                            logger.info(f"  • {key}: {value[:5]}... ({len(value)} items total)")
                        else:
                            logger.info(f"  • {key}: {value}")
                    logger.info("")
                    logger.info("▶️  Executing...")

                    # Execute the requested function
                    if function_name == "run_quantum_espresso":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("quantum_espresso", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "Quantum ESPRESSO"}
                        logger.info(f"🚀 Result: Quantum ESPRESSO job submitted with ID: {job_id}")

                    elif function_name == "run_cp2k":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("cp2k", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "CP2K"}
                        logger.info(f"🚀 Result: CP2K job submitted with ID: {job_id}")

                    elif function_name == "run_gpaw":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("gpaw", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "GPAW"}
                        logger.info(f"🚀 Result: GPAW job submitted with ID: {job_id}")

                    elif function_name == "run_lammps":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("lammps", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "LAMMPS"}
                        logger.info(f"🚀 Result: LAMMPS job submitted with ID: {job_id}")

                    elif function_name == "run_gromacs":
                        experiment_name = arguments.pop("experiment_name")
                        job_id = self.submit_app_job("gromacs", arguments, experiment_name)
                        result = {"job_id": job_id, "status": "QUEUED", "application": "GROMACS"}
                        logger.info(f"🚀 Result: GROMACS job submitted with ID: {job_id}")

                    elif function_name == "check_simulation_status":
                        result = self.check_status(arguments["job_id"])
                        logger.info(f"✅ Result: Status = {result['status']}")
                        if result.get('result'):
                            logger.info("")
                            logger.info("📊 Simulation Results:")
                            logger.info(format_result_for_display(result['result'], indent=1))

                    elif function_name == "mace_predict_energy":
                        smiles = arguments["smiles"]
                        model_type = arguments.get("model_type", "off")
                        model_size = arguments.get("model_size", "medium")

                        logger.info(f"Calling MACE for energy prediction: {smiles}")
                        logger.info(f"  Model: {model_type}, Size: {model_size}")

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
                                logger.info(f"⚡ Result: MACE predicted energy = {mace_result.get('energy')} {mace_result.get('unit', 'eV')}")
                            else:
                                result = {"status": "FAILED", "error": "MACE prediction failed"}
                                logger.warning("Result: MACE prediction FAILED")
                        except Exception as e:
                            result = {"status": "FAILED", "error": str(e)}
                            logger.error(f"Result: Error - {e}")

                    elif function_name == "mace_rapid_screening":
                        smiles_list = arguments["smiles_list"]
                        model_type = arguments.get("model_type", "off")
                        top_n = arguments.get("top_n", None)

                        logger.info(f"Calling MACE for rapid screening: {len(smiles_list)} molecules")
                        logger.info(f"  Model: {model_type}, Top N: {top_n or 'all'}")

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
                                logger.info(f"🔍 Result: Screened {len(smiles_list)} molecules, returning top {len(screening_results)}")
                                logger.info(f"🏆 Best candidate: {screening_results[0]['smiles']} ({screening_results[0]['energy']:.3f} eV)")
                            else:
                                result = {"status": "FAILED", "error": "MACE screening failed"}
                                logger.warning("Result: MACE screening FAILED")
                        except Exception as e:
                            result = {"status": "FAILED", "error": str(e)}
                            logger.error(f"Result: Error - {e}")

                    elif function_name == "mace_optimize_geometry":
                        smiles = arguments["smiles"]
                        model_type = arguments.get("model_type", "off")
                        fmax = arguments.get("fmax", 0.05)
                        steps = arguments.get("steps", 200)

                        logger.info(f"Calling MACE for geometry optimization: {smiles}")
                        logger.info(f"  Model: {model_type}, fmax: {fmax}, max steps: {steps}")

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
                                steps = opt_result.get('steps')
                                initial_e = opt_result.get('initial_energy')
                                final_e = opt_result.get('final_energy')

                                if steps is not None:
                                    logger.info(f"🔄 Result: Optimized in {steps} steps")
                                else:
                                    logger.info("🔄 Result: Optimization completed (steps not reported)")

                                if initial_e is not None and final_e is not None:
                                    logger.info(f"  📈 Initial: {initial_e:.3f} eV → Final: {final_e:.3f} eV")
                                elif final_e is not None:
                                    logger.info(f"  📊 Final energy: {final_e:.3f} eV")
                            else:
                                result = {"status": "FAILED", "error": "MACE optimization failed"}
                                logger.warning("Result: MACE optimization FAILED")
                        except Exception as e:
                            result = {"status": "FAILED", "error": str(e)}
                            logger.error(f"Result: Error - {e}")

                    else:
                        logger.error(f"Unknown function: {function_name}")
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
                log_section(logger, "AGENT FINAL RESPONSE")
                for line in message.content.split('\n'):
                    logger.info(f"  {line}")
                logger.info("="*80)

                # If we have a job_id and agent is done, return the final result
                if job_id:
                    logger.info("Retrieving final job status...")
                    final_status = self.check_status(job_id)
                    logger.info(f"Final status: {final_status.get('status')}")
                    return final_status

                logger.info("Workflow complete!")
                return {"message": message.content}

        logger.error("Agent exceeded maximum iterations without completing task!")
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
