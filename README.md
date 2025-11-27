# Scientific AI Backplane

An asynchronous job submission system for HPC/supercomputer workloads, designed to prevent AI agents from timing out on long-running scientific simulations.

## Table of Contents

- [Architecture](#architecture)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [API Endpoints](#api-endpoints)
- [Agentic Workflows](#agentic-workflows)
- [Simulation Results](#simulation-results)
- [Examples](#examples)

## Architecture

- **FastAPI Server** (`main.py`) - REST API with job submission and status endpoints
- **Celery Worker** (`tasks.py`) - Processes jobs asynchronously with progress tracking
- **Redis** - Message broker and result backend for Celery queue
- **Pydantic Models** (`models.py`) - Data validation schemas
- **AI Agents** (`agent.py`, `agent_demo.py`) - Autonomous agents that submit and monitor simulations

## Features

- Immediate job submission with unique job IDs
- Asynchronous processing of long-running simulations
- Progress tracking and status polling
- Support for 10+ simulation codes (DFT, MD, ML)
- SLURM/HPC scheduler integration ready
- **Agentic AI workflows** - Autonomous agents that submit and monitor jobs
- **Detailed scientific results** - Energy levels, convergence, molecular properties
- **OpenAI-compatible** - Works with any OpenAI API-compatible model (Ollama, etc.)
- **Multiple input formats** - SMILES, XYZ, CIF for diverse chemical systems
- **Machine Learning** - MACE force field for fast pre-screening

## Supported Simulation Codes

### Quantum Chemistry (DFT/Wavefunction)
- **Quantum ESPRESSO** - Plane-wave DFT for periodic systems
- **CP2K** - Mixed Gaussian/plane-wave DFT and MD
- **GPAW** - Real-space DFT for molecules and surfaces
- **ORCA 6.0** - Molecular DFT and multireference methods ✨ NEW
- **NWChem 7.3** - Parallel quantum chemistry ✨ NEW
- **PySCF** - Python-based quantum chemistry with CASSCF/NEVPT2 ✨ NEW

### Classical Molecular Dynamics
- **LAMMPS** - Large-scale MD simulations
- **GROMACS** - Biomolecular MD

### Machine Learning
- **MACE** - ML force field for fast energy predictions ✨ NEW

### Workflow Automation
- **ASE** - Atomic Simulation Environment for NEB, optimization ✨ NEW

## Installation

**For complete installation instructions including all simulation codes, see [INSTALLATION.md](INSTALLATION.md)**

### Quick Start - Core Dependencies

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/backplane.git
cd backplane

# Install Python dependencies
pip install -r requirements.txt

# Install and start Redis (macOS)
brew install redis
brew services start redis

# Or run Redis via Docker
docker run -d -p 6379:6379 redis
```

### Simulation Codes

See [INSTALLATION.md](INSTALLATION.md) for detailed instructions on installing:
- Quantum chemistry codes (ORCA, NWChem, PySCF, Quantum ESPRESSO, CP2K, GPAW)
- Classical MD codes (LAMMPS, GROMACS)
- Machine learning models (MACE)
- Workflow tools (ASE)

## Quick Start

```bash
# Terminal 1: Start Celery worker
celery -A tasks.celery_app worker --loglevel=info

# Terminal 2: Start FastAPI server
uvicorn main:app --reload

# Terminal 3: Run the demo agent workflow
python agent_demo.py
```

## API Endpoints

### Submit Job
```bash
POST /submit_job
Content-Type: application/json

{
  "experiment_name": "Li-Ion Battery Test",
  "molecule_smiles": "Li[C6]",
  "simulation_type": "dft_optimization",
  "parameters": {"temperature": 300.0}
}
```

Returns:
```json
{
  "job_id": "0d27a6f3-a117-4bfa-8b0e-3fc00bb51c86",
  "status": "QUEUED"
}
```

### Check Status
```bash
GET /job_status/{job_id}
```

Returns detailed simulation results (see [Simulation Results](#simulation-results) section for full schema).

## Agentic Workflows

This system includes two autonomous agent implementations that can understand natural language requests, submit simulations, monitor progress, and interpret results.

### Agent Demo (`agent_demo.py`)

A demonstration agent that simulates AI reasoning without requiring an external LLM. Perfect for testing and understanding the workflow.

**Features:**
- Simulates step-by-step agent reasoning
- Parses natural language requests
- Submits jobs to the API
- Polls for completion
- Interprets scientific results

**Usage:**
```bash
python agent_demo.py
```

**Example Output:**
```
================================================================================
[AGENT REASONING] Step 1: Understanding User Request
================================================================================
I received a request to run a DFT optimization simulation.
Parsing the request...
- Molecule: Li[C6] (Lithium-ion battery compound)
- Simulation type: DFT optimization
- Experiment name: Battery Performance Analysis
- Temperature: 298.15 K
...

[MONITORING] Polling job status...
  [0.0s] Job 26bce258... status: PENDING
  [2.0s] Job 26bce258... status: PROGRESS
  [8.0s] Job 26bce258... status: PROGRESS

================================================================================
[AGENT FINAL RESPONSE]
================================================================================
I've successfully completed the DFT optimization simulation for your lithium-ion
battery molecule (Li[C6])...

**Energy Analysis:**
- Total Energy: -420.5 Hartree (-11442.7 eV)
- The calculation converged successfully after 5 optimization steps

**Electronic Properties:**
- HOMO-LUMO Gap: 4.82 eV
- Dipole Moment: 2.45 Debye
```

### Full Agent (`agent.py`)

A complete agentic implementation using OpenAI-compatible APIs with function calling.

**Features:**
- Connects to models from `spark_servers.yaml`
- Uses function calling for tool use
- Autonomous decision making
- Real LLM reasoning

**Configuration:**

Edit `spark_servers.yaml` to configure your LLM endpoint:
```yaml
servers:
  - server: "spark-36ae"
    shortname: "qwen3-coder:latest"
    openai_api_key: "CELS"
    openai_api_base: "http://192.168.1.212:11434/v1"
    openai_model: "qwen3-coder:latest"
```

**Usage:**
```bash
python agent.py
```

**How It Works:**

1. **Agent receives user request** (natural language)
   ```
   "Please run a DFT optimization simulation on a lithium-ion battery
    molecule (Li[C6]) for an experiment called 'Battery Performance Analysis'."
   ```

2. **Agent decides which tools to use**
   - Available tools: `submit_simulation`, `check_simulation_status`
   - Agent calls `submit_simulation` with parsed parameters

3. **Agent submits the job**
   ```python
   submit_simulation(
       experiment_name="Battery Performance Analysis",
       molecule_smiles="Li[C6]",
       simulation_type="dft_optimization",
       parameters={"temperature": 298.15}
   )
   ```

4. **Agent monitors progress**
   - Repeatedly calls `check_simulation_status(job_id)`
   - Continues until status is `SUCCESS`

5. **Agent interprets results**
   - Receives detailed scientific data
   - Explains findings in human-readable format

### Agent Workflow Diagram

```
User Request (Natural Language)
        ↓
   [AI Agent]
   - Parse request
   - Identify parameters
        ↓
   submit_simulation()
        ↓
   Job ID returned
        ↓
   [Polling Loop]
   ├─ check_simulation_status()
   ├─ Wait 2 seconds
   └─ Repeat until SUCCESS
        ↓
   Detailed Results
        ↓
   [AI Agent]
   - Interpret results
   - Format for user
        ↓
   Human-readable summary
```

## Simulation Results

The simulation worker returns comprehensive scientific data:

```json
{
  "job_id": "26bce258-a349-4ac1-9491-4a32526a7ced",
  "status": "SUCCESS",
  "result": {
    "experiment_name": "Battery Performance Analysis",
    "molecule_smiles": "Li[C6]",
    "simulation_type": "dft_optimization",
    "parameters": {
      "temperature": 298.15
    },
    "results": {
      "energy_level_hartree": -420.5,
      "energy_level_ev": -11442.7,
      "convergence_achieved": true,
      "optimization_steps": 5,
      "final_geometry_rmsd": 0.0023,
      "dipole_moment_debye": 2.45,
      "total_charge": 0,
      "spin_multiplicity": 1,
      "vibrational_frequencies": [345.2, 567.8, 892.1, 1024.5, 1456.3],
      "homo_lumo_gap_ev": 4.82,
      "computation_time_seconds": 10
    },
    "output_files": {
      "checkpoint": "/scratch/users/ai_factory/26bce258-a349-4ac1-9491-4a32526a7ced.chk",
      "log_file": "/scratch/users/ai_factory/26bce258-a349-4ac1-9491-4a32526a7ced.log",
      "geometry": "/scratch/users/ai_factory/26bce258-a349-4ac1-9491-4a32526a7ced.xyz"
    },
    "status": "completed_successfully",
    "completion_time": "2025-11-21 12:59:36"
  }
}
```

### Result Fields Explained

| Field | Description |
|-------|-------------|
| `energy_level_hartree` | Total energy in Hartree atomic units |
| `energy_level_ev` | Total energy in electron volts |
| `convergence_achieved` | Whether the optimization converged |
| `optimization_steps` | Number of optimization iterations |
| `final_geometry_rmsd` | Root mean square deviation of final geometry (Å) |
| `dipole_moment_debye` | Molecular dipole moment (Debye) |
| `total_charge` | Total molecular charge |
| `spin_multiplicity` | Spin multiplicity of the system |
| `vibrational_frequencies` | Calculated vibrational modes (cm⁻¹) |
| `homo_lumo_gap_ev` | HOMO-LUMO energy gap (eV) |
| `computation_time_seconds` | Wall-clock time for computation |

## Examples

### Example 1: Direct API Usage

```bash
# Submit a test job
curl -X POST http://127.0.0.1:8000/submit_job \
  -H "Content-Type: application/json" \
  -d @json_ex.json

# Poll for status
curl http://127.0.0.1:8000/job_status/{job_id}
```

### Example 2: Python Client

```python
import requests
import time

# Submit job
response = requests.post(
    "http://127.0.0.1:8000/submit_job",
    json={
        "experiment_name": "Water Molecule Analysis",
        "molecule_smiles": "O",
        "simulation_type": "dft_optimization",
        "parameters": {"temperature": 298.15}
    }
)
job_id = response.json()["job_id"]

# Poll until complete
while True:
    status_response = requests.get(f"http://127.0.0.1:8000/job_status/{job_id}")
    status = status_response.json()

    if status["status"] == "SUCCESS":
        print("Results:", status["result"])
        break

    time.sleep(2)
```

### Example 3: Agentic Workflow

```bash
# Run the demo agent
python agent_demo.py

# Or run the full agent (requires Ollama/LLM server)
python agent.py
```

The agent will:
1. Parse your natural language request
2. Submit the simulation automatically
3. Monitor progress
4. Interpret and explain results

## System Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                         User / AI Agent                         │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             │ Natural Language Request
                             ↓
                    ┌────────────────┐
                    │   AI Agent     │
                    │  (agent.py)    │
                    │                │
                    │ - Parse request│
                    │ - Call tools   │
                    └───────┬────────┘
                            │
                            │ POST /submit_job
                            ↓
                   ┌─────────────────┐
                   │  FastAPI Server │
                   │   (main.py)     │
                   │                 │
                   │ - Validate data │
                   │ - Return job_id │
                   └────────┬────────┘
                            │
                            │ Queue job
                            ↓
                    ┌───────────────┐
                    │     Redis     │
                    │ Message Broker│
                    └───────┬───────┘
                            │
                            │ Dequeue job
                            ↓
                   ┌─────────────────┐
                   │ Celery Worker   │
                   │   (tasks.py)    │
                   │                 │
                   │ - Run simulation│
                   │ - Update status │
                   │ - Store results │
                   └────────┬────────┘
                            │
                            │ Results stored in Redis
                            ↓
                    ┌───────────────┐
                    │     Redis     │
                    │ Result Backend│
                    └───────┬───────┘
                            │
                            │ GET /job_status/{job_id}
                            ↓
                   ┌─────────────────┐
                   │  FastAPI Server │
                   │ Returns results │
                   └────────┬────────┘
                            │
                            │ Detailed scientific results
                            ↓
                    ┌────────────────┐
                    │   AI Agent     │
                    │ Interprets &   │
                    │ explains       │
                    └───────┬────────┘
                            │
                            │ Human-readable summary
                            ↓
                ┌───────────────────────┐
                │        User           │
                └───────────────────────┘
```

## Key Benefits

1. **No Timeouts** - Long-running simulations don't block the API
2. **Async Processing** - Multiple jobs can run concurrently
3. **Progress Tracking** - Monitor simulation progress in real-time
4. **AI-Friendly** - Designed for autonomous agents with function calling
5. **Scalable** - Add more workers to handle increased load
6. **Rich Results** - Comprehensive scientific data returned
7. **Human-Readable** - AI agents interpret and explain results

## Use Cases

- **AI Research Assistants** - Autonomous agents that run simulations based on natural language requests
- **High-Throughput Screening** - Submit thousands of molecular simulations
- **Interactive Notebooks** - Jupyter/Python environments that need async job submission
- **Workflow Orchestration** - Part of larger scientific pipelines
- **HPC Integration** - Interface between AI agents and supercomputers

## File Structure

```
.
├── agent.py              # Full agentic workflow with LLM
├── agent_demo.py         # Demo agentic workflow (no LLM required)
├── main.py               # FastAPI REST API server
├── tasks.py              # Celery worker tasks
├── models.py             # Pydantic data models
├── spark_servers.yaml    # LLM server configuration
├── json_ex.json          # Example job submission
├── run_first.sh          # Start Celery worker
├── run_second.sh         # Start FastAPI server
└── README.md             # This file
```

## Contributing

Contributions welcome! Areas for improvement:
- Real SLURM integration
- Additional simulation types
- Authentication/authorization
- Result visualization
- Monitoring dashboard
- Webhook notifications

## License

MIT
