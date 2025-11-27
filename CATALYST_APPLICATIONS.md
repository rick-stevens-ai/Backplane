# Computational Chemistry Applications for Catalyst Design

This document describes the integration of five major computational chemistry applications into the Scientific AI Backplane system for catalyst design workflows.

## Overview

The Backplane now supports five computational chemistry codes through a unified agentic workflow:

1. **Quantum ESPRESSO (QE)** - Plane-wave DFT for electronic structure
2. **CP2K** - Mixed Gaussian/plane-wave DFT, excellent for QM/MM
3. **GPAW** - Real-space grid DFT in Python
4. **LAMMPS** - Classical molecular dynamics with ReaxFF
5. **GROMACS** - Biomolecular molecular dynamics

## Architecture

```
User Request (Natural Language)
         ↓
gpt-oss:120b Agent (agent_apps.py)
         ↓
FastAPI Endpoint (/submit_app_job)
         ↓
Celery Task Queue (run_simulation)
         ↓
Application Wrapper (Python)
         ↓
Computational Chemistry Code
         ↓
Results Returned to Agent
```

### Key Components

- **`agent_apps.py`**: Extended agentic workflow with tools for all five applications
- **`wrappers/`**: Python wrappers providing unified interface for each application
- **`tasks.py`**: Celery tasks for async job execution
- **`main.py`**: FastAPI endpoints for job submission and status

## Application Wrappers

All wrappers inherit from `SimulationWrapper` base class in `wrappers/base_wrapper.py`.

### Common Interface

Each wrapper provides:
- **Input file generation**: Application-specific input from JSON parameters
- **Job execution**: Subprocess management with timeouts
- **Output parsing**: Extract scientific results from output files
- **Error handling**: Capture and report failures

### Application-Specific Features

#### 1. Quantum ESPRESSO (`quantum_espresso.py`)

**Use Cases:**
- Electronic structure calculations
- Band structure and density of states
- Geometry optimization
- Phonon calculations

**Key Parameters:**
```python
{
    "calculation": "scf" | "relax" | "vc-relax" | "md",
    "molecule_smiles": "CCO",  # Ethanol
    "cutoff_wfc": 50.0,        # Rydberg
    "cutoff_rho": 400.0        # Rydberg
}
```

**Output:**
- Energy (Hartree and eV)
- Forces on atoms
- Stress tensor
- Convergence status

#### 2. CP2K (`cp2k_wrapper.py`)

**Use Cases:**
- Mixed quantum/classical (QM/MM) simulations
- Large system DFT calculations
- Surface catalysis
- Condensed phase systems

**Key Parameters:**
```python
{
    "run_type": "ENERGY" | "GEO_OPT" | "MD" | "CELL_OPT",
    "molecule_smiles": "CCO",
    "functional": "PBE" | "BLYP" | "B3LYP",
    "cutoff": 400.0  # Rydberg
}
```

**Output:**
- Total energy (Hartree and eV)
- Forces (atomic units)
- SCF convergence information
- Timing data

#### 3. GPAW (`gpaw_wrapper.py`)

**Use Cases:**
- Rapid DFT prototyping
- Python-integrated workflows
- Real-space calculations
- Nanostructures and molecules

**Key Parameters:**
```python
{
    "calculation_type": "energy" | "relaxation" | "md",
    "molecule_smiles": "O",  # Water
    "mode": "fd" | "pw" | "lcao",  # Finite difference, plane wave, or LCAO
    "xc": "PBE",
    "h": 0.2  # Grid spacing in Angstrom
}
```

**Output:**
- Energy (eV and Hartree)
- Forces (eV/Angstrom)
- Dipole moment
- Final atomic positions

#### 4. LAMMPS (`lammps_wrapper.py`)

**Use Cases:**
- Large-scale molecular dynamics
- Reactive force fields (ReaxFF)
- Catalyst surface dynamics
- Microkinetic modeling

**Key Parameters:**
```python
{
    "simulation_type": "energy" | "minimize" | "md" | "reaxff",
    "molecule_smiles": "CCO",
    "force_field": "lj" | "reaxff" | "eam",
    "temperature": 300.0,  # Kelvin
    "timestep": 1.0,       # Femtoseconds
    "steps": 1000
}
```

**Output:**
- Potential energy (kcal/mol and eV)
- Kinetic energy
- Temperature
- Pressure and volume
- Thermodynamic data

#### 5. GROMACS (`gromacs_wrapper.py`)

**Use Cases:**
- Biomolecular simulations
- Enzyme catalysis
- Solvated reactions
- Protein-ligand interactions

**Key Parameters:**
```python
{
    "simulation_type": "em" | "nvt" | "npt" | "md",
    "molecule_smiles": "O",
    "temperature": 300.0,    # Kelvin
    "pressure": 1.0,         # Bar (for NPT)
    "timestep": 0.002,       # Picoseconds
    "steps": 50000,
    "force_field": "oplsaa"
}
```

**Output:**
- Potential energy (kJ/mol and kcal/mol)
- Temperature and pressure
- Convergence status (for minimization)
- Simulation time

## Agent Tool Definitions

The agent (`agent_apps.py`) provides five `run_*` functions that the LLM can call:

### Tool 1: `run_quantum_espresso`
```json
{
  "experiment_name": "Ethanol Electronic Structure",
  "calculation": "relax",
  "molecule_smiles": "CCO",
  "cutoff_wfc": 50,
  "cutoff_rho": 400
}
```

### Tool 2: `run_cp2k`
```json
{
  "experiment_name": "Ethanol QM/MM Study",
  "run_type": "GEO_OPT",
  "molecule_smiles": "CCO",
  "functional": "PBE",
  "cutoff": 400
}
```

### Tool 3: `run_gpaw`
```json
{
  "experiment_name": "Water Molecule DFT",
  "calculation_type": "relaxation",
  "molecule_smiles": "O",
  "mode": "fd",
  "xc": "PBE",
  "h": 0.2
}
```

### Tool 4: `run_lammps`
```json
{
  "experiment_name": "Ethanol MD Simulation",
  "simulation_type": "minimize",
  "molecule_smiles": "CCO",
  "force_field": "lj",
  "temperature": 300,
  "steps": 1000
}
```

### Tool 5: `run_gromacs`
```json
{
  "experiment_name": "Water Biomolecular MD",
  "simulation_type": "em",
  "molecule_smiles": "O",
  "temperature": 300,
  "force_field": "oplsaa"
}
```

### Tool 6: `check_simulation_status`
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000"
}
```

## Usage Examples

### Example 1: Basic Usage

```python
from agent_apps import ComputationalChemistryAgent

# Initialize agent with gpt-oss:120b
agent = ComputationalChemistryAgent(
    server_config_path="spark_servers.yaml",
    server_name="spark-container-03"
)

# Submit natural language request
result = agent.run_agentic_workflow("""
Please analyze ethanol (CCO) for catalyst design using quantum chemistry.
Calculate its electronic structure and molecular properties.
""")

print(result)
```

### Example 2: Application-Specific Request

```python
# Request specific application
result = agent.run_agentic_workflow("""
Use CP2K to perform a geometry optimization on ethanol (CCO)
using the PBE functional. I need accurate energies for catalyst screening.
""")
```

### Example 3: Comparative Study

```python
# Agent will choose appropriate applications
result = agent.run_agentic_workflow("""
Compare the electronic properties of ethanol (CCO) using both
DFT and classical MD methods. Use appropriate codes for each approach.
""")
```

## Testing

Run the comprehensive test suite:

```bash
python test_catalyst_applications.py
```

This tests all five applications with gpt-oss:120b, demonstrating:
- Autonomous application selection
- Parameter setting
- Job submission and monitoring
- Results interpretation

## Installation Status

Check which applications are installed:

```python
from wrappers.quantum_espresso import QuantumEspressoWrapper
from wrappers.cp2k_wrapper import CP2KWrapper
from wrappers.gpaw_wrapper import GPAWWrapper
from wrappers.lammps_wrapper import LAMMPSWrapper
from wrappers.gromacs_wrapper import GROMACSWrapper

wrappers = [
    QuantumEspressoWrapper(),
    CP2KWrapper(),
    GPAWWrapper(),
    LAMMPSWrapper(),
    GROMACSWrapper()
]

for wrapper in wrappers:
    print(f"{wrapper.app_name}: version {wrapper.version}")
```

## Scientific Workflow Recommendations

### For Catalyst Design:

1. **Initial Screening (LAMMPS)**: Use classical MD with force fields for rapid screening of large catalyst libraries

2. **Electronic Structure (QE or CP2K)**: Perform DFT calculations on promising candidates for accurate energies and band gaps

3. **Reaction Mechanisms (GPAW)**: Use Python integration for custom reaction coordinate calculations

4. **QM/MM Studies (CP2K)**: Model active sites with quantum accuracy while treating environment classically

5. **Enzyme Catalysis (GROMACS)**: Simulate protein dynamics and ligand binding for biocatalysis

### Best Practices:

- **Start simple**: Use small molecules for testing before scaling up
- **Validate**: Compare results across multiple codes when possible
- **Monitor resources**: Check job execution times and adjust parameters
- **Iterate**: Use agent's ability to adapt based on initial results

## Troubleshooting

### Common Issues:

**1. Application not found:**
```
version: "not_installed"
```
Solution: Check APPS directory or system PATH for executables

**2. Simulation timeout:**
```
status: "timeout"
```
Solution: Increase timeout parameter or reduce system size

**3. Convergence failure:**
```
converged: False
```
Solution: Adjust convergence criteria, cutoffs, or starting geometry

**4. Import errors:**
```
ModuleNotFoundError: No module named 'wrappers'
```
Solution: Ensure wrappers directory is in Python path

## API Reference

### FastAPI Endpoints

**POST `/submit_app_job`**
```json
{
  "application": "quantum_espresso",
  "job_params": {
    "calculation": "scf",
    "molecule_smiles": "O"
  },
  "experiment_name": "Water DFT"
}
```

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "QUEUED"
}
```

**GET `/job_status/{job_id}`**

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "SUCCESS",
  "result": {
    "application": "quantum_espresso",
    "results": {
      "energy_ry": -34.12,
      "energy_ev": -464.68,
      "converged": true
    }
  }
}
```

## Performance Considerations

### Resource Requirements:

| Application | Typical Time | Memory | CPU Cores |
|------------|--------------|--------|-----------|
| QE (small) | 1-10 min | 1-4 GB | 1-8 |
| CP2K (small) | 2-20 min | 2-8 GB | 1-16 |
| GPAW (small) | 1-5 min | 1-2 GB | 1-4 |
| LAMMPS (small) | <1 min | <1 GB | 1-4 |
| GROMACS (small) | 1-5 min | <1 GB | 1-8 |

### Scaling:

- **Molecules**: Up to ~100 atoms per application
- **QM/MM**: Use CP2K for systems >100 atoms
- **MD**: LAMMPS and GROMACS can handle millions of atoms
- **Production**: Deploy to HPC with SLURM integration

## Future Enhancements

1. **SLURM Integration**: Replace subprocess with actual HPC job submission
2. **RDKit Integration**: Automatic SMILES to 3D structure conversion
3. **Result Analysis**: Automated property extraction and visualization
4. **Multi-step Workflows**: Chain calculations across applications
5. **ML Potentials**: Integration with neural network potentials

## References

- Quantum ESPRESSO: https://www.quantum-espresso.org
- CP2K: https://www.cp2k.org
- GPAW: https://wiki.fysik.dtu.dk/gpaw
- LAMMPS: https://lammps.org
- GROMACS: https://www.gromacs.org

## Last Updated

November 21, 2025
