import time
import subprocess
from celery import Celery
from pathlib import Path
import sys

# Import application wrappers
from wrappers.quantum_espresso import QuantumEspressoWrapper
from wrappers.cp2k_wrapper import CP2KWrapper
from wrappers.gpaw_wrapper import GPAWWrapper
from wrappers.lammps_wrapper import LAMMPSWrapper
from wrappers.gromacs_wrapper import GROMACSWrapper

# Connect to Redis (assuming localhost for this demo)
celery_app = Celery("hpc_worker", broker="redis://localhost:6379/0", backend="redis://localhost:6379/0")

# Initialize wrappers (cached for reuse)
WRAPPERS = {
    "quantum_espresso": None,
    "qe": None,
    "cp2k": None,
    "gpaw": None,
    "lammps": None,
    "gromacs": None
}

def get_wrapper(app_name: str):
    """Get or initialize a wrapper for the specified application"""
    app_name_lower = app_name.lower()

    # Map variations to canonical names
    if app_name_lower in ["quantum_espresso", "qe", "quantum-espresso"]:
        canonical_name = "quantum_espresso"
        wrapper_class = QuantumEspressoWrapper
    elif app_name_lower == "cp2k":
        canonical_name = "cp2k"
        wrapper_class = CP2KWrapper
    elif app_name_lower == "gpaw":
        canonical_name = "gpaw"
        wrapper_class = GPAWWrapper
    elif app_name_lower == "lammps":
        canonical_name = "lammps"
        wrapper_class = LAMMPSWrapper
    elif app_name_lower == "gromacs":
        canonical_name = "gromacs"
        wrapper_class = GROMACSWrapper
    else:
        raise ValueError(f"Unknown application: {app_name}")

    # Initialize wrapper if not already cached
    if WRAPPERS[canonical_name] is None:
        WRAPPERS[canonical_name] = wrapper_class()

    return WRAPPERS[canonical_name]

@celery_app.task(bind=True)
def submit_slurm_job(self, experiment_data: dict):
    """
    This task simulates submitting a job to a supercomputer (e.g., NERSC).
    """
    job_id = self.request.id
    sim_type = experiment_data.get("simulation_type")
    
    # ---------------------------------------------------------
    # REAL IMPLEMENTATION:
    # slurm_cmd = f"sbatch --job-name={job_id} run_sim.sh"
    # result = subprocess.run(slurm_cmd, shell=True, capture_output=True)
    # ---------------------------------------------------------
    
    # MOCK IMPLEMENTATION (for demonstration):
    print(f"[{job_id}] connecting to Slurm scheduler...")
    print(f"[{job_id}] allocating nodes for {sim_type}...")
    
    # Simulate long-running science (e.g., 10 seconds)
    total_steps = 5
    for step in range(total_steps):
        time.sleep(2) # Simulating computation
        # Update status so the API can query it
        self.update_state(state='PROGRESS', meta={'current': step, 'total': total_steps})
        print(f"[{job_id}] Simulation step {step+1}/{total_steps} complete.")

    # Return final "scientific" result with detailed information
    experiment_name = experiment_data.get("experiment_name", "Unknown")
    molecule = experiment_data.get("molecule_smiles", "N/A")
    params = experiment_data.get("parameters", {})

    return {
        "experiment_name": experiment_name,
        "molecule_smiles": molecule,
        "simulation_type": sim_type,
        "parameters": params,
        "results": {
            "energy_level_hartree": -420.5,
            "energy_level_ev": -11442.7,
            "convergence_achieved": True,
            "optimization_steps": total_steps,
            "final_geometry_rmsd": 0.0023,
            "dipole_moment_debye": 2.45,
            "total_charge": 0,
            "spin_multiplicity": 1,
            "vibrational_frequencies": [345.2, 567.8, 892.1, 1024.5, 1456.3],
            "homo_lumo_gap_ev": 4.82,
            "computation_time_seconds": total_steps * 2
        },
        "output_files": {
            "checkpoint": f"/scratch/users/ai_factory/{job_id}.chk",
            "log_file": f"/scratch/users/ai_factory/{job_id}.log",
            "geometry": f"/scratch/users/ai_factory/{job_id}.xyz"
        },
        "status": "completed_successfully",
        "completion_time": time.strftime("%Y-%m-%d %H:%M:%S")
    }


@celery_app.task(bind=True)
def run_simulation(self, job_data: dict):
    """
    Run a real computational chemistry simulation using application wrappers

    Expected job_data:
    - application: Name of the application ('quantum_espresso', 'cp2k', 'gpaw', 'lammps', 'gromacs')
    - job_params: Dictionary of application-specific parameters
    - experiment_name: Optional experiment name
    """
    job_id = self.request.id
    application = job_data.get("application")
    job_params = job_data.get("job_params", {})
    experiment_name = job_data.get("experiment_name", "simulation")

    print(f"[{job_id}] Starting {application} simulation: {experiment_name}")

    try:
        # Get the appropriate wrapper
        wrapper = get_wrapper(application)

        print(f"[{job_id}] Using {wrapper.app_name} (version: {wrapper.version})")

        # Update initial state
        self.update_state(
            state='PROGRESS',
            meta={
                'status': 'initializing',
                'application': application,
                'version': wrapper.version
            }
        )

        # Prepare job directory
        job_dir = wrapper.prepare_job(job_id, job_params)
        print(f"[{job_id}] Job directory created: {job_dir}")

        self.update_state(
            state='PROGRESS',
            meta={
                'status': 'running',
                'job_directory': str(job_dir)
            }
        )

        # Run the simulation
        result = wrapper.run_job(job_dir)

        # Check if successful
        if result['status'] == 'success':
            print(f"[{job_id}] Simulation completed successfully")
            result['job_id'] = job_id
            result['application'] = application
            result['experiment_name'] = experiment_name
            result['completion_time'] = time.strftime("%Y-%m-%d %H:%M:%S")
            return result
        else:
            print(f"[{job_id}] Simulation failed: {result.get('error', 'Unknown error')}")
            return {
                "job_id": job_id,
                "application": application,
                "experiment_name": experiment_name,
                "status": "failed",
                "error": result.get('error', 'Unknown error'),
                "error_details": result.get('error_details', ''),
                "completion_time": time.strftime("%Y-%m-%d %H:%M:%S")
            }

    except Exception as e:
        print(f"[{job_id}] Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return {
            "job_id": job_id,
            "application": application,
            "experiment_name": experiment_name,
            "status": "error",
            "error": str(e),
            "completion_time": time.strftime("%Y-%m-%d %H:%M:%S")
        }
