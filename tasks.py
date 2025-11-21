import time
import subprocess
from celery import Celery

# Connect to Redis (assuming localhost for this demo)
celery_app = Celery("hpc_worker", broker="redis://localhost:6379/0", backend="redis://localhost:6379/0")

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

    # Return final "scientific" result
    return {
        "energy_level": -420.5,
        "convergence": True,
        "output_path": f"/scratch/users/ai_factory/{job_id}.chk"
    }
