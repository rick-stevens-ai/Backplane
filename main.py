from fastapi import FastAPI
from celery.result import AsyncResult
from pydantic import BaseModel
from typing import Dict, Any
from models import SimulationRequest, JobStatus
from tasks import submit_slurm_job, run_simulation

app = FastAPI(title="Scientific AI Backplane")


class AppJobRequest(BaseModel):
    """Request model for application-specific simulation jobs"""
    application: str
    job_params: Dict[str, Any]
    experiment_name: str = "simulation"

@app.post("/submit_job", response_model=JobStatus)
async def submit_job(request: SimulationRequest):
    """
    Endpoint for the AI Agent.
    Returns a Job ID immediately, preventing the Agent from timing out.
    """
    # 1. Validate data (handled by Pydantic)
    
    # 2. Offload to Celery Worker
    # .delay() is the magic method that sends the task to Redis
    task = submit_slurm_job.delay(request.model_dump())
    
    return JobStatus(job_id=task.id, status="QUEUED")

@app.get("/job_status/{job_id}", response_model=JobStatus)
async def get_status(job_id: str):
    """
    Agent polls this endpoint to check if the simulation is done.
    """
    task_result = AsyncResult(job_id)
    
    response = JobStatus(
        job_id=job_id,
        status=task_result.status
    )
    
    if task_result.successful():
        response.result = task_result.result

    return response


@app.post("/submit_app_job", response_model=JobStatus)
async def submit_app_job(request: AppJobRequest):
    """
    Endpoint for submitting application-specific simulation jobs.
    Supports Quantum ESPRESSO, CP2K, GPAW, LAMMPS, and GROMACS.
    Returns a Job ID immediately.
    """
    # Offload to Celery Worker using run_simulation task
    task = run_simulation.delay(request.model_dump())

    return JobStatus(job_id=task.id, status="QUEUED")
