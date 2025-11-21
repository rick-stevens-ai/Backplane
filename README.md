# Scientific AI Backplane

An asynchronous job submission system for HPC/supercomputer workloads, designed to prevent AI agents from timing out on long-running scientific simulations.

## Architecture

- **FastAPI Server** - REST API with job submission and status endpoints
- **Celery Worker** - Processes jobs asynchronously with progress tracking
- **Redis** - Message broker and result backend for Celery queue
- **Pydantic Models** - Data validation schemas

## Features

- Immediate job submission with unique job IDs
- Asynchronous processing of long-running simulations
- Progress tracking and status polling
- Support for DFT, MD, and docking simulations
- SLURM/HPC scheduler integration ready

## Installation

```bash
# Install dependencies
pip install fastapi uvicorn celery redis pydantic

# Install and start Redis (macOS)
brew install redis
brew services start redis

# Or run Redis via Docker
docker run -d -p 6379:6379 redis
```

## Usage

```bash
# Terminal 1: Start Celery worker
celery -A tasks.celery_app worker --loglevel=info

# Terminal 2: Start FastAPI server
uvicorn main:app --reload
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

Returns:
```json
{
  "job_id": "0d27a6f3-a117-4bfa-8b0e-3fc00bb51c86",
  "status": "SUCCESS",
  "result": {
    "energy_level": -420.5,
    "convergence": true,
    "output_path": "/scratch/users/ai_factory/0d27a6f3-a117-4bfa-8b0e-3fc00bb51c86.chk"
  }
}
```

## Example

```bash
# Submit a test job
curl -X POST http://127.0.0.1:8000/submit_job \
  -H "Content-Type: application/json" \
  -d @json_ex.json

# Poll for status
curl http://127.0.0.1:8000/job_status/{job_id}
```

## Workflow

```
AI Agent → POST /submit_job → Celery Queue (Redis) → Worker executes → AI Agent polls /job_status → Gets results
```

This prevents AI agents from timing out on long-running scientific simulations by returning immediately with a job ID and allowing asynchronous processing.

## License

MIT
