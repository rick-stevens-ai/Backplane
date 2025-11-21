from pydantic import BaseModel, Field
from typing import Dict, Optional

class SimulationRequest(BaseModel):
    """
    The blueprint for a scientific simulation job.
    The Agent must strictly adhere to this schema.
    """
    experiment_name: str = Field(..., description="Name of the experiment for logging")
    molecule_smiles: str = Field(..., description="SMILES string of the molecule to analyze")
    simulation_type: str = Field("dft_optimization", description="Type of simulation: 'dft', 'md', or 'docking'")
    parameters: Dict[str, float] = Field(default_factory=dict, description="Physical parameters (e.g., {'temperature': 300.0})")

class JobStatus(BaseModel):
    job_id: str
    status: str
    result: Optional[Dict] = None
