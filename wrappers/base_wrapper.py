#!/usr/bin/env python3
"""
Base wrapper class for computational chemistry simulation codes
Provides unified interface for all application wrappers
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List
from pathlib import Path
import subprocess
import tempfile
import shutil
import json
import time
import os


class SimulationWrapper(ABC):
    """Abstract base class for simulation code wrappers"""

    def __init__(self, app_path: Optional[str] = None, scratch_dir: Optional[str] = None):
        """
        Initialize simulation wrapper

        Args:
            app_path: Path to application installation directory
            scratch_dir: Directory for temporary files and job execution
        """
        self.app_path = Path(app_path) if app_path else self._detect_app_path()
        self.scratch_dir = Path(scratch_dir) if scratch_dir else Path(tempfile.gettempdir()) / "backplane_jobs"
        self.scratch_dir.mkdir(parents=True, exist_ok=True)

        # Application metadata
        self.app_name = self.__class__.__name__.replace("Wrapper", "")
        self.version = None
        self._validate_installation()

    @abstractmethod
    def _detect_app_path(self) -> Path:
        """Detect application installation path"""
        pass

    @abstractmethod
    def _validate_installation(self) -> bool:
        """Validate that the application is properly installed"""
        pass

    @abstractmethod
    def _generate_input_file(self, job_params: Dict[str, Any], input_file: Path) -> None:
        """
        Generate application-specific input file

        Args:
            job_params: Job parameters from user request
            input_file: Path where input file should be written
        """
        pass

    @abstractmethod
    def _parse_output_file(self, output_file: Path, job_dir: Path) -> Dict[str, Any]:
        """
        Parse application output and extract results

        Args:
            output_file: Path to main output file
            job_dir: Path to job directory (for accessing other output files)

        Returns:
            Dictionary containing parsed results
        """
        pass

    @abstractmethod
    def _get_run_command(self, input_file: Path, output_file: Path, job_dir: Path) -> List[str]:
        """
        Get command to execute the simulation

        Args:
            input_file: Path to input file
            output_file: Path to output file
            job_dir: Path to job directory

        Returns:
            Command as list of strings (for subprocess)
        """
        pass

    def prepare_job(self, job_id: str, job_params: Dict[str, Any]) -> Path:
        """
        Prepare job directory and input files

        Args:
            job_id: Unique job identifier
            job_params: Job parameters from user request

        Returns:
            Path to job directory
        """
        # Create job directory
        job_dir = self.scratch_dir / job_id
        job_dir.mkdir(parents=True, exist_ok=True)

        # Generate input file
        input_file = job_dir / self._get_input_filename()
        self._generate_input_file(job_params, input_file)

        # Save job metadata
        metadata = {
            "job_id": job_id,
            "app_name": self.app_name,
            "app_version": self.version,
            "parameters": job_params,
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "status": "prepared"
        }
        metadata_file = job_dir / "job_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        return job_dir

    def run_job(self, job_dir: Path, timeout: int = 3600) -> Dict[str, Any]:
        """
        Execute the simulation job

        Args:
            job_dir: Path to job directory
            timeout: Maximum execution time in seconds (default: 1 hour)

        Returns:
            Dictionary containing job status and results
        """
        input_file = job_dir / self._get_input_filename()
        output_file = job_dir / self._get_output_filename()
        log_file = job_dir / "job.log"
        error_file = job_dir / "job.err"

        # Get run command
        cmd = self._get_run_command(input_file, output_file, job_dir)

        # Update metadata
        metadata_file = job_dir / "job_metadata.json"
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        metadata["status"] = "running"
        metadata["started_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        # Execute simulation
        start_time = time.time()
        try:
            result = subprocess.run(
                cmd,
                cwd=job_dir,
                stdout=open(log_file, 'w'),
                stderr=open(error_file, 'w'),
                timeout=timeout,
                check=False
            )
            execution_time = time.time() - start_time

            # Check if simulation succeeded
            if result.returncode != 0:
                error_msg = error_file.read_text() if error_file.exists() else "Unknown error"
                return {
                    "status": "failed",
                    "error": f"Simulation failed with return code {result.returncode}",
                    "error_details": error_msg[:1000],  # First 1000 chars
                    "execution_time": execution_time
                }

            # Parse output
            parsed_results = self._parse_output_file(output_file, job_dir)

            # Update metadata
            metadata["status"] = "completed"
            metadata["completed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
            metadata["execution_time"] = execution_time
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)

            return {
                "status": "success",
                "execution_time": execution_time,
                "results": parsed_results,
                "output_files": {
                    "input": str(input_file),
                    "output": str(output_file),
                    "log": str(log_file),
                    "metadata": str(metadata_file)
                }
            }

        except subprocess.TimeoutExpired:
            execution_time = time.time() - start_time
            return {
                "status": "timeout",
                "error": f"Simulation exceeded timeout of {timeout} seconds",
                "execution_time": execution_time
            }
        except Exception as e:
            execution_time = time.time() - start_time
            return {
                "status": "error",
                "error": str(e),
                "execution_time": execution_time
            }

    def submit_job(self, job_params: Dict[str, Any], job_id: Optional[str] = None) -> Dict[str, Any]:
        """
        High-level method to prepare and run a job

        Args:
            job_params: Job parameters from user request
            job_id: Optional job ID (generated if not provided)

        Returns:
            Dictionary containing job results
        """
        if job_id is None:
            import uuid
            job_id = str(uuid.uuid4())

        # Prepare job
        job_dir = self.prepare_job(job_id, job_params)

        # Run job
        result = self.run_job(job_dir)
        result["job_id"] = job_id
        result["job_directory"] = str(job_dir)

        return result

    def _get_input_filename(self) -> str:
        """Get default input filename (can be overridden)"""
        return "input.in"

    def _get_output_filename(self) -> str:
        """Get default output filename (can be overridden)"""
        return "output.out"

    def cleanup_job(self, job_dir: Path) -> None:
        """Clean up job directory"""
        if job_dir.exists():
            shutil.rmtree(job_dir)

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Get status of a job

        Args:
            job_id: Job identifier

        Returns:
            Dictionary containing job status
        """
        job_dir = self.scratch_dir / job_id
        metadata_file = job_dir / "job_metadata.json"

        if not metadata_file.exists():
            return {"status": "not_found", "job_id": job_id}

        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        return metadata

    def __repr__(self) -> str:
        return f"{self.app_name}Wrapper(app_path={self.app_path}, version={self.version})"
