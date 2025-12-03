#!/usr/bin/env python3
"""
Test script for water (H2O) calculation with comprehensive logging demonstration.
This script shows all logging output for screen recording purposes.
"""

import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent
from backplane_logging import setup_logging, get_logger, log_section
import time

# Setup logging at INFO level for demonstration
setup_logging(level=20)  # INFO level

# Get logger for this test script
logger = get_logger('test')

def main():
    """Run H2O water molecule test with full logging"""

    log_section(logger, "WATER MOLECULE (H2O) DFT CALCULATION TEST")
    logger.info("This test demonstrates the agentic workflow with comprehensive logging")
    logger.info("Perfect for screen recording demonstrations!")
    logger.info("")

    # Initialize agent
    logger.info("Step 1: Initializing Computational Chemistry Agent")
    logger.info("-" * 80)
    agent = ComputationalChemistryAgent(
        server_config_path='spark_servers.yaml',
        server_name='spark-container-03'
    )
    logger.info("")

    # Create request
    logger.info("Step 2: Creating user request")
    logger.info("-" * 80)
    request = """Calculate the energy and geometry of a water molecule (H2O).
SMILES: O
Use DFT to compute the total energy and optimize the geometry.
This is a simple test molecule."""

    logger.info(f"Request: {request}")
    logger.info("")

    # Run workflow
    logger.info("Step 3: Running agentic workflow")
    logger.info("-" * 80)
    logger.info("The agent will:")
    logger.info("  1. Parse the user request")
    logger.info("  2. Decide which simulation tool to use")
    logger.info("  3. Submit the job to the Celery worker")
    logger.info("  4. Monitor the job status")
    logger.info("  5. Return results")
    logger.info("")

    start_time = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed_time = time.time() - start_time

    # Display results
    log_section(logger, "FINAL RESULTS")
    logger.info(f"Total workflow time: {elapsed_time:.2f}s ({elapsed_time/60:.2f} minutes)")
    logger.info("")
    logger.info("Result summary:")
    if isinstance(result, dict):
        for key, value in result.items():
            if key != 'result':  # Skip nested result dict for clarity
                logger.info(f"  {key}: {value}")

    logger.info("")
    log_section(logger, "TEST COMPLETE")
    logger.info("All logging demonstrated successfully!")
    logger.info("Ready for screen recording!")

    return result


if __name__ == "__main__":
    try:
        result = main()
        sys.exit(0)
    except KeyboardInterrupt:
        logger.warning("Test interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
