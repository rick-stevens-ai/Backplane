#!/usr/bin/env python3
"""
Test script for nitric oxide (NO) calculation with comprehensive logging demonstration.
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
    """Run NO nitric oxide molecule test with full logging"""

    log_section(logger, "NITRIC OXIDE (NO) DFT CALCULATION TEST")
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
    request = """Calculate the energy and properties of nitric oxide (NO).
SMILES: [N]=O

Nitric oxide is a radical molecule with one unpaired electron. It plays crucial
roles in biological signaling, atmospheric chemistry, and industrial processes.
The N-O bond has partial multiple bond character.

WORKFLOW REQUIREMENTS:
1. First, run MACE ML prediction for rapid baseline assessment
2. Then, run DFT calculation(s) for high-accuracy validation
3. Finally, compare and contrast the MACE and DFT results

Please compute:
- Total energy (MACE baseline + DFT validation)
- Optimized geometry
- Electronic structure properties
- N-O bond length and bond order
- Spin state and magnetic properties

At the end, provide a detailed comparison of:
- MACE vs DFT energy predictions
- Geometry differences
- Speed vs accuracy trade-offs

This molecule is important for understanding radical chemistry, NO₂ pollution,
and biological NO signaling pathways."""

    logger.info(f"Request: {request}")
    logger.info("")

    # Run workflow
    logger.info("Step 3: Running agentic workflow")
    logger.info("-" * 80)
    logger.info("The agent will:")
    logger.info("  1. Analyze the request and deliberate on tool selection")
    logger.info("  2. Consider whether to use MACE (fast ML) or DFT (accurate)")
    logger.info("  3. Note that NO is a radical with spin considerations")
    logger.info("  4. Execute the calculation with appropriate parameters")
    logger.info("  5. Return results with scientific interpretation")
    logger.info("")

    start_time = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed_time = time.time() - start_time

    # Display results
    log_section(logger, "FINAL RESULTS")
    logger.info(f"Total workflow time: {elapsed_time:.2f}s ({elapsed_time/60:.2f} minutes)")
    logger.info("")
    logger.info("Result summary:")
    logger.info(f"  Result: {result}")
    logger.info("")

    log_section(logger, "TEST COMPLETE")
    logger.info("Nitric oxide (NO) calculation complete!")
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
