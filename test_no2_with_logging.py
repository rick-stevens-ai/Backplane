#!/usr/bin/env python3
"""
Test script for nitrogen dioxide (NO2) calculation with comprehensive logging demonstration.
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
    """Run NO2 nitrogen dioxide molecule test with full logging"""

    log_section(logger, "NITROGEN DIOXIDE (NO2) DFT CALCULATION TEST")
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
    request = """Calculate the energy and properties of nitrogen dioxide (NO2).
SMILES: N(=O)=O

Nitrogen dioxide is a bent radical molecule with one unpaired electron. It is a major
air pollutant formed from combustion processes and is a key component of photochemical
smog. The molecule has resonance structures with delocalized bonding.

WORKFLOW REQUIREMENTS:
1. First, run MACE ML prediction for rapid baseline assessment
2. Then, run DFT calculation(s) for high-accuracy validation
3. Finally, compare and contrast the MACE and DFT results

Please compute:
- Total energy (MACE baseline + DFT validation)
- Optimized geometry (bent structure)
- Electronic structure properties
- N-O bond lengths and O-N-O bond angle
- Spin state and magnetic properties
- Molecular orbital structure showing π-system delocalization

At the end, provide a detailed comparison of:
- MACE vs DFT energy predictions
- Geometry differences
- Speed vs accuracy trade-offs

This molecule is important for understanding:
- Atmospheric chemistry and air pollution
- Photochemical smog formation
- NOₓ emissions and catalytic converters
- Radical chemistry in combustion"""

    logger.info(f"Request: {request}")
    logger.info("")

    # Run workflow
    logger.info("Step 3: Running agentic workflow")
    logger.info("-" * 80)
    logger.info("The agent will:")
    logger.info("  1. Analyze the request and deliberate on tool selection")
    logger.info("  2. Consider whether to use MACE (fast ML) or DFT (accurate)")
    logger.info("  3. Note that NO2 is a radical with bent geometry")
    logger.info("  4. Handle spin state and electronic structure properly")
    logger.info("  5. Execute the calculation with appropriate parameters")
    logger.info("  6. Return results with scientific interpretation")
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
    logger.info("Nitrogen dioxide (NO2) calculation complete!")
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
