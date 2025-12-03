#!/usr/bin/env python3
"""
Test parallel agentic workflows for multiple molecules
Demonstrates concurrent calculations with enhanced logging
"""
import sys
import time
import threading
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from agent_apps import ComputationalChemistryAgent
from backplane_logging import setup_logging, get_logger

# Setup logging
setup_logging(level=20)  # INFO level
logger = get_logger('test')


def run_workflow(molecule_name, request, results_dict):
    """
    Run a workflow in a separate thread

    Args:
        molecule_name: Name for logging
        request: User request string
        results_dict: Shared dictionary to store results
    """
    thread_logger = get_logger(f'workflow-{molecule_name}')

    thread_logger.info("="*80)
    thread_logger.info(f"STARTING WORKFLOW: {molecule_name.upper()}")
    thread_logger.info("="*80)

    try:
        # Create agent instance
        agent = ComputationalChemistryAgent(
            server_config_path='spark_servers.yaml',
            server_name='spark-container-03'
        )

        # Run workflow
        start_time = time.time()
        result = agent.run_agentic_workflow(request)
        elapsed = time.time() - start_time

        # Store results
        results_dict[molecule_name] = {
            'result': result,
            'elapsed': elapsed,
            'success': True
        }

        thread_logger.info("")
        thread_logger.info("="*80)
        thread_logger.info(f"COMPLETED: {molecule_name.upper()} in {elapsed:.2f}s")
        thread_logger.info("="*80)

    except Exception as e:
        thread_logger.error(f"ERROR in {molecule_name}: {e}")
        results_dict[molecule_name] = {
            'result': None,
            'elapsed': 0,
            'success': False,
            'error': str(e)
        }


def main():
    logger.info("="*80)
    logger.info("    PARALLEL AGENTIC WORKFLOWS TEST")
    logger.info("="*80)
    logger.info("")
    logger.info("Testing concurrent computational chemistry workflows:")
    logger.info("  1. Ammonia (NH3) - Simple inorganic molecule")
    logger.info("  2. Ethanol (C2H5OH) - Small organic molecule")
    logger.info("")
    logger.info("Each workflow will:")
    logger.info("  - Deliberate on tool selection")
    logger.info("  - Choose between MACE (fast) and DFT (accurate)")
    logger.info("  - Execute calculations")
    logger.info("  - Report results")
    logger.info("")

    # Define workflows
    workflows = {
        'ammonia': """Calculate the energy and properties of ammonia (NH3, SMILES: N).
This is a simple nitrogen hydride molecule commonly used in fertilizers and as a chemical building block.
Please provide accurate energy calculations.""",

        'ethanol': """Calculate the energy and properties of ethanol (C2H5OH, SMILES: CCO).
This is a simple alcohol molecule used as a solvent and biofuel.
Please provide accurate energy calculations."""
    }

    # Shared results dictionary
    results = {}

    # Create threads
    threads = []
    logger.info("Launching parallel workflows...")
    logger.info("")

    start_time = time.time()

    for molecule, request in workflows.items():
        thread = threading.Thread(
            target=run_workflow,
            args=(molecule, request, results),
            name=f"workflow-{molecule}"
        )
        threads.append(thread)
        thread.start()
        logger.info(f"  ✓ Started workflow for {molecule}")

    logger.info("")
    logger.info("Waiting for workflows to complete...")
    logger.info("")

    # Wait for all threads to complete
    for thread in threads:
        thread.join()

    total_elapsed = time.time() - start_time

    # Print summary
    logger.info("")
    logger.info("="*80)
    logger.info("                    WORKFLOW SUMMARY")
    logger.info("="*80)
    logger.info("")

    for molecule, result_data in results.items():
        if result_data['success']:
            logger.info(f"{molecule.upper()}:")
            logger.info(f"  Status: ✓ SUCCESS")
            logger.info(f"  Time: {result_data['elapsed']:.2f}s")
            logger.info(f"  Result: {result_data['result'][:150]}..." if len(result_data['result']) > 150 else f"  Result: {result_data['result']}")
        else:
            logger.info(f"{molecule.upper()}:")
            logger.info(f"  Status: ✗ FAILED")
            logger.info(f"  Error: {result_data.get('error', 'Unknown error')}")
        logger.info("")

    logger.info(f"Total parallel execution time: {total_elapsed:.2f}s")
    logger.info("")

    # Calculate efficiency
    sequential_time = sum(r['elapsed'] for r in results.values() if r['success'])
    if total_elapsed > 0 and sequential_time > 0:
        speedup = sequential_time / total_elapsed
        logger.info(f"Sequential time would be: ~{sequential_time:.2f}s")
        logger.info(f"Speedup from parallelization: {speedup:.2f}x")

    logger.info("")
    logger.info("="*80)

    # Return success if all workflows succeeded
    all_success = all(r['success'] for r in results.values())
    return 0 if all_success else 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
