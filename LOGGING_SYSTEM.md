# Comprehensive Logging System

## Overview

The Scientific AI Backplane now has comprehensive, color-coded logging throughout all components. This logging system is specifically designed for:

- **Screen recording demonstrations** - Every step is logged with clear, readable output
- **Debugging** - Detailed execution traces at every level
- **Monitoring** - Track workflow progress in real-time
- **Understanding** - See exactly what the agent is thinking and doing

## Features

### Color-Coded Output
- **INFO** (Cyan): Normal operations, workflow progress
- **WARNING** (Yellow): Non-fatal issues
- **ERROR** (Red): Failures and exceptions
- **DEBUG** (Gray): Detailed technical information

### Logged Components

1. **Agent Workflow** (`agent_apps.py`)
   - Agent initialization and configuration
   - User request logging
   - LLM API calls with timing
   - LLM prompts and responses (full content)
   - Tool calls with arguments
   - Tool execution results
   - Final workflow results

2. **Celery Worker** (`tasks.py`)
   - Job initialization
   - Wrapper selection and initialization
   - Job directory creation
   - Simulation execution with timing
   - Success/failure status
   - Error handling

3. **Logging Module** (`backplane_logging.py`)
   - Centralized configuration
   - Color formatting
   - Section headers for visual organization
   - Dictionary logging helpers

## Testing

### Simple H2O Test (Recommended for First Demo)

```bash
# Make sure Redis and Celery worker are running
redis-server &
celery -A tasks worker --loglevel=info &

# Run the comprehensive logging test
python3 test_h2o_with_logging.py
```

This test will show:
- Agent initialization
- LLM request/response cycle
- Tool selection (will choose CP2K, GPAW, or QE for water)
- Job submission
- Status monitoring
- Final results

### Expected Output

You'll see color-coded output like:

```
[INFO    ] [test        ] ================================================================================
[INFO    ] [test        ]          WATER MOLECULE (H2O) DFT CALCULATION TEST
[INFO    ] [test        ] ================================================================================
[INFO    ] [agent       ] Initializing Computational Chemistry Agent
[INFO    ] [agent       ]   Server config path: spark_servers.yaml
[INFO    ] [agent       ]   LLM model: gpt-oss:120b
[INFO    ] [agent       ] ================================================================================
[INFO    ] [agent       ]     COMPUTATIONAL CHEMISTRY AGENTIC WORKFLOW
[INFO    ] [agent       ] ================================================================================
[INFO    ] [agent       ] User Request:
[INFO    ] [agent       ]   Calculate the energy and geometry of a water molecule...
[INFO    ] [agent       ]
[INFO    ] [agent       ] ITERATION 1
[INFO    ] [agent       ] --------------------------------------------------------------------------------
[INFO    ] [agent       ] Sending request to LLM (gpt-oss:120b)
[INFO    ] [agent       ] LLM response received in 8.45s
[INFO    ] [agent       ] Tool Call: run_cp2k
[INFO    ] [agent       ] Arguments:
[INFO    ] [agent       ]   experiment_name: water_energy_calculation
[INFO    ] [agent       ]   run_type: ENERGY
[INFO    ] [agent       ]   molecule_smiles: O
[INFO    ] [celery.worker] ================================================================================
[INFO    ] [celery.worker]              SIMULATION JOB: abc123...
[INFO    ] [celery.worker] ================================================================================
[INFO    ] [celery.worker] Application: cp2k
[INFO    ] [celery.worker] Using CP2K version 2023.2
[INFO    ] [celery.worker] Job directory created: /tmp/job_abc123
[INFO    ] [celery.worker] Starting simulation execution...
[INFO    ] [celery.worker] Simulation completed successfully in 15.23s
```

## Log Levels

### INFO (Default for Demos)
Perfect for screen recordings - shows all important steps without overwhelming detail.

```python
from backplane_logging import setup_logging
setup_logging(level=20)  # INFO level
```

### DEBUG (For Development)
Shows everything including internal details.

```python
setup_logging(level=10)  # DEBUG level
```

## File Locations

- **Logging Module**: `backplane_logging.py`
- **Agent Logging**: `agent_apps.py` (lines 14-17, throughout)
- **Worker Logging**: `tasks.py` (lines 14-22, throughout)
- **Test Script**: `test_h2o_with_logging.py`

## What's Logged

### Agent Workflow
- ✅ Agent initialization (server, model, MACE client)
- ✅ User request
- ✅ Each iteration number
- ✅ LLM API call timing
- ✅ LLM response content (full text)
- ✅ Tool calls with full arguments
- ✅ Tool execution timing and results
- ✅ Final workflow completion

### Simulation Jobs
- ✅ Job ID and metadata
- ✅ Application selection
- ✅ Wrapper initialization
- ✅ Job directory creation
- ✅ Execution timing
- ✅ Success/failure status
- ✅ Error messages and stack traces

## Screen Recording Tips

1. **Use INFO level** - Clean, readable output
2. **Terminal width**: Set to 120+ columns for best formatting
3. **Color scheme**: Use a dark terminal theme for best contrast
4. **Test run first**: The H2O example takes 1-2 minutes total
5. **Window title**: Shows "COMPUTATIONAL CHEMISTRY AGENTIC WORKFLOW" headers

## Troubleshooting

If you don't see colors:
- Check terminal supports ANSI colors
- Try: `export TERM=xterm-256color`

If logging is too verbose:
- Reduce level: `setup_logging(level=30)`  # WARNING level

If you need file logging:
- Add file: `setup_logging(level=20, log_file='workflow.log')`

## Next Steps

Once you've confirmed the logging works with the H2O test:
1. Try with larger molecules
2. Test MACE ML predictions (very fast, good for demos)
3. Run Fe(PNP)(N₂)(H)₂ catalyst workflow (comprehensive example)

All workflows now have comprehensive logging suitable for demonstrations!
