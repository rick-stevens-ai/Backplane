#!/usr/bin/env python3
"""Minimal agent test"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

print("Step 1: Starting test", flush=True)

from agent_apps import ComputationalChemistryAgent
print("Step 2: Imported agent", flush=True)

agent = ComputationalChemistryAgent(
    server_config_path="spark_servers.yaml",
    server_name="spark-container-03"
)
print(f"Step 3: Agent initialized - Model: {agent.model}", flush=True)

# Try a simple test
print("Step 4: Sending test request to agent...", flush=True)
result = agent.run_agentic_workflow("What is 2+2?")
print(f"Step 5: Got result: {result}", flush=True)

print("SUCCESS: Agent is working!", flush=True)
