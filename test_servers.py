#!/usr/bin/env python3
"""Test script to validate spark_servers.yaml configuration"""

import yaml
import requests
from typing import Dict, List

def load_server_config(yaml_file: str = "spark_servers.yaml") -> List[Dict]:
    """Load and parse the server configuration file"""
    with open(yaml_file, 'r') as f:
        data = yaml.safe_load(f)
    return data['servers']

def validate_server_config(server: Dict) -> bool:
    """Validate that a server configuration has all required fields"""
    required_fields = ['server', 'shortname', 'openai_api_key', 'openai_api_base', 'openai_model']

    for field in required_fields:
        if field not in server:
            print(f"❌ Missing required field: {field}")
            return False

    return True

def test_server_connectivity(server: Dict) -> bool:
    """Test if the server endpoint is reachable"""
    try:
        # Test if the base URL is accessible (models endpoint)
        url = f"{server['openai_api_base']}/models"
        response = requests.get(url, timeout=5)

        if response.status_code == 200:
            print(f"✅ {server['shortname']}: Server is reachable")
            return True
        else:
            print(f"⚠️  {server['shortname']}: Server returned status {response.status_code}")
            return False
    except requests.exceptions.ConnectionError:
        print(f"⚠️  {server['shortname']}: Connection refused (server may be offline)")
        return False
    except requests.exceptions.Timeout:
        print(f"⚠️  {server['shortname']}: Connection timeout")
        return False
    except Exception as e:
        print(f"⚠️  {server['shortname']}: Error - {str(e)}")
        return False

def main():
    print("=" * 60)
    print("Testing spark_servers.yaml Configuration")
    print("=" * 60)

    # Load configuration
    try:
        servers = load_server_config()
        print(f"\n✅ Successfully loaded {len(servers)} server(s)\n")
    except Exception as e:
        print(f"❌ Failed to load configuration: {e}")
        return

    # Validate each server
    print("-" * 60)
    print("Configuration Validation:")
    print("-" * 60)
    for i, server in enumerate(servers, 1):
        print(f"\nServer {i}: {server.get('shortname', 'Unknown')}")
        print(f"  Server ID: {server.get('server')}")
        print(f"  API Base: {server.get('openai_api_base')}")
        print(f"  Model: {server.get('openai_model')}")

        if validate_server_config(server):
            print(f"  ✅ Configuration is valid")
        else:
            print(f"  ❌ Configuration is invalid")

    # Test connectivity
    print("\n" + "-" * 60)
    print("Connectivity Testing:")
    print("-" * 60)
    print()

    for server in servers:
        test_server_connectivity(server)

    print("\n" + "=" * 60)
    print("Test Complete")
    print("=" * 60)

if __name__ == "__main__":
    main()
