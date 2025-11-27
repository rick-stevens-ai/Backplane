#!/usr/bin/env python3
"""Comprehensive test for the containerized server endpoint"""

import yaml
import requests
import json
from typing import Dict

def load_container_server() -> Dict:
    """Load the containerized server configuration"""
    with open("spark_servers.yaml", 'r') as f:
        data = yaml.safe_load(f)

    # Find the container server
    for server in data['servers']:
        if 'container' in server['server']:
            return server

    raise ValueError("No containerized server found in configuration")

def test_basic_connectivity(server: Dict) -> bool:
    """Test basic server connectivity"""
    print("\n" + "=" * 60)
    print("Test 1: Basic Connectivity")
    print("=" * 60)

    try:
        url = f"{server['openai_api_base']}/models"
        print(f"Testing: {url}")
        response = requests.get(url, timeout=10)

        print(f"Status Code: {response.status_code}")

        if response.status_code == 200:
            print("✅ Server is reachable")
            try:
                data = response.json()
                print(f"Response: {json.dumps(data, indent=2)}")
            except:
                print(f"Response (text): {response.text[:200]}")
            return True
        else:
            print(f"⚠️  Server returned status {response.status_code}")
            print(f"Response: {response.text[:200]}")
            return False
    except Exception as e:
        print(f"❌ Connection failed: {str(e)}")
        return False

def test_model_listing(server: Dict) -> bool:
    """Test listing available models"""
    print("\n" + "=" * 60)
    print("Test 2: Model Listing")
    print("=" * 60)

    try:
        url = f"{server['openai_api_base']}/models"
        print(f"Testing: {url}")
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            print("✅ Models endpoint accessible")
            try:
                data = response.json()
                if 'data' in data:
                    print(f"Available models: {len(data['data'])}")
                    for model in data['data']:
                        print(f"  - {model.get('id', 'Unknown')}")
                else:
                    print(f"Response: {json.dumps(data, indent=2)}")
            except:
                print(f"Response (text): {response.text[:500]}")
            return True
        else:
            print(f"⚠️  Status {response.status_code}")
            return False
    except Exception as e:
        print(f"❌ Failed: {str(e)}")
        return False

def test_completion_request(server: Dict) -> bool:
    """Test a simple completion request"""
    print("\n" + "=" * 60)
    print("Test 3: Completion Request")
    print("=" * 60)

    try:
        url = f"{server['openai_api_base']}/chat/completions"
        print(f"Testing: {url}")

        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {server['openai_api_key']}"
        }

        payload = {
            "model": server['openai_model'],
            "messages": [
                {"role": "user", "content": "Say 'Hello from the containerized server!' and nothing else."}
            ],
            "max_tokens": 50,
            "temperature": 0.7
        }

        print(f"Request payload: {json.dumps(payload, indent=2)}")

        response = requests.post(url, headers=headers, json=payload, timeout=30)

        print(f"Status Code: {response.status_code}")

        if response.status_code == 200:
            print("✅ Completion request successful")
            try:
                data = response.json()
                if 'choices' in data and len(data['choices']) > 0:
                    message = data['choices'][0].get('message', {}).get('content', '')
                    print(f"Model response: {message}")
                else:
                    print(f"Response: {json.dumps(data, indent=2)}")
            except:
                print(f"Response (text): {response.text[:500]}")
            return True
        else:
            print(f"⚠️  Status {response.status_code}")
            print(f"Response: {response.text[:500]}")
            return False
    except requests.exceptions.Timeout:
        print("⚠️  Request timed out (this is common for first request)")
        return False
    except Exception as e:
        print(f"❌ Failed: {str(e)}")
        return False

def test_health_check(server: Dict) -> bool:
    """Test various health check endpoints"""
    print("\n" + "=" * 60)
    print("Test 4: Health Check Endpoints")
    print("=" * 60)

    base_url = server['openai_api_base']
    endpoints = [
        "/health",
        "/v1/health",
        "",  # Base URL
    ]

    for endpoint in endpoints:
        url = f"{base_url}{endpoint}"
        try:
            print(f"\nTrying: {url}")
            response = requests.get(url, timeout=5)
            print(f"  Status: {response.status_code}")
            if response.status_code == 200:
                print(f"  ✅ Endpoint accessible")
                try:
                    print(f"  Response: {response.json()}")
                except:
                    print(f"  Response: {response.text[:100]}")
        except Exception as e:
            print(f"  ⚠️  {str(e)}")

    return True

def main():
    print("=" * 60)
    print("Containerized Server Comprehensive Test")
    print("=" * 60)

    try:
        server = load_container_server()
        print(f"\n📦 Testing containerized server:")
        print(f"   Server ID: {server['server']}")
        print(f"   Shortname: {server['shortname']}")
        print(f"   API Base: {server['openai_api_base']}")
        print(f"   Model: {server['openai_model']}")
    except Exception as e:
        print(f"❌ Failed to load configuration: {e}")
        return

    # Run tests
    results = []
    results.append(("Basic Connectivity", test_basic_connectivity(server)))
    results.append(("Model Listing", test_model_listing(server)))
    results.append(("Completion Request", test_completion_request(server)))
    results.append(("Health Check", test_health_check(server)))

    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    for test_name, passed in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status}: {test_name}")

    passed = sum(1 for _, p in results if p)
    total = len(results)
    print(f"\nResults: {passed}/{total} tests passed")
    print("=" * 60)

if __name__ == "__main__":
    main()
