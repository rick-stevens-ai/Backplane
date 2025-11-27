#!/usr/bin/env python3
"""Test script to discover the correct API endpoint for the containerized server"""

import requests
import json

BASE_HOST = "http://100.94.58.120:12000"
UUID = "5831e6f2-98c7-42b8-b75e-918d07923706"

# Try various endpoint patterns
ENDPOINT_PATTERNS = [
    f"{BASE_HOST}/v1",
    f"{BASE_HOST}/api/v1",
    f"{BASE_HOST}/api",
    f"{BASE_HOST}/openai/v1",
    f"{BASE_HOST}",
    f"{BASE_HOST}/ollama/v1",
    f"{BASE_HOST}/ollama/api",
]

def test_endpoint(base_url: str) -> dict:
    """Test an endpoint to see if it's a valid OpenAI-compatible API"""
    results = {
        'base_url': base_url,
        'models_endpoint': None,
        'chat_endpoint': None,
        'working': False
    }

    print(f"\n{'='*60}")
    print(f"Testing: {base_url}")
    print('='*60)

    # Test /models endpoint
    try:
        models_url = f"{base_url}/models"
        print(f"  Trying: {models_url}")
        response = requests.get(models_url, timeout=5)
        print(f"    Status: {response.status_code}")

        if response.status_code == 200:
            try:
                data = response.json()
                if 'data' in data or 'models' in data:
                    print(f"    ✅ Valid models endpoint (JSON)")
                    print(f"    Response keys: {list(data.keys())}")
                    results['models_endpoint'] = 'valid'
                    results['working'] = True

                    # Show available models
                    if 'data' in data:
                        print(f"    Models found: {len(data['data'])}")
                        for model in data['data'][:3]:  # Show first 3
                            print(f"      - {model.get('id', 'unknown')}")
                    elif 'models' in data:
                        print(f"    Models: {data['models'][:3]}")
                else:
                    print(f"    ⚠️  JSON but unexpected format: {list(data.keys())}")
                    results['models_endpoint'] = 'unexpected_json'
            except json.JSONDecodeError:
                print(f"    ⚠️  Not JSON (probably HTML)")
                results['models_endpoint'] = 'html'
        else:
            print(f"    ❌ Status {response.status_code}")
            results['models_endpoint'] = f'status_{response.status_code}'

    except requests.exceptions.ConnectionError:
        print(f"    ❌ Connection refused")
        results['models_endpoint'] = 'connection_refused'
    except requests.exceptions.Timeout:
        print(f"    ⚠️  Timeout")
        results['models_endpoint'] = 'timeout'
    except Exception as e:
        print(f"    ❌ Error: {str(e)}")
        results['models_endpoint'] = 'error'

    # If models endpoint worked, try chat completion
    if results['working']:
        try:
            chat_url = f"{base_url}/chat/completions"
            print(f"\n  Trying: {chat_url}")

            headers = {
                "Content-Type": "application/json",
                "Authorization": "Bearer CELS"
            }

            payload = {
                "model": "gpt-oss:20b",
                "messages": [{"role": "user", "content": "test"}],
                "max_tokens": 5
            }

            response = requests.post(chat_url, headers=headers, json=payload, timeout=10)
            print(f"    Status: {response.status_code}")

            if response.status_code == 200:
                print(f"    ✅ Chat completions working!")
                results['chat_endpoint'] = 'valid'
            elif response.status_code == 401:
                print(f"    ⚠️  Authentication required")
                results['chat_endpoint'] = 'auth_required'
            elif response.status_code == 404:
                print(f"    ❌ Chat endpoint not found")
                results['chat_endpoint'] = 'not_found'
            else:
                print(f"    ⚠️  Status {response.status_code}")
                try:
                    print(f"    Response: {response.json()}")
                except:
                    print(f"    Response: {response.text[:200]}")
                results['chat_endpoint'] = f'status_{response.status_code}'

        except Exception as e:
            print(f"    ❌ Error: {str(e)}")
            results['chat_endpoint'] = 'error'

    return results

def main():
    print("="*60)
    print("Container Server Endpoint Discovery")
    print("="*60)
    print(f"Base Host: {BASE_HOST}")
    print(f"Testing {len(ENDPOINT_PATTERNS)} endpoint patterns...")

    all_results = []
    for pattern in ENDPOINT_PATTERNS:
        result = test_endpoint(pattern)
        all_results.append(result)

    # Summary
    print("\n" + "="*60)
    print("Summary of Results")
    print("="*60)

    working_endpoints = [r for r in all_results if r['working']]

    if working_endpoints:
        print("\n✅ Found working endpoints:\n")
        for result in working_endpoints:
            print(f"  Base URL: {result['base_url']}")
            print(f"    Models: {result['models_endpoint']}")
            print(f"    Chat: {result['chat_endpoint']}")
            print()
    else:
        print("\n❌ No working OpenAI-compatible endpoints found")
        print("\nAll tested endpoints:")
        for result in all_results:
            print(f"  {result['base_url']}")
            print(f"    Models: {result['models_endpoint']}")
            if result['chat_endpoint']:
                print(f"    Chat: {result['chat_endpoint']}")

    print("="*60)

if __name__ == "__main__":
    main()
