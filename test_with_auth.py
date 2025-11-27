#!/usr/bin/env python3
"""Test endpoints with authentication to find the working API"""

import requests
import json

BASE_HOST = "http://100.94.58.120:12000"

# Endpoints that returned 401 (need auth)
AUTH_REQUIRED_ENDPOINTS = [
    f"{BASE_HOST}/api/v1",
    f"{BASE_HOST}/api",
    f"{BASE_HOST}/openai/v1",
    f"{BASE_HOST}/ollama/v1",
]

def test_with_auth(base_url: str) -> dict:
    """Test endpoint with authentication"""
    results = {
        'base_url': base_url,
        'models_works': False,
        'chat_works': False
    }

    print(f"\n{'='*60}")
    print(f"Testing: {base_url}")
    print('='*60)

    headers = {
        "Authorization": "Bearer CELS",
        "Content-Type": "application/json"
    }

    # Test /models endpoint with auth
    try:
        models_url = f"{base_url}/models"
        print(f"  GET {models_url}")
        response = requests.get(models_url, headers=headers, timeout=5)
        print(f"    Status: {response.status_code}")

        if response.status_code == 200:
            try:
                data = response.json()
                print(f"    ✅ Models endpoint working!")
                print(f"    Response keys: {list(data.keys())}")

                if 'data' in data:
                    print(f"    Found {len(data['data'])} models:")
                    for model in data['data'][:5]:
                        print(f"      - {model.get('id', 'unknown')}")
                    results['models_works'] = True
                elif 'models' in data:
                    print(f"    Models: {data['models'][:5]}")
                    results['models_works'] = True
                else:
                    print(f"    Full response: {json.dumps(data, indent=2)[:500]}")

            except json.JSONDecodeError:
                print(f"    ⚠️  Response is not JSON")
                print(f"    {response.text[:200]}")
        elif response.status_code == 401:
            print(f"    ❌ Still 401 - authentication failed")
        elif response.status_code == 403:
            print(f"    ❌ 403 Forbidden - may need different credentials")
        else:
            print(f"    ⚠️  Status {response.status_code}")
            print(f"    {response.text[:200]}")

    except Exception as e:
        print(f"    ❌ Error: {str(e)}")

    # If models worked, test chat completions
    if results['models_works']:
        try:
            chat_url = f"{base_url}/chat/completions"
            print(f"\n  POST {chat_url}")

            payload = {
                "model": "gpt-oss:20b",
                "messages": [{"role": "user", "content": "Say 'API works!' and nothing else."}],
                "max_tokens": 10,
                "temperature": 0.1
            }

            print(f"    Model: {payload['model']}")
            response = requests.post(chat_url, headers=headers, json=payload, timeout=30)
            print(f"    Status: {response.status_code}")

            if response.status_code == 200:
                try:
                    data = response.json()
                    print(f"    ✅ Chat completions working!")
                    if 'choices' in data and len(data['choices']) > 0:
                        message = data['choices'][0].get('message', {}).get('content', '')
                        print(f"    Model response: {message}")
                        results['chat_works'] = True
                    else:
                        print(f"    Response: {json.dumps(data, indent=2)[:300]}")
                except:
                    print(f"    Response: {response.text[:200]}")
            else:
                print(f"    ⚠️  Status {response.status_code}")
                print(f"    {response.text[:300]}")

        except requests.exceptions.Timeout:
            print(f"    ⚠️  Request timed out (model may be loading)")
        except Exception as e:
            print(f"    ❌ Error: {str(e)}")

    return results

def main():
    print("="*60)
    print("Testing Endpoints with Authentication")
    print("="*60)
    print(f"Base Host: {BASE_HOST}")
    print(f"Auth Token: Bearer CELS")
    print(f"\nTesting {len(AUTH_REQUIRED_ENDPOINTS)} endpoints...")

    all_results = []
    for endpoint in AUTH_REQUIRED_ENDPOINTS:
        result = test_with_auth(endpoint)
        all_results.append(result)

    # Summary
    print("\n" + "="*60)
    print("Summary")
    print("="*60)

    working = [r for r in all_results if r['models_works'] or r['chat_works']]

    if working:
        print("\n✅ Working endpoints found:\n")
        for result in working:
            print(f"  {result['base_url']}")
            print(f"    Models API: {'✅' if result['models_works'] else '❌'}")
            print(f"    Chat API: {'✅' if result['chat_works'] else '❌'}")
            if result['models_works'] and result['chat_works']:
                print(f"    👉 Use this in spark_servers.yaml")
            print()
    else:
        print("\n❌ No working endpoints found with current authentication")
        print("\nPossible issues:")
        print("  - API key 'CELS' may not be valid for this server")
        print("  - Server may require a different authentication method")
        print("  - This might be a web UI, not an API server")

    print("="*60)

if __name__ == "__main__":
    main()
