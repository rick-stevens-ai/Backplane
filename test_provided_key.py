#!/usr/bin/env python3
"""Test the provided API key with discovered endpoints"""

import requests
import json

BASE_HOST = "http://100.94.58.120:12000"
API_KEY = "sk-aad91216be8b44a194672560b4dde7b4"

# Promising endpoints we discovered
ENDPOINTS = [
    f"{BASE_HOST}/api/v1",
    f"{BASE_HOST}/openai/v1",
    f"{BASE_HOST}/api",
]

def test_endpoint_full(base_url: str, api_key: str):
    """Test an endpoint with the provided API key"""
    print(f"\n{'='*60}")
    print(f"Testing: {base_url}")
    print('='*60)

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }

    results = {
        'base_url': base_url,
        'models_works': False,
        'chat_works': False,
        'models': []
    }

    # Test /models endpoint
    try:
        models_url = f"{base_url}/models"
        print(f"\n1. Testing Models Endpoint")
        print(f"   GET {models_url}")
        response = requests.get(models_url, headers=headers, timeout=10)
        print(f"   Status: {response.status_code}")

        if response.status_code == 200:
            print(f"   ✅ SUCCESS!")
            data = response.json()

            if 'data' in data:
                results['models_works'] = True
                results['models'] = [m.get('id') for m in data['data']]
                print(f"   Found {len(data['data'])} models:")
                for model in data['data'][:10]:
                    model_id = model.get('id', 'unknown')
                    print(f"     - {model_id}")
            else:
                print(f"   Response keys: {list(data.keys())}")
                print(f"   Data: {json.dumps(data, indent=6)[:500]}")
        else:
            print(f"   ❌ Failed with status {response.status_code}")
            try:
                error = response.json()
                print(f"   Error: {json.dumps(error, indent=6)}")
            except:
                print(f"   Response: {response.text[:300]}")
            return results

    except Exception as e:
        print(f"   ❌ Error: {str(e)}")
        return results

    # If models worked, test chat completions
    if results['models_works']:
        try:
            chat_url = f"{base_url}/chat/completions"
            print(f"\n2. Testing Chat Completions Endpoint")
            print(f"   POST {chat_url}")

            # Use first available model
            test_model = results['models'][0] if results['models'] else "gpt-oss:20b"
            print(f"   Using model: {test_model}")

            payload = {
                "model": test_model,
                "messages": [
                    {"role": "user", "content": "Say 'Hello from containerized server!' and nothing else."}
                ],
                "max_tokens": 20,
                "temperature": 0.1
            }

            response = requests.post(chat_url, headers=headers, json=payload, timeout=60)
            print(f"   Status: {response.status_code}")

            if response.status_code == 200:
                print(f"   ✅ SUCCESS!")
                data = response.json()
                if 'choices' in data and len(data['choices']) > 0:
                    message = data['choices'][0].get('message', {}).get('content', '')
                    print(f"   Model response: {message}")
                    results['chat_works'] = True
                else:
                    print(f"   Unexpected response format:")
                    print(f"   {json.dumps(data, indent=6)[:500]}")
            else:
                print(f"   ⚠️  Failed with status {response.status_code}")
                try:
                    error = response.json()
                    print(f"   Error: {json.dumps(error, indent=6)[:300]}")
                except:
                    print(f"   Response: {response.text[:300]}")

        except requests.exceptions.Timeout:
            print(f"   ⚠️  Request timed out (model may be loading/slow)")
        except Exception as e:
            print(f"   ❌ Error: {str(e)}")

    return results

def main():
    print("="*60)
    print("Testing Provided API Key")
    print("="*60)
    print(f"API Key: {API_KEY[:20]}...")
    print(f"Testing {len(ENDPOINTS)} endpoints...")

    all_results = []
    for endpoint in ENDPOINTS:
        result = test_endpoint_full(endpoint, API_KEY)
        all_results.append(result)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    working = [r for r in all_results if r['models_works']]

    if working:
        print("\n✅ Working API endpoints found!\n")
        for result in working:
            print(f"Base URL: {result['base_url']}")
            print(f"  Models API: {'✅' if result['models_works'] else '❌'}")
            print(f"  Chat API: {'✅' if result['chat_works'] else '❌'}")
            print(f"  Available models: {len(result['models'])}")

            if result['models_works'] and result['chat_works']:
                print(f"\n  👉 RECOMMENDED for spark_servers.yaml:")
                print(f"     openai_api_base: \"{result['base_url']}\"")
                print(f"     openai_api_key: \"{API_KEY}\"")
                if result['models']:
                    print(f"     Available models to use:")
                    for model in result['models'][:5]:
                        print(f"       - {model}")
            print()
    else:
        print("\n❌ No working endpoints found")
        print("API key may be invalid or expired")

    print("="*60)

if __name__ == "__main__":
    main()
