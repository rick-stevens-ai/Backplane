#!/usr/bin/env python3
"""Final comprehensive test of spark_servers.yaml configuration"""

import yaml
import requests
import json

def test_server_full(server: dict) -> dict:
    """Comprehensively test a server configuration"""
    results = {
        'server': server['shortname'],
        'config_valid': True,
        'connectivity': False,
        'models_api': False,
        'chat_api': False,
        'error': None
    }

    print(f"\n{'='*60}")
    print(f"Testing: {server['shortname']}")
    print('='*60)
    print(f"  Server ID: {server['server']}")
    print(f"  API Base: {server['openai_api_base']}")
    print(f"  Model: {server['openai_model']}")
    print(f"  API Key: {server['openai_api_key'][:15]}...")

    headers = {
        "Authorization": f"Bearer {server['openai_api_key']}",
        "Content-Type": "application/json"
    }

    # Test 1: Models endpoint
    try:
        models_url = f"{server['openai_api_base']}/models"
        print(f"\n  [1/2] Testing models endpoint...")
        response = requests.get(models_url, headers=headers, timeout=10)

        if response.status_code == 200:
            data = response.json()
            if 'data' in data:
                print(f"        ✅ Models API working ({len(data['data'])} models)")
                results['models_api'] = True
                results['connectivity'] = True
            else:
                print(f"        ⚠️  Unexpected response format")
        elif response.status_code == 401:
            print(f"        ❌ Authentication failed (401)")
            results['error'] = "Authentication failed"
        else:
            print(f"        ⚠️  Status {response.status_code}")
            results['error'] = f"Status {response.status_code}"

    except requests.exceptions.ConnectionError:
        print(f"        ⚠️  Connection refused (server offline)")
        results['error'] = "Connection refused"
    except requests.exceptions.Timeout:
        print(f"        ⚠️  Timeout")
        results['error'] = "Timeout"
    except Exception as e:
        print(f"        ❌ Error: {str(e)}")
        results['error'] = str(e)

    # Test 2: Chat completions (only if models worked)
    if results['models_api']:
        try:
            chat_url = f"{server['openai_api_base']}/chat/completions"
            print(f"  [2/2] Testing chat completions...")

            payload = {
                "model": server['openai_model'],
                "messages": [{"role": "user", "content": "Reply with just 'OK'"}],
                "max_tokens": 10,
                "temperature": 0.1
            }

            response = requests.post(chat_url, headers=headers, json=payload, timeout=30)

            if response.status_code == 200:
                data = response.json()
                if 'choices' in data and len(data['choices']) > 0:
                    message = data['choices'][0].get('message', {}).get('content', '')
                    print(f"        ✅ Chat API working (response: {message[:50]})")
                    results['chat_api'] = True
                else:
                    print(f"        ⚠️  Unexpected response format")
            else:
                print(f"        ⚠️  Status {response.status_code}")

        except requests.exceptions.Timeout:
            print(f"        ⚠️  Timeout (model may be loading)")
        except Exception as e:
            print(f"        ❌ Error: {str(e)}")

    return results

def main():
    print("="*60)
    print("Final Configuration Test - spark_servers.yaml")
    print("="*60)

    # Load configuration
    try:
        with open("spark_servers.yaml", 'r') as f:
            data = yaml.safe_load(f)
        servers = data['servers']
        print(f"\n✅ Loaded {len(servers)} server configurations")
    except Exception as e:
        print(f"❌ Failed to load configuration: {e}")
        return

    # Test each server
    all_results = []
    for server in servers:
        result = test_server_full(server)
        all_results.append(result)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    working = [r for r in all_results if r['models_api'] and r['chat_api']]
    partially_working = [r for r in all_results if r['connectivity'] and not r['chat_api']]
    offline = [r for r in all_results if not r['connectivity']]

    print(f"\n✅ Fully Working: {len(working)}/{len(all_results)}")
    for r in working:
        print(f"   - {r['server']}")

    if partially_working:
        print(f"\n⚠️  Partially Working: {len(partially_working)}")
        for r in partially_working:
            print(f"   - {r['server']} (models OK, chat needs testing)")

    if offline:
        print(f"\n⚠️  Offline/Unreachable: {len(offline)}")
        for r in offline:
            error = f" ({r['error']})" if r['error'] else ""
            print(f"   - {r['server']}{error}")

    print("\n" + "="*60)

    # Show working configurations
    if working:
        print("Ready to use configurations:")
        print("="*60)
        for r in working:
            print(f"\n{r['server']}:")
            server = next(s for s in servers if s['shortname'] == r['server'])
            print(f"  openai_api_base: {server['openai_api_base']}")
            print(f"  openai_model: {server['openai_model']}")
        print("\n" + "="*60)

if __name__ == "__main__":
    main()
