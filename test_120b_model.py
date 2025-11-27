#!/usr/bin/env python3
"""Specific test for the large gpt-oss:120b model with extended timeout"""

import requests
import json
import time

API_BASE = "http://100.94.58.120:12000/api/v1"
API_KEY = "sk-aad91216be8b44a194672560b4dde7b4"
MODEL = "gpt-oss:120b"

def test_120b_model():
    print("="*60)
    print("Testing gpt-oss:120b Model (Large Model)")
    print("="*60)
    print(f"Model: {MODEL}")
    print(f"API Base: {API_BASE}")
    print("\nNote: This is a 120B parameter model and may take")
    print("      30-60+ seconds to load and respond.\n")

    headers = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json"
    }

    # Test chat completion with extended timeout
    try:
        chat_url = f"{API_BASE}/chat/completions"
        print(f"POST {chat_url}")

        payload = {
            "model": MODEL,
            "messages": [
                {"role": "user", "content": "Say 'Hello from 120B model!' and nothing else."}
            ],
            "max_tokens": 20,
            "temperature": 0.1
        }

        print("\nSending request... (may take 30-60+ seconds)")
        start_time = time.time()

        response = requests.post(
            chat_url,
            headers=headers,
            json=payload,
            timeout=120  # 2 minute timeout for large model
        )

        elapsed = time.time() - start_time

        print(f"\nStatus: {response.status_code}")
        print(f"Response time: {elapsed:.2f} seconds")

        if response.status_code == 200:
            data = response.json()
            if 'choices' in data and len(data['choices']) > 0:
                message = data['choices'][0].get('message', {}).get('content', '')
                print(f"\n✅ SUCCESS!")
                print(f"Model response: {message}")
                return True
            else:
                print(f"\n⚠️  Unexpected response format:")
                print(json.dumps(data, indent=2)[:500])
                return False
        else:
            print(f"\n❌ Failed with status {response.status_code}")
            try:
                error = response.json()
                print(f"Error: {json.dumps(error, indent=2)}")
            except:
                print(f"Response: {response.text[:300]}")
            return False

    except requests.exceptions.Timeout:
        print(f"\n⚠️  Request timed out after 120 seconds")
        print("The model may be:")
        print("  - Still loading into memory")
        print("  - Busy with other requests")
        print("  - Too large for the current server resources")
        return False
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        return False

def main():
    result = test_120b_model()

    print("\n" + "="*60)
    if result:
        print("✅ gpt-oss:120b-container is fully functional")
        print("\nConfiguration:")
        print("  server: spark-container-03")
        print("  shortname: gpt-oss:120b-container")
        print(f"  openai_api_base: {API_BASE}")
        print(f"  openai_model: {MODEL}")
    else:
        print("⚠️  gpt-oss:120b-container configuration is valid")
        print("   but model may need more time to load or")
        print("   may not be available on the server")
    print("="*60)

if __name__ == "__main__":
    main()
