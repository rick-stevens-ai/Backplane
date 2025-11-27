#!/usr/bin/env python3
"""Detailed authentication testing to understand what the server expects"""

import requests
import json

BASE_HOST = "http://100.94.58.120:12000"

def test_endpoint_detailed(base_url: str, auth_variants: list):
    """Test an endpoint with various authentication approaches"""
    print(f"\n{'='*60}")
    print(f"Testing: {base_url}/models")
    print('='*60)

    for auth_name, headers in auth_variants:
        print(f"\n  Trying: {auth_name}")
        try:
            response = requests.get(f"{base_url}/models", headers=headers, timeout=5)
            print(f"    Status: {response.status_code}")

            # Try to get detailed error information
            if response.status_code == 401:
                print(f"    Headers: {dict(response.headers)}")
                try:
                    error_data = response.json()
                    print(f"    Error: {json.dumps(error_data, indent=6)}")
                except:
                    print(f"    Body: {response.text[:300]}")

            elif response.status_code == 200:
                print(f"    ✅ SUCCESS!")
                try:
                    data = response.json()
                    print(f"    Response keys: {list(data.keys())}")
                    if 'data' in data:
                        print(f"    Models: {[m.get('id') for m in data['data'][:3]]}")
                    return True
                except:
                    print(f"    Body: {response.text[:200]}")

            else:
                print(f"    Body: {response.text[:200]}")

        except Exception as e:
            print(f"    Error: {str(e)}")

    return False

def main():
    print("="*60)
    print("Detailed Authentication Testing")
    print("="*60)

    # Different authentication approaches to try
    auth_variants = [
        ("No authentication", {}),
        ("Bearer CELS", {"Authorization": "Bearer CELS"}),
        ("Bearer sk-CELS", {"Authorization": "Bearer sk-CELS"}),
        ("API Key header", {"api-key": "CELS"}),
        ("X-API-Key header", {"X-API-Key": "CELS"}),
        ("Basic auth", {"Authorization": "Basic CELS"}),
    ]

    # Test the most promising endpoints
    endpoints = [
        f"{BASE_HOST}/api/v1",
        f"{BASE_HOST}/api",
        f"{BASE_HOST}/openai/v1",
    ]

    print(f"\nTesting {len(endpoints)} endpoints with {len(auth_variants)} auth methods...")

    for endpoint in endpoints:
        found = test_endpoint_detailed(endpoint, auth_variants)
        if found:
            print(f"\n🎉 Found working configuration!")
            print(f"   Endpoint: {endpoint}")
            break

    # Also check if we can get any info from the base URL
    print(f"\n{'='*60}")
    print("Checking base server info")
    print('='*60)
    try:
        response = requests.get(BASE_HOST, timeout=5)
        print(f"Status: {response.status_code}")
        print(f"Server: {response.headers.get('Server', 'Unknown')}")
        print(f"Content-Type: {response.headers.get('Content-Type', 'Unknown')}")

        # Check if there's an API documentation link
        if 'html' in response.headers.get('Content-Type', ''):
            text = response.text.lower()
            if 'docs' in text or 'api' in text or 'swagger' in text or 'openapi' in text:
                print("\n💡 Server appears to have documentation")
                print("   Check the web interface for API key generation")

    except Exception as e:
        print(f"Error: {e}")

    print("\n" + "="*60)
    print("Recommendations")
    print("="*60)
    print("\n1. This appears to be an OpenWebUI or similar interface")
    print("2. You may need to:")
    print("   - Log into the web UI at http://100.94.58.120:12000")
    print("   - Generate an API key in the settings/profile")
    print("   - Use that API key instead of 'CELS'")
    print("3. The URL /c/{uuid} is typically a chat session, not an API endpoint")
    print("="*60)

if __name__ == "__main__":
    main()
