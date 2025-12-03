#!/usr/bin/env python3
"""
Debug script to inspect the oss120b API response structure
and identify where reasoning/thinking content is stored.
"""
import sys
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from openai import OpenAI
import json

# Initialize client (using the same API key as spark_servers.yaml)
client = OpenAI(
    base_url="http://100.94.58.120:12000/api/v1",
    api_key="sk-aad91216be8b44a194672560b4dde7b4"
)

print("="*80)
print("OSS120B REASONING MODEL - API RESPONSE INSPECTION")
print("="*80)
print()

# Simple test request
messages = [
    {
        "role": "system",
        "content": "You are a helpful assistant. Think step by step and show your reasoning."
    },
    {
        "role": "user",
        "content": "What is 15 * 23? Please think through this step by step."
    }
]

print("Sending test request to oss120b...")
print("Model: gpt-oss:120b")
print()

try:
    response = client.chat.completions.create(
        model="gpt-oss:120b",
        messages=messages
    )

    print("="*80)
    print("RESPONSE OBJECT STRUCTURE")
    print("="*80)
    print()

    # Inspect the response object
    message = response.choices[0].message

    print("Available attributes in response:")
    print(f"  - response type: {type(response)}")
    print(f"  - response attributes: {dir(response)}")
    print()

    print("Available attributes in response.choices[0]:")
    choice = response.choices[0]
    print(f"  - choice type: {type(choice)}")
    print(f"  - choice attributes: {dir(choice)}")
    print()

    print("Available attributes in message:")
    print(f"  - message type: {type(message)}")
    print(f"  - message attributes: {dir(message)}")
    print()

    # Check for reasoning-related fields
    print("="*80)
    print("CHECKING FOR REASONING FIELDS")
    print("="*80)
    print()

    reasoning_fields = [
        'reasoning', 'thinking', 'reasoning_content', 'thought',
        'internal_reasoning', 'chain_of_thought', 'thoughts'
    ]

    print("Checking message object:")
    for field in reasoning_fields:
        if hasattr(message, field):
            value = getattr(message, field)
            print(f"  ✓ Found: message.{field}")
            print(f"    Type: {type(value)}")
            print(f"    Value: {value[:200] if isinstance(value, str) and len(value) > 200 else value}")
            print()

    print("Checking choice object:")
    for field in reasoning_fields:
        if hasattr(choice, field):
            value = getattr(choice, field)
            print(f"  ✓ Found: choice.{field}")
            print(f"    Type: {type(value)}")
            print(f"    Value: {value[:200] if isinstance(value, str) and len(value) > 200 else value}")
            print()

    # Try to get raw dict representation
    print("="*80)
    print("RAW RESPONSE DATA")
    print("="*80)
    print()

    # Try model_dump if it's a Pydantic model
    if hasattr(response, 'model_dump'):
        print("Response as dict (via model_dump):")
        response_dict = response.model_dump()
        print(json.dumps(response_dict, indent=2))
    elif hasattr(response, 'dict'):
        print("Response as dict (via dict):")
        response_dict = response.dict()
        print(json.dumps(response_dict, indent=2))
    else:
        print("Could not convert response to dict")
        print(f"Response: {response}")

    print()
    print("="*80)
    print("MESSAGE CONTENT")
    print("="*80)
    print()
    print(message.content)

except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

print()
print("="*80)
print("INSPECTION COMPLETE")
print("="*80)
