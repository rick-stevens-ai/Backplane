# Model Compatibility Notes

This document tracks which models work with the agentic workflow and which don't.

## Requirements

The agentic workflow requires models that support:
1. **Function Calling / Tool Use** - Essential for calling simulation APIs
2. **OpenAI-compatible API** - For standardized integration
3. **JSON argument parsing** - For structured tool inputs

## Compatible Models (Tested & Working)

### ✓ gpt-oss:20b
- **Parameters**: 20.9B
- **Quantization**: MXFP4
- **Function Calling**: Yes
- **Server**: Local Ollama (localhost:11434)
- **Status**: Fully compatible
- **Test Result**: SUCCESS
- **Notes**: Excellent performance, formal output style

### ✓ qwen3:32b
- **Parameters**: 32.8B
- **Quantization**: Q4_K_M
- **Function Calling**: Yes
- **Server**: Local Ollama (localhost:11434)
- **Status**: Fully compatible
- **Test Result**: SUCCESS
- **Notes**: Excellent performance, conversational output style

### ✓ gpt-oss:120b
- **Parameters**: 120B
- **Quantization**: Unknown
- **Function Calling**: Yes
- **Server**: Remote Containerized (100.94.58.120:12000)
- **Status**: Fully compatible
- **Test Result**: SUCCESS
- **Notes**: Outstanding performance, most sophisticated output with interpretation and scientific recommendations

## Incompatible Models

### ✗ llama3:70b
- **Parameters**: 70.6B
- **Quantization**: Q4_0
- **Function Calling**: No
- **Status**: Incompatible
- **Error**: `registry.ollama.ai/library/llama3:70b does not support tools`
- **Notes**: Model available locally but lacks function calling support required for agentic workflows

### ✗ llama3:latest (8B)
- **Parameters**: 8.0B
- **Quantization**: Q4_0
- **Function Calling**: No (assumed, not tested)
- **Status**: Likely incompatible
- **Notes**: Same family as llama3:70b, likely lacks tool support

## Untested Models (Available Locally)

### deepseek-r1:70b
- **Parameters**: 70.6B
- **Quantization**: Q4_K_M
- **Function Calling**: Unknown
- **Status**: Available but untested
- **Notes**: Large model, function calling support unknown

## Recommendations

For agentic workflows with scientific simulations:

1. **Best Practice**: Use models explicitly advertised as supporting function calling/tool use
2. **Tested Options**: gpt-oss:20b, qwen3:32b
3. **Model Size**: Both 20B and 32B models perform excellently
4. **Avoid**: llama3 family models (no function calling support)

## Testing New Models

To test if a model supports the agentic workflow:

```bash
python test_agent_local.py  # For gpt-oss models
python test_qwen3.py        # For qwen3 models
python test_llama3_70b.py   # Example test script
```

Expected error if incompatible:
```
openai.BadRequestError: Error code: 400 - {'error': {'message': 'registry.ollama.ai/library/MODEL_NAME does not support tools', ...}}
```

## Function Calling Test

Quick test to check if a model supports function calling:

```python
from openai import OpenAI

client = OpenAI(
    api_key="ollama",
    base_url="http://localhost:11434/v1"
)

try:
    response = client.chat.completions.create(
        model="YOUR_MODEL_NAME",
        messages=[{"role": "user", "content": "Test"}],
        tools=[{
            "type": "function",
            "function": {
                "name": "test_function",
                "parameters": {"type": "object", "properties": {}}
            }
        }]
    )
    print("✓ Model supports function calling")
except Exception as e:
    print(f"✗ Model does NOT support function calling: {e}")
```

## Last Updated

November 21, 2025
