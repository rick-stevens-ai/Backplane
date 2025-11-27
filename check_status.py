#!/usr/bin/env python3
"""
Quick status check for computational chemistry testing environment
"""
import sys
import subprocess
from pathlib import Path

def check_services():
    """Check if required services are running"""
    print("\n" + "="*70)
    print("SERVICE STATUS CHECK")
    print("="*70)

    # Check Redis
    try:
        result = subprocess.run(['redis-cli', 'ping'], capture_output=True, text=True, timeout=2)
        if result.returncode == 0 and 'PONG' in result.stdout:
            print("✓ Redis is running")
        else:
            print("✗ Redis is not responding")
            print("  Start with: brew services start redis")
    except Exception as e:
        print(f"✗ Redis check failed: {e}")
        print("  Start with: brew services start redis")

    # Check FastAPI
    try:
        import requests
        response = requests.get("http://127.0.0.1:8000/docs", timeout=2)
        if response.status_code == 200:
            print("✓ FastAPI server is running (http://127.0.0.1:8000)")
        else:
            print("✗ FastAPI server returned error")
    except Exception as e:
        print("✗ FastAPI server is not running")
        print("  Start with: uvicorn main:app --reload")

    # Check Celery (by checking if we can import and connect)
    try:
        from celery import Celery
        test_app = Celery("test", broker="redis://localhost:6379/0")
        # Try to ping
        result = test_app.control.inspect().ping()
        if result:
            print(f"✓ Celery worker is running ({len(result)} worker(s))")
        else:
            print("✗ No Celery workers responding")
            print("  Start with: celery -A tasks.celery_app worker --loglevel=info")
    except Exception as e:
        print("✗ Celery check failed")
        print("  Start with: celery -A tasks.celery_app worker --loglevel=info")


def check_applications():
    """Check if applications are installed"""
    print("\n" + "="*70)
    print("APPLICATION INSTALLATION STATUS")
    print("="*70)

    apps_dir = Path("/Users/stevens/Dropbox/Backplane/APPS")

    if not apps_dir.exists():
        print("✗ APPS directory not found")
        return

    apps = {
        "Quantum ESPRESSO": ["q-e", "pw.x"],
        "CP2K": ["cp2k", "cp2k.ssmp"],
        "GPAW": ["gpaw", "gpaw"],
        "LAMMPS": ["lammps", "lmp"],
        "GROMACS": ["gromacs", "gmx"]
    }

    for app_name, (dir_name, executable) in apps.items():
        app_path = apps_dir / dir_name
        if app_path.exists():
            # Check if there are files
            file_count = len(list(app_path.rglob("*")))
            if file_count > 5:
                print(f"✓ {app_name:20s} installed ({file_count} files)")
            else:
                print(f"⚠ {app_name:20s} directory exists but may be incomplete ({file_count} files)")
        else:
            print(f"✗ {app_name:20s} not found in APPS/")

    # Also check system PATH
    print("\nChecking system PATH for executables:")
    for app_name, (dir_name, executable) in apps.items():
        try:
            result = subprocess.run(['which', executable], capture_output=True, text=True)
            if result.returncode == 0:
                print(f"  ✓ {executable} found: {result.stdout.strip()}")
        except:
            pass


def check_python_deps():
    """Check Python dependencies"""
    print("\n" + "="*70)
    print("PYTHON DEPENDENCIES")
    print("="*70)

    deps = [
        "fastapi",
        "celery",
        "redis",
        "requests",
        "yaml",
        "openai",
        "pydantic",
        "docx"
    ]

    for dep in deps:
        try:
            if dep == "yaml":
                __import__("yaml")
            elif dep == "docx":
                __import__("docx")
            else:
                __import__(dep)
            print(f"✓ {dep:15s} installed")
        except ImportError:
            print(f"✗ {dep:15s} NOT installed")
            if dep == "docx":
                print(f"  Install with: pip install python-docx")
            else:
                print(f"  Install with: pip install {dep}")


def check_model_access():
    """Check if we can access gpt-oss:120b"""
    print("\n" + "="*70)
    print("MODEL ACCESS CHECK")
    print("="*70)

    try:
        import yaml
        from openai import OpenAI

        # Load config
        with open("spark_servers.yaml", 'r') as f:
            config = yaml.safe_load(f)

        # Find spark-container-03
        server_config = None
        for server in config['servers']:
            if server['server'] == "spark-container-03":
                server_config = server
                break

        if not server_config:
            print("✗ spark-container-03 not found in spark_servers.yaml")
            return

        print(f"  Model: {server_config['openai_model']}")
        print(f"  API Base: {server_config['openai_api_base']}")

        # Try to connect
        client = OpenAI(
            api_key=server_config["openai_api_key"],
            base_url=server_config["openai_api_base"]
        )

        # Test with a simple request (with short timeout)
        print("  Testing connection...")
        response = client.chat.completions.create(
            model=server_config["openai_model"],
            messages=[{"role": "user", "content": "Say 'test successful'"}],
            max_tokens=10,
            timeout=10
        )

        print("✓ Successfully connected to gpt-oss:120b")
        print(f"  Response: {response.choices[0].message.content}")

    except Exception as e:
        print(f"✗ Model access failed: {e}")
        print("  Check network connection and server status")


def main():
    """Run all checks"""
    print("\n" + "="*70)
    print("COMPUTATIONAL CHEMISTRY TESTING ENVIRONMENT STATUS")
    print("="*70)

    check_services()
    check_applications()
    check_python_deps()
    check_model_access()

    print("\n" + "="*70)
    print("STATUS CHECK COMPLETE")
    print("="*70)

    print("""
Next steps:
1. If services are down, start them:
   - Terminal 1: brew services start redis
   - Terminal 2: celery -A tasks.celery_app worker --loglevel=info
   - Terminal 3: uvicorn main:app --reload

2. Once all services are running, run tests:
   python test_caffeine_energy.py

3. For more information:
   cat TESTING_GUIDE.md
""")


if __name__ == "__main__":
    main()
