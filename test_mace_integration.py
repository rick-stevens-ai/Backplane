#!/usr/bin/env python3
"""
Test MACE Integration with Agent - Phase 1 Verification
Tests that agent can call MACE tools and get ML predictions
"""
import sys
import time
sys.path.insert(0, '/Users/stevens/Dropbox/Backplane')

from agent_apps import ComputationalChemistryAgent


def test_single_energy_prediction():
    """Test 1: Single molecule MACE energy prediction"""
    print("\n" + "="*80)
    print("TEST 1: Single Molecule MACE Energy Prediction")
    print("="*80)

    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    request = """Use MACE to quickly predict the energy of ammonia (NH3, SMILES: N).
This is a fast ML prediction for screening purposes."""

    start = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed = time.time() - start

    print(f"\n{'='*80}")
    print(f"RESULT:")
    print(f"{'='*80}")
    print(f"Time elapsed: {elapsed:.2f}s")
    print(f"Expected: <1s (MACE is fast!)")

    if result:
        print(f"\n✓ Test 1 PASSED: Agent successfully used MACE")
        return True
    else:
        print(f"\n✗ Test 1 FAILED")
        return False


def test_rapid_screening():
    """Test 2: MACE rapid screening of 5 molecules"""
    print("\n" + "="*80)
    print("TEST 2: MACE Rapid Screening (5 molecules)")
    print("="*80)

    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    request = """Use MACE to rapidly screen these 5 small nitrogen-containing molecules:
1. Ammonia: N
2. Methylamine: CN
3. Dimethylamine: CNC
4. Trimethylamine: CN(C)C
5. Hydrazine: NN

Rank them by predicted energy from MACE ML screening."""

    start = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed = time.time() - start

    print(f"\n{'='*80}")
    print(f"RESULT:")
    print(f"{'='*80}")
    print(f"Time elapsed: {elapsed:.2f}s")
    print(f"Expected: <5s (MACE batch screening)")

    if result:
        print(f"\n✓ Test 2 PASSED: Agent successfully used MACE batch screening")
        return True
    else:
        print(f"\n✗ Test 2 FAILED")
        return False


def test_hybrid_workflow():
    """Test 3: Hybrid MACE + DFT workflow"""
    print("\n" + "="*80)
    print("TEST 3: Hybrid MACE → DFT Workflow")
    print("="*80)

    agent = ComputationalChemistryAgent(
        server_config_path="spark_servers.yaml",
        server_name="spark-container-03"
    )

    request = """Screen these 3 molecules with MACE first, then validate
the best candidate with GPAW DFT:
1. Ammonia: N
2. Methylamine: CN
3. Ethylamine: CCN

Use MACE for rapid screening, then run GPAW (h=0.2, PBE) on the top candidate."""

    start = time.time()
    result = agent.run_agentic_workflow(request)
    elapsed = time.time() - start

    print(f"\n{'='*80}")
    print(f"RESULT:")
    print(f"{'='*80}")
    print(f"Time elapsed: {elapsed:.2f}s")
    print(f"Expected: ~5s (MACE) + ~2min (GPAW) = ~2-3 minutes")

    if result:
        print(f"\n✓ Test 3 PASSED: Agent successfully used MACE + DFT hybrid workflow")
        return True
    else:
        print(f"\n✗ Test 3 FAILED")
        return False


def main():
    """Run all MACE integration tests"""
    print("\n" + "="*80)
    print("MACE INTEGRATION TEST SUITE")
    print("Phase 1: Verify agent can call MACE tools")
    print("="*80)

    tests = [
        ("Single Energy Prediction", test_single_energy_prediction),
        ("Rapid Screening (5 molecules)", test_rapid_screening),
        ("Hybrid MACE + DFT Workflow", test_hybrid_workflow)
    ]

    results = []
    for test_name, test_func in tests:
        try:
            print(f"\n\nRunning: {test_name}...")
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"\n✗ Test CRASHED: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))

    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    passed = sum(1 for _, success in results if success)
    total = len(results)

    for test_name, success in results:
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{status}: {test_name}")

    print(f"\n{passed}/{total} tests passed")

    if passed == total:
        print("\n✓✓✓ ALL TESTS PASSED! MACE integration working! ✓✓✓")
        return True
    else:
        print(f"\n✗✗✗ {total - passed} test(s) failed ✗✗✗")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
