# MACE Integration with Scientific AI Backplane

## Executive Summary

**MACE** (Machine learning Atomic Cluster Expansion) is a foundation model for molecular calculations that provides **fast ML predictions** (~0.3s/molecule) to complement the existing **high-accuracy DFT simulations** (2-5 min/molecule) in the Backplane workflow.

## Current State Analysis

### Backplane Workflow (Existing)
1. **gpt-oss:120b** → Suggests molecules, interprets results
2. **DFT Simulations** → High-accuracy calculations
   - Quantum ESPRESSO (~7 min)
   - CP2K (~3 min)
   - GPAW (~2 min)
   - LAMMPS/GROMACS (~3 min)
3. **Bottleneck**: Screening many molecules is slow (5 molecules = 15-30 minutes)

### MACE Capabilities (Available)
- **Location**: `~/Dropbox/MACE/`
- **MCP Server**: Already implemented (`mace_mcp_server.py`)
- **Speed**: ~0.3s per molecule (100x faster than DFT!)
- **Accuracy**: Near-quantum accuracy (trained on DFT data)
- **Models**:
  - **MACE-OFF**: H, C, N, O, P, S, F, Cl, Br, I (perfect for organic catalysts!)
  - **MACE-MP**: 89 elements (for materials/metal clusters)
- **5 Tools**:
  1. `calculate_energy` - Fast energy predictions
  2. `calculate_forces` - Atomic forces
  3. `predict_properties` - Energy + forces together
  4. `batch_calculate_energies` - Process multiple molecules
  5. `optimize_geometry` - ML-based geometry optimization

## Integration Strategy

### Proposed Workflow Enhancement

```
User Request: "Screen 100 NH3 catalyst candidates"
                    ↓
         gpt-oss:120b (suggest molecules)
                    ↓
    ┌───────────────┴───────────────┐
    │                               │
MACE Rapid Screening            (optional)
  • 100 molecules                Direct to DFT
  • ~30 seconds total            if <10 molecules
  • Predict energies
  • Rank by ML predictions
    ↓
Top 10 candidates selected
    │
    └───────────────┬─────────────────┘
                    ↓
        DFT Validation (GPAW/CP2K/QE)
          • 10 molecules
          • ~20-50 minutes
          • High-accuracy results
                    ↓
        Final ranked results
```

**Key Benefit**: 100 molecules screened in ~30 seconds + validation of top 10 in 20 min
**vs. Previous**: 100 molecules × 3 min = 5 hours!

## Implementation Options

### Option 1: MACE as Agent Tool (RECOMMENDED)
**Approach**: Add MACE to gpt-oss:120b's tool set via MCP

**Pros:**
- Agent decides when to use MACE vs DFT intelligently
- No job queue overhead (MACE is fast, runs synchronously)
- Fits agentic workflow perfectly
- Already MCP-ready

**Implementation**:
1. Add MACE client to `agent_apps.py`
2. Define tool functions for agent to call
3. Agent orchestrates: MACE screening → DFT validation

### Option 2: MACE as Separate Wrapper
**Approach**: Create `mace_wrapper.py` like other simulation codes

**Pros:**
- Consistent with existing architecture
- Could use job queue for batch processing

**Cons:**
- Overhead not needed (MACE is instant)
- Doesn't fit MACE's strength (rapid screening)

### Option 3: Hybrid Approach
**Approach**: Both MCP access AND wrapper

**Use Case**: Direct Python scripts can use wrapper, agent uses MCP

**Complexity**: Higher, may be overkill

## Recommended Implementation (Option 1)

### Architecture

```python
# agent_apps.py (enhanced)

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

class ComputationalChemistryAgent:
    def __init__(self, ...):
        self.mace_server_path = "/Users/stevens/Dropbox/MACE/mace_mcp_server.py"
        self.tools = [
            # Existing DFT tools
            {
                "name": "run_quantum_espresso",
                "description": "High-accuracy plane-wave DFT...",
                ...
            },
            # ... other DFT tools ...

            # NEW MACE tools
            {
                "name": "mace_rapid_screening",
                "description": "Fast ML screening of multiple molecules (~0.3s each). Use for large-scale screening before DFT validation.",
                "parameters": {
                    "molecules": "List of SMILES or XYZ structures",
                    "model_type": "'off' for organic, 'mp' for materials"
                },
                "function": self.mace_rapid_screening
            },
            {
                "name": "mace_predict_energy",
                "description": "Fast ML energy prediction for single molecule. Use for quick estimates.",
                "parameters": {
                    "smiles_or_xyz": "Molecule structure",
                    "model_type": "'off' or 'mp'"
                },
                "function": self.mace_predict_energy
            },
            {
                "name": "mace_optimize_structure",
                "description": "ML-based geometry optimization (~2s). Faster than DFT but less accurate.",
                "parameters": {
                    "smiles_or_xyz": "Initial structure",
                    "model_type": "'off' or 'mp'"
                },
                "function": self.mace_optimize_structure
            }
        ]

    async def mace_predict_energy(self, smiles_or_xyz, model_type="off"):
        """Fast ML energy prediction via MCP."""
        server_params = StdioServerParameters(
            command="python3.10",
            args=[self.mace_server_path]
        )

        async with stdio_client(server_params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()

                result = await session.call_tool(
                    "calculate_energy",
                    arguments={
                        "atoms_data": self._convert_to_xyz(smiles_or_xyz),
                        "format": "xyz",
                        "model_type": model_type
                    }
                )

                return self._parse_mace_result(result)

    async def mace_rapid_screening(self, molecules, model_type="off"):
        """Batch screening with MACE."""
        server_params = StdioServerParameters(
            command="python3.10",
            args=[self.mace_server_path]
        )

        async with stdio_client(server_params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()

                # Convert to structures
                structures = [
                    {
                        "atoms_data": self._convert_to_xyz(mol),
                        "format": "xyz"
                    }
                    for mol in molecules
                ]

                result = await session.call_tool(
                    "batch_calculate_energies",
                    arguments={
                        "structures": structures,
                        "model_type": model_type
                    }
                )

                return self._parse_mace_batch_result(result)
```

### Integration Points

1. **Fast Pre-Screening**
   - User: "Screen 50 NH3 catalysts"
   - Agent: Uses MACE to rank → Top 10 to DFT

2. **Structure Preparation**
   - Agent: Use MACE to optimize geometry
   - Then: Use optimized structure for DFT

3. **Property Estimation**
   - Agent: Quick MACE estimate to check if molecule is viable
   - If promising: Run full DFT

4. **Iterative Design**
   - Agent: Generate variant → MACE check → If good → DFT validate
   - Rapid design-test cycles

## Implementation Files

### New Files to Create

1. **`mace_client.py`** - Helper module for MACE MCP calls
   ```python
   """Helper module for calling MACE MCP server."""
   class MACEClient:
       async def predict_energy(self, structure, model="off"):
           ...
       async def batch_predict(self, structures, model="off"):
           ...
       async def optimize_geometry(self, structure, model="off"):
           ...
   ```

2. **`agent_apps.py`** (modifications) - Add MACE tools to agent

3. **`test_mace_integration.py`** - Integration tests
   ```python
   """Test MACE integration with Backplane."""
   # Test 1: Agent uses MACE for screening
   # Test 2: MACE + DFT workflow
   # Test 3: Rapid screening of 20 molecules
   ```

### Modified Files

1. **`agent_apps.py`**
   - Import MACEClient
   - Add 3 MACE tool definitions
   - Add helper methods for SMILES→XYZ conversion

2. **`requirements.txt`** (if needed)
   - Ensure MCP SDK is included

## Use Case Examples

### Use Case 1: Large-Scale Catalyst Screening

```
User: "Screen 100 potential NH3 catalysts and validate the top 5 with DFT"

Agent workflow:
1. Generate 100 catalyst SMILES from knowledge base
2. Call mace_rapid_screening(molecules=100, model="off")
   → 30 seconds, ranked by ML energy
3. Select top 5 candidates
4. For each top 5:
   - Call run_gpaw(smiles=candidate)
   - Validate with high-accuracy DFT
5. Return final ranked list

Time: 30s (MACE) + 15min (5× DFT) = ~16 minutes
vs. Previous: 100× 3min = 5 hours!
```

### Use Case 2: Structure Optimization + Validation

```
User: "Optimize this distorted pyridine structure and calculate its accurate energy"

Agent workflow:
1. mace_optimize_structure(distorted_pyridine) → 2s
2. run_gpaw(optimized_structure) → 3min
3. Return optimized structure + accurate DFT energy

Time: 2s + 3min = ~3 minutes
Benefit: Good starting structure for DFT (converges faster)
```

### Use Case 3: Rapid Property Estimates

```
User: "Estimate energies for these 10 molecules quickly"

Agent workflow:
1. Recognize "quickly" → Use MACE not DFT
2. mace_rapid_screening(10 molecules)
3. Return ML predictions with confidence

Time: ~3 seconds (vs 30 minutes for DFT)
Use: Quick feasibility checks
```

## Performance Characteristics

| Method | Speed | Accuracy | Use Case |
|--------|-------|----------|----------|
| **MACE** | 0.3s/mol | ~98% DFT | Screening, pre-filtering |
| **GPAW** | 2-3 min/mol | DFT | Final validation, small molecules |
| **CP2K** | 3-5 min/mol | Hybrid DFT | Metal clusters, high accuracy |
| **QE** | 5-10 min/mol | Plane-wave DFT | Validation, convergence studies |

**Optimal Strategy**: MACE screens 100 → DFT validates top 10

## Technical Requirements

### Dependencies (Already Met)
- ✓ Python 3.10
- ✓ MACE-torch installed (`~/Dropbox/MACE/`)
- ✓ MCP SDK installed
- ✓ MCP server tested and working

### Additional Setup Needed
- Minimal: Just import MCP client in `agent_apps.py`
- MACE server runs on-demand (stdio communication)

## Testing Strategy

### Test Suite

1. **Basic Connectivity Test**
   ```bash
   python3 test_mace_integration.py --test basic
   ```
   - Verify MCP connection works
   - Single energy calculation

2. **Screening Performance Test**
   ```bash
   python3 test_mace_integration.py --test screening
   ```
   - 20 molecules via MACE
   - Measure total time
   - Compare to DFT baseline

3. **Agent Integration Test**
   ```bash
   python3 test_mace_integration.py --test agent
   ```
   - Agent uses MACE autonomously
   - Hybrid MACE + DFT workflow
   - Verify correct tool selection

4. **End-to-End Workflow Test**
   ```bash
   python3 test_mace_integration.py --test e2e
   ```
   - Full catalyst screening pipeline
   - MACE pre-screen → DFT validation
   - Generate report

## Rollout Plan

### Phase 1: Basic Integration (1 hour)
- [ ] Create `mace_client.py` helper module
- [ ] Add 3 MACE tools to `agent_apps.py`
- [ ] Test single MACE call from agent
- [ ] Deliverable: Agent can call MACE for energy

### Phase 2: Screening Workflow (30 min)
- [ ] Implement `mace_rapid_screening` batch function
- [ ] Test with 20-molecule screening
- [ ] Deliverable: Batch screening working

### Phase 3: Integration Tests (30 min)
- [ ] Create `test_mace_integration.py`
- [ ] Run all test cases
- [ ] Document performance improvements
- [ ] Deliverable: Full test suite passing

### Phase 4: Documentation & Examples (30 min)
- [ ] Add MACE section to TESTING_GUIDE.md
- [ ] Create example workflow script
- [ ] Update CATALYST_APPLICATIONS.md
- [ ] Deliverable: Complete documentation

**Total Estimated Time**: 2.5 hours to full integration

## Success Criteria

✅ **Functional**:
- Agent can call MACE tools
- Batch screening works (20+ molecules)
- MACE + DFT hybrid workflow demonstrated

✅ **Performance**:
- MACE: <1s per molecule
- Screening 20 molecules: <30s total
- 10x speed improvement demonstrated

✅ **Integration Quality**:
- Tests passing
- Documentation complete
- Example workflows provided

## Risks & Mitigation

### Risk 1: MCP Connection Overhead
**Mitigation**: Profile first call vs subsequent calls. May need session pooling.

### Risk 2: SMILES → XYZ Conversion Issues
**Mitigation**: Already have RDKit integration (from GPAW fix!). Reuse same code.

### Risk 3: Agent Doesn't Use MACE Appropriately
**Mitigation**: Tool descriptions guide agent to use MACE for screening, DFT for validation.

## Future Enhancements

### Near-Term
- Confidence scores on MACE predictions
- Automatic MACE→DFT handoff for low-confidence results
- Property-based filtering (not just energy)

### Long-Term
- Active learning: DFT results retrain MACE
- Uncertainty quantification
- Multi-fidelity optimization

## Conclusion

**MACE integration provides 10-100x speedup for molecular screening** while maintaining the option for high-accuracy DFT validation of promising candidates. The MCP-based integration allows the agent to intelligently choose between fast ML predictions and accurate quantum calculations based on the task requirements.

**Recommended Action**: Proceed with Phase 1 implementation immediately while NH3 screening runs in background.

---

**Next Steps**: Implement `mace_client.py` and integrate into `agent_apps.py`
