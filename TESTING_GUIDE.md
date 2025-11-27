# Testing Guide for Computational Chemistry Applications

## Overview

This guide explains how to test the five integrated computational chemistry applications using the caffeine molecule energy calculation benchmark.

## Test Files

### 1. `test_caffeine_energy.py` - Comprehensive Energy Calculations

**What it does:**
- Calculates the total energy of caffeine molecule (C8H10N4O2) using all five applications
- Uses gpt-oss:120b to autonomously set parameters
- Compares results across quantum (DFT) and classical (force field) methods
- Provides detailed analysis and timing information

**Caffeine molecule:**
- SMILES: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
- 24 atoms: 8 Carbon, 10 Hydrogen, 4 Nitrogen, 2 Oxygen
- Common stimulant drug, excellent test case

**Calculations performed:**
1. **Quantum ESPRESSO**: Plane-wave DFT SCF calculation
2. **CP2K**: Mixed Gaussian/plane-wave DFT energy
3. **GPAW**: Real-space grid DFT energy
4. **LAMMPS**: Classical force field energy minimization
5. **GROMACS**: Biomolecular force field energy minimization

### 2. `test_catalyst_applications.py` - General Application Tests

**What it does:**
- Tests each application with simple requests
- Validates agentic workflow with natural language
- Demonstrates autonomous application selection

## Prerequisites

### 1. Services Running

You need three services running:

**Terminal 1: Redis**
```bash
brew services start redis
# Or check status: brew services list
```

**Terminal 2: Celery Worker**
```bash
cd /Users/stevens/Dropbox/Backplane
celery -A tasks.celery_app worker --loglevel=info
```

**Terminal 3: FastAPI Server**
```bash
cd /Users/stevens/Dropbox/Backplane
uvicorn main:app --reload
```

### 2. Applications Installed

Check installation status in `APPS/` directory:
```bash
ls -la APPS/
# Should show: q-e/, cp2k/, gpaw/, lammps/, gromacs/
```

### 3. Python Dependencies

```bash
pip install python-docx  # Already installed
# Other dependencies should be installed from earlier setup
```

## Running Tests

### Quick Test (Caffeine Energy Calculations)

```bash
cd /Users/stevens/Dropbox/Backplane
python test_caffeine_energy.py
```

**What to expect:**
- Interactive prompts before starting
- Detailed output for each calculation
- Total time: 5-30 minutes depending on application installation
- Comparative analysis at the end

**Output includes:**
- Energy values in multiple units (eV, Hartree, kcal/mol)
- Convergence status
- Execution timing
- Output file locations
- Cross-method comparison

### Full Application Test

```bash
python test_catalyst_applications.py
```

**What to expect:**
- Tests all five applications with simple molecules
- Natural language requests to gpt-oss:120b
- Validates end-to-end workflow

## Understanding Results

### Energy Values

Different methods report energies in different units:

**Quantum methods (DFT):**
- **Hartree (Ha)**: Atomic units of energy (1 Ha = 27.2114 eV)
- **Electronvolts (eV)**: Common in physics/chemistry
- **Rydberg (Ry)**: 1 Ry = 13.6057 eV = 0.5 Ha

Example: -11,000 eV = -404.6 Ha

**Classical methods (force fields):**
- **kcal/mol**: Kilocalories per mole (1 kcal/mol = 0.0434 eV)
- **kJ/mol**: Kilojoules per mole (1 kJ/mol = 0.0104 eV)

Example: 100 kcal/mol = 4.34 eV

### Expected Behavior

**DFT methods should give:**
- Large negative energies (e.g., -10,000 to -15,000 eV for caffeine)
- Similar values across QE, CP2K, and GPAW (within 1-5%)
- Convergence after 5-20 SCF iterations

**Classical methods should give:**
- Smaller magnitude energies (e.g., -100 to +500 kcal/mol)
- Dependent on force field parameters
- Quick minimization (10-100 steps)

### Success Criteria

✓ **Calculation completes** without errors
✓ **Energy is reported** in appropriate units
✓ **Convergence achieved** (for iterative methods)
✓ **Output files generated** in job directory

## Troubleshooting

### Problem: Application not found

**Error message:**
```
version: "not_installed"
```

**Solution:**
- Check APPS/ directory: `ls APPS/`
- Wait for installation to complete
- Check system PATH: `which pw.x` (for QE), `which cp2k.ssmp`, etc.

### Problem: Celery worker not responding

**Symptoms:**
- Job stays in QUEUED status
- No output from worker terminal

**Solution:**
```bash
# Restart Celery worker
# Stop current worker (Ctrl+C)
celery -A tasks.celery_app worker --loglevel=info --purge
```

### Problem: Import errors

**Error message:**
```
ModuleNotFoundError: No module named 'wrappers'
```

**Solution:**
```bash
cd /Users/stevens/Dropbox/Backplane
export PYTHONPATH="${PYTHONPATH}:/Users/stevens/Dropbox/Backplane"
python test_caffeine_energy.py
```

### Problem: Simulation timeout

**Error message:**
```
status: "timeout"
```

**Solution:**
- Increase timeout in wrapper (default: 1 hour)
- Use smaller molecule for initial tests
- Check if application is actually running

### Problem: Convergence failure

**Symptoms:**
```
converged: False
```

**Solution:**
- Normal for first attempts with default parameters
- Agent may adjust parameters on retry
- Check input file for reasonable settings
- May need better initial geometry

## Interpreting Output

### Successful Run Example

```
================================================================================
TEST 1: QUANTUM ESPRESSO - DFT Energy Calculation
================================================================================

Quantum ESPRESSO is a plane-wave DFT code excellent for accurate electronic
structure calculations. We'll perform an SCF calculation to determine the
total energy of caffeine.

Submitting job to gpt-oss:120b agent...

────────────────────────────────────────────────────────────────────────────────
Results for Quantum ESPRESSO:
────────────────────────────────────────────────────────────────────────────────
✓ Calculation completed successfully
  Application: quantum_espresso
  Experiment: Caffeine_Energy_QE
  Execution time: 127.3 seconds

Energy Results:
  Energy: -423.456789 Hartree
        = -11520.34 eV
  Converged: True
  SCF Iterations: 12

Output Files:
  input: /tmp/backplane_jobs/abc123/pw.in
  output: /tmp/backplane_jobs/abc123/pw.out
────────────────────────────────────────────────────────────────────────────────

Total time: 143.2 seconds
```

### Comparative Analysis

At the end, you'll see:

```
================================================================================
                        COMPARATIVE ANALYSIS
================================================================================

Caffeine Molecule: C8H10N4O2 (24 atoms)
SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C

Energy Comparison Across Methods:

Application          Method     Energy                        Converged    Time (s)
──────────────────────────────────────────────────────────────────────────────────
Quantum ESPRESSO     DFT        -423.46 Ha (-11520.34 eV)   True         127.3
CP2K                 DFT        -11502.12 eV                 True         98.7
GPAW                 DFT        -11535.89 eV                 True         65.4
LAMMPS               Classical  -245.67 kcal/mol             True         12.1
GROMACS              Classical  -312.45 kcal/mol             True         8.9
```

## Advanced Usage

### Custom Calculations

You can modify the test script to:

1. **Change molecule**: Edit `CAFFEINE_SMILES` variable
2. **Adjust parameters**: Modify request strings
3. **Add properties**: Request forces, dipoles, etc.
4. **Change methods**: Request different functionals or force fields

### Direct Wrapper Usage

For development and debugging:

```python
from wrappers.quantum_espresso import QuantumEspressoWrapper

# Initialize wrapper
qe = QuantumEspressoWrapper()
print(f"Version: {qe.version}")

# Run calculation
job_params = {
    "calculation": "scf",
    "molecule_smiles": "O",  # Water
    "cutoff_wfc": 50.0
}

result = qe.submit_job(job_params, job_id="test_001")
print(result)
```

### Monitoring Jobs

```python
# Check job status
from wrappers.quantum_espresso import QuantumEspressoWrapper

qe = QuantumEspressoWrapper()
status = qe.get_job_status("job_id_here")
print(status)
```

## Performance Expectations

### By Application

| Application | Setup Time | Calculation Time | Total Time |
|------------|------------|------------------|------------|
| QE         | 2-5s       | 1-10 min         | ~5 min     |
| CP2K       | 2-5s       | 1-15 min         | ~8 min     |
| GPAW       | 2-5s       | 30s-5 min        | ~3 min     |
| LAMMPS     | 1-2s       | 10s-1 min        | ~1 min     |
| GROMACS    | 1-3s       | 30s-2 min        | ~2 min     |

**Total for all five**: ~20-40 minutes

### By Molecule Size

| Atoms | DFT Time | Classical Time |
|-------|----------|----------------|
| 3-5   | 10s-1min | <10s           |
| 10-20 | 1-10min  | 10s-1min       |
| 20-50 | 5-30min  | 30s-5min       |
| 50+   | 30min+   | 1-10min        |

## Next Steps

After successful tests:

1. **Try different molecules**: Test with your catalyst molecules
2. **Adjust parameters**: Fine-tune cutoffs, functionals, etc.
3. **Production runs**: Deploy to HPC with SLURM
4. **Workflow chains**: Combine multiple calculations
5. **Result analysis**: Add visualization and property extraction

## Files Generated

After running tests, check:

```bash
# Job directories (temporary files)
ls /tmp/backplane_jobs/

# Or check configured scratch directory
ls ~/scratch/backplane_jobs/  # If configured

# Each job creates:
# - input files (.in, .inp, .py, .mdp, .gro)
# - output files (.out, .log)
# - result files (.json)
# - metadata (job_metadata.json)
```

## References

- **Main documentation**: `CATALYST_APPLICATIONS.md`
- **Implementation details**: `IMPLEMENTATION_SUMMARY.md`
- **Application manuals**: See application websites

## Support

For issues:
1. Check application logs in output files
2. Review Celery worker output
3. Check FastAPI logs (`uvicorn` terminal)
4. Verify application installation: `APPS/*/`

---

**Last Updated**: November 21, 2025
