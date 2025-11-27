# Scientific AI Backplane - Installation Guide

Complete installation guide for all simulation codes, machine learning models, and dependencies.

**Last Updated:** 2025-11-27
**System Requirements:** macOS (Apple Silicon), Linux, or HPC cluster

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Core Dependencies](#core-dependencies)
3. [Quantum Chemistry Codes](#quantum-chemistry-codes)
4. [Classical MD Codes](#classical-md-codes)
5. [Machine Learning Models](#machine-learning-models)
6. [Python Dependencies](#python-dependencies)
7. [Verification](#verification)
8. [Troubleshooting](#troubleshooting)

---

## Quick Start

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/backplane.git
cd backplane

# Install Python dependencies
pip install -r requirements.txt

# Install simulation codes (see sections below)
# Then test installation
python test_all_apps.py
```

---

## Core Dependencies

### 1. Python 3.9+
```bash
# Check Python version
python3 --version

# Should be 3.9 or higher
```

### 2. Redis (for Celery task queue)
```bash
# macOS
brew install redis
brew services start redis

# Linux
sudo apt-get install redis-server
sudo systemctl start redis

# Test Redis
redis-cli ping  # Should return "PONG"
```

### 3. Celery (task queue)
```bash
pip install celery redis
```

---

## Quantum Chemistry Codes

### 1. Quantum ESPRESSO (Plane-wave DFT)

**Use case:** Bulk crystals, surfaces, electrides, periodic systems

#### macOS:
```bash
brew install quantum-espresso
```

#### Linux/HPC:
```bash
# Download from https://www.quantum-espresso.org/
wget https://github.com/QEF/q-e/releases/download/qe-7.2/qe-7.2-ReleasePack.tar.gz
tar -xzf qe-7.2-ReleasePack.tar.gz
cd qe-7.2
./configure
make all
sudo make install
```

#### Verification:
```bash
pw.x -v
# Should display: Program PWSCF v.7.2 starts...
```

---

### 2. CP2K (Mixed Gaussian/Plane-wave DFT)

**Use case:** Molecular dynamics, hybrid DFT, large periodic systems

#### macOS:
```bash
brew install cp2k
```

#### Linux/HPC:
```bash
# Download from https://www.cp2k.org/
git clone --recursive https://github.com/cp2k/cp2k.git
cd cp2k
# Follow compilation instructions in INSTALL.md
```

#### Verification:
```bash
cp2k.ssmp --version
# Should display: CP2K version 2023.2
```

---

### 3. GPAW (Real-space DFT)

**Use case:** Molecules, nanoparticles, surfaces

#### Installation:
```bash
pip install gpaw
# Download pseudopotentials
gpaw install-data
```

#### Verification:
```python
python3 -c "import gpaw; print(gpaw.__version__)"
```

---

### 4. ORCA 6.0 (Molecular DFT/Wavefunction)

**Use case:** Molecular organometallics, coordination complexes, multireference

**IMPORTANT:** ORCA requires manual download (free for academic use)

#### Installation Steps:

1. **Register for Academic Account:**
   - Visit: https://orcaforum.kofo.mpg.de/
   - Register with your academic email
   - Wait for confirmation email

2. **Download ORCA:**
   - Log in to ORCA forum
   - Navigate to Downloads
   - Select version based on your system:
     - **macOS Apple Silicon:** `orca_6_0_0_macosx_arm64_openmpi411.run`
     - **Linux x86_64:** `orca_6_0_0_linux_x86-64_openmpi411.tar.xz`

3. **Install ORCA:**
   ```bash
   # macOS (self-extracting installer)
   chmod +x orca_6_0_0_macosx_arm64_openmpi411.run
   ./orca_6_0_0_macosx_arm64_openmpi411.run

   # Or extract to specific location
   mkdir -p ~/APPS/orca
   cp -r orca_6_0_0/ ~/APPS/orca/

   # Linux (tar archive)
   tar -xf orca_6_0_0_linux_x86-64_openmpi411.tar.xz
   mv orca_6_0_0 ~/APPS/orca

   # Make binaries executable
   chmod +x ~/APPS/orca/orca
   chmod +x ~/APPS/orca/orca_*
   ```

4. **Add to PATH (optional):**
   ```bash
   # Add to ~/.bashrc or ~/.zshrc
   export PATH=$PATH:~/APPS/orca
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/APPS/orca
   ```

#### Verification:
```bash
~/APPS/orca/orca
# Should display: This program requires the name of a parameterfile as argument

# Test calculation
cat > test.inp << 'EOF'
! B3LYP def2-SVP
* xyz 0 1
O   0.0   0.0   0.0
H   0.76  0.59  0.0
H  -0.76  0.59  0.0
*
EOF

~/APPS/orca/orca test.inp > test.out
grep "FINAL SINGLE POINT ENERGY" test.out
# Should display: FINAL SINGLE POINT ENERGY       -76.321089...
```

---

### 5. NWChem (Parallel Quantum Chemistry)

**Use case:** Large clusters (100-500 atoms), parallel DFT/MCSCF

#### macOS:
```bash
brew install nwchem
```

#### Linux/HPC:
```bash
# Download from https://nwchemgit.github.io/
wget https://github.com/nwchemgit/nwchem/releases/download/v7.3.1-release/nwchem-7.3.1-release.tar.gz
tar -xzf nwchem-7.3.1-release.tar.gz
cd nwchem-7.3.1
# Follow INSTALL instructions
```

#### Verification:
```bash
which nwchem
# Should display: /opt/homebrew/bin/nwchem (or similar)
```

---

### 6. PySCF (Python Quantum Chemistry)

**Use case:** Multireference calculations (CASSCF, NEVPT2), method development

#### Installation:
```bash
pip install pyscf
```

#### Verification:
```python
python3 << 'EOF'
from pyscf import gto, scf, mcscf
from pyscf import mp, cc, fci
from pyscf.mrpt import nevpt2
print("✓ PySCF installed successfully")
print(f"  Version: {pyscf.__version__ if hasattr(pyscf, '__version__') else 'OK'}")
EOF
```

---

### 7. ASE (Atomic Simulation Environment)

**Use case:** Workflow automation, structure manipulation, NEB calculations

#### Installation:
```bash
pip install ase
```

#### Verification:
```python
python3 << 'EOF'
from ase import Atoms
from ase.optimize import BFGS
from ase.neb import NEB
print("✓ ASE installed successfully")
EOF
```

---

## Classical MD Codes

### 1. LAMMPS (Large-scale Molecular Dynamics)

**Use case:** Classical MD, force fields, large systems

#### macOS:
```bash
brew install lammps
```

#### Linux/HPC:
```bash
# Download from https://www.lammps.org/
git clone -b stable https://github.com/lammps/lammps.git lammps
cd lammps
mkdir build && cd build
cmake ../cmake
make -j4
sudo make install
```

#### Verification:
```bash
lmp -help
# Should display LAMMPS version and compilation info
```

---

### 2. GROMACS (MD for Biomolecules)

**Use case:** Biomolecular MD, proteins, explicit solvent

#### macOS:
```bash
brew install gromacs
```

#### Linux/HPC:
```bash
wget http://ftp.gromacs.org/pub/gromacs/gromacs-2023.2.tar.gz
tar -xzf gromacs-2023.2.tar.gz
cd gromacs-2023.2
mkdir build && cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON
make -j4
sudo make install
```

#### Verification:
```bash
gmx --version
# Should display GROMACS version 2023.2
```

---

## Machine Learning Models

### 1. MACE (Machine Learning Force Field)

**Use case:** Fast ML predictions, pre-screening, molecular dynamics

#### Installation:

1. **Install MACE Python Package:**
   ```bash
   pip install mace-torch
   ```

2. **Download MACE-MP-0 Model (Materials Project):**
   ```bash
   # Create models directory
   mkdir -p ~/MODELS/mace
   cd ~/MODELS/mace

   # Download MACE-MP-0 (medium, balanced speed/accuracy)
   wget https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0/MACE-MP-0-medium.model
   ```

3. **Alternative: Download in Python:**
   ```python
   from mace.calculators import mace_mp

   # This will auto-download the model
   calc = mace_mp(model="medium", device="cpu")
   ```

4. **Set Environment Variable (optional):**
   ```bash
   # Add to ~/.bashrc or ~/.zshrc
   export MACE_MODEL_PATH=~/MODELS/mace/MACE-MP-0-medium.model
   ```

#### Available MACE Models:
- **small:** 0.5M parameters, fastest, lower accuracy
- **medium:** 5.5M parameters, balanced (RECOMMENDED)
- **large:** 88M parameters, highest accuracy, slower

#### Verification:
```python
python3 << 'EOF'
from mace.calculators import mace_mp
from ase import Atoms

# Create simple molecule
atoms = Atoms('H2O', positions=[[0,0,0], [0.76,0.59,0], [-0.76,0.59,0]])

# Setup MACE calculator
calc = mace_mp(model="medium", device="cpu")
atoms.calc = calc

# Predict energy
energy = atoms.get_potential_energy()
print(f"✓ MACE working! H2O energy: {energy:.3f} eV")
EOF
```

---

## Python Dependencies

### Core Requirements

Create `requirements.txt`:
```txt
# Core dependencies
celery>=5.3.0
redis>=4.6.0
requests>=2.31.0
pyyaml>=6.0
fastapi>=0.104.0
uvicorn>=0.24.0

# Scientific computing
numpy>=1.24.0
scipy>=1.11.0
pandas>=2.0.0
matplotlib>=3.7.0

# Quantum chemistry
ase>=3.22.0
pyscf>=2.3.0

# Machine learning
torch>=2.0.0
mace-torch>=0.3.0

# File parsing
pymatgen>=2023.0.0  # For CIF parsing (optional)
```

### Installation:
```bash
pip install -r requirements.txt
```

---

## Verification

### Test All Installations

Run the comprehensive test suite:
```bash
python test_all_apps.py
```

Expected output:
```
========================================
Testing Simulation Code Installations
========================================

✓ Quantum ESPRESSO: PASSED
✓ CP2K: PASSED
✓ GPAW: PASSED
✓ ORCA: PASSED
✓ NWChem: PASSED
✓ PySCF: PASSED
✓ ASE: PASSED
✓ LAMMPS: PASSED
✓ GROMACS: PASSED
✓ MACE: PASSED

Summary: 10/10 codes installed and working
```

### Test Individual Components

```bash
# Test Quantum ESPRESSO
python test_qe_only.py

# Test MACE
python test_mace_integration.py

# Test agent system
python test_agent_minimal.py
```

---

## Directory Structure

```
backplane/
├── main.py                 # FastAPI server
├── tasks.py                # Celery tasks
├── agent_apps.py           # Agentic workflows
├── requirements.txt        # Python dependencies
│
├── wrappers/               # Simulation code wrappers
│   ├── quantum_espresso.py
│   ├── cp2k_wrapper.py
│   ├── gpaw_wrapper.py
│   ├── orca_wrapper.py     # NEW
│   ├── nwchem_wrapper.py   # NEW
│   ├── pyscf_wrapper.py    # NEW
│   ├── lammps_wrapper.py
│   └── gromacs_wrapper.py
│
├── cif_parser.py           # Crystallographic structure parser
├── smiles_converter.py     # SMILES to 3D structure
├── mace_client.py          # MACE ML predictions
│
└── spark_servers.yaml      # HPC cluster configuration
```

---

## Troubleshooting

### ORCA Not Found
```bash
# Check if ORCA is in the correct location
ls -lh ~/APPS/orca/orca

# If not found, re-download and extract
# Make sure all binaries are executable
chmod +x ~/APPS/orca/*
```

### NWChem MPI Errors
```bash
# NWChem requires MPI for parallel execution
# For single-core testing, this is usually OK

# Check NWChem location
which nwchem
```

### MACE Model Not Found
```bash
# Download model manually
mkdir -p ~/MODELS/mace
cd ~/MODELS/mace
wget https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0/MACE-MP-0-medium.model

# Set environment variable
export MACE_MODEL_PATH=~/MODELS/mace/MACE-MP-0-medium.model
```

### Redis Connection Failed
```bash
# Start Redis server
redis-server

# Or on macOS with Homebrew
brew services start redis

# Test connection
redis-cli ping
```

### Celery Worker Not Starting
```bash
# Make sure Redis is running first
redis-cli ping

# Start Celery worker
celery -A tasks worker --loglevel=info

# Check for Python import errors
python3 -c "import tasks"
```

### Python Import Errors
```bash
# Reinstall dependencies
pip install --upgrade -r requirements.txt

# Check specific packages
python3 -c "import ase; print(ase.__version__)"
python3 -c "import pyscf; print('PySCF OK')"
python3 -c "from mace.calculators import mace_mp; print('MACE OK')"
```

---

## Platform-Specific Notes

### macOS Apple Silicon (M1/M2/M3)

- Quantum ESPRESSO: Use Homebrew version (ARM64 native)
- CP2K: Use Homebrew version
- ORCA: Download ARM64-specific version
- LAMMPS: Homebrew installs ARM64 native version

### Linux x86_64

- All codes have native x86_64 support
- Compilation from source recommended for HPC systems
- Enable AVX2/AVX512 optimizations during compilation

### HPC Clusters

- Use module system if available:
  ```bash
  module load quantum-espresso
  module load cp2k
  module load orca
  ```
- Configure `spark_servers.yaml` with cluster details
- Set up SSH key authentication for job submission

---

## Performance Optimization

### Quantum ESPRESSO
```bash
# Use OpenMP threading
export OMP_NUM_THREADS=4

# Use k-point parallelization for large systems
mpirun -np 16 pw.x -nk 4 -i input.in
```

### ORCA
```bash
# Enable parallel execution in input file
%pal nprocs 8 end

# Set memory limit
%maxcore 2000  # MB per core
```

### MACE
```python
# Use GPU if available
calc = mace_mp(model="medium", device="cuda")

# Or use CPU
calc = mace_mp(model="medium", device="cpu")
```

---

## Additional Resources

- **Quantum ESPRESSO:** https://www.quantum-espresso.org/Doc/user_guide/
- **CP2K:** https://manual.cp2k.org/
- **ORCA:** https://www.orcasoftware.de/tutorials_orca/
- **NWChem:** https://nwchemgit.github.io/Home.html
- **PySCF:** https://pyscf.org/user.html
- **ASE:** https://wiki.fysik.dtu.dk/ase/
- **MACE:** https://github.com/ACEsuit/mace
- **LAMMPS:** https://docs.lammps.org/
- **GROMACS:** https://manual.gromacs.org/

---

## Support

For issues specific to the Backplane system:
- Open an issue on GitHub
- Check existing documentation in `/docs`

For simulation code issues:
- Refer to official documentation (links above)
- Check code-specific forums and mailing lists

---

**Installation Guide Version:** 1.0
**Last Updated:** 2025-11-27
**Maintained By:** Backplane Development Team
