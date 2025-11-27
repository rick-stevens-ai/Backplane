#!/usr/bin/env python3
"""
Simple MACE Client - Communicates with MACE MCP server via subprocess without MCP SDK.

Uses direct JSON-RPC communication over stdin/stdout to avoid requiring MCP SDK
in the Python 3.9 Backplane environment. The MACE server runs in its own
Python 3.10 environment with full MCP support.
"""

import subprocess
import json
import uuid
from typing import List, Dict, Any, Optional


class MACEClient:
    """Client for calling MACE MCP server via subprocess JSON-RPC."""

    def __init__(self, server_path: str = "/Users/stevens/Dropbox/MACE/mace_mcp_server.py"):
        """
        Initialize MACE client.

        Args:
            server_path: Path to mace_mcp_server.py
        """
        self.server_path = server_path
        self.server_python = "python3.10"  # MACE requires Python 3.10+

    def _call_mace_tool(self, tool_name: str, arguments: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Call MACE MCP tool via subprocess JSON-RPC.

        Args:
            tool_name: Name of MACE tool
            arguments: Tool arguments

        Returns:
            Parsed result dictionary or None on error
        """
        try:
            # Start MACE MCP server as subprocess
            process = subprocess.Popen(
                [self.server_python, self.server_path],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # MCP initialization sequence
            init_id = str(uuid.uuid4())
            init_request = {
                "jsonrpc": "2.0",
                "id": init_id,
                "method": "initialize",
                "params": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {},
                    "clientInfo": {
                        "name": "backplane-mace-client",
                        "version": "1.0.0"
                    }
                }
            }

            # Send initialization
            process.stdin.write(json.dumps(init_request) + "\n")
            process.stdin.flush()

            # Read initialization response
            init_response = process.stdout.readline()
            if init_response:
                init_data = json.loads(init_response)
                if "error" in init_data:
                    print(f"MACE server initialization error: {init_data['error']}")
                    process.terminate()
                    return None

            # Call the actual tool
            tool_id = str(uuid.uuid4())
            tool_request = {
                "jsonrpc": "2.0",
                "id": tool_id,
                "method": "tools/call",
                "params": {
                    "name": tool_name,
                    "arguments": arguments
                }
            }

            # Send tool call
            process.stdin.write(json.dumps(tool_request) + "\n")
            process.stdin.flush()

            # Read tool response
            tool_response = process.stdout.readline()
            if not tool_response:
                stderr_output = process.stderr.read()
                print(f"No response from MACE server. Stderr: {stderr_output}")
                process.terminate()
                return None

            tool_data = json.loads(tool_response)

            # Check for errors
            if "error" in tool_data:
                print(f"MACE tool error: {tool_data['error']}")
                process.terminate()
                return None

            # Extract result
            if "result" in tool_data and "content" in tool_data["result"]:
                for content in tool_data["result"]["content"]:
                    if content.get("type") == "text":
                        result = json.loads(content["text"])
                        process.terminate()
                        return result

            process.terminate()
            return None

        except Exception as e:
            print(f"Error calling MACE tool '{tool_name}': {e}")
            import traceback
            traceback.print_exc()
            return None

    def predict_energy(
        self,
        structure: str,
        format: str = "xyz",
        model_type: str = "off",
        model_size: str = "medium"
    ) -> Optional[Dict[str, Any]]:
        """
        Predict molecular energy using MACE.

        Args:
            structure: Molecular structure (XYZ string, SMILES, etc.)
            format: Format of structure ("xyz", "smiles", "pdb", etc.)
            model_type: "off" for organic molecules, "mp" for materials
            model_size: "small", "medium", or "large"

        Returns:
            Dictionary with 'energy', 'unit', 'model' keys or None
        """
        arguments = {
            "atoms_data": structure,
            "format": format,
            "model_type": model_type,
            "model_size": model_size
        }

        return self._call_mace_tool("calculate_energy", arguments)

    def predict_forces(
        self,
        structure: str,
        format: str = "xyz",
        model_type: str = "off",
        model_size: str = "medium"
    ) -> Optional[Dict[str, Any]]:
        """
        Predict atomic forces using MACE.

        Args:
            structure: Molecular structure
            format: Format of structure
            model_type: "off" or "mp"
            model_size: "small", "medium", or "large"

        Returns:
            Dictionary with 'forces', 'max_force', etc. or None
        """
        arguments = {
            "atoms_data": structure,
            "format": format,
            "model_type": model_type,
            "model_size": model_size
        }

        return self._call_mace_tool("calculate_forces", arguments)

    def predict_properties(
        self,
        structure: str,
        format: str = "xyz",
        model_type: str = "off",
        model_size: str = "medium"
    ) -> Optional[Dict[str, Any]]:
        """
        Predict multiple properties (energy + forces) in one call.

        Args:
            structure: Molecular structure
            format: Format of structure
            model_type: "off" or "mp"
            model_size: "small", "medium", or "large"

        Returns:
            Dictionary with 'energy', 'forces', etc. or None
        """
        arguments = {
            "atoms_data": structure,
            "format": format,
            "model_type": model_type,
            "model_size": model_size
        }

        return self._call_mace_tool("predict_properties", arguments)

    def batch_predict_energies(
        self,
        structures: List[str],
        format: str = "xyz",
        model_type: str = "off",
        model_size: str = "medium"
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Predict energies for multiple structures in batch (fast!).

        Args:
            structures: List of molecular structures
            format: Format of structures
            model_type: "off" or "mp"
            model_size: "small", "medium", or "large"

        Returns:
            List of energy dictionaries or None
        """
        # Convert to structure list format
        structure_list = [
            {"atoms_data": struct, "format": format}
            for struct in structures
        ]

        arguments = {
            "structures": structure_list,
            "model_type": model_type,
            "model_size": model_size
        }

        return self._call_mace_tool("batch_calculate_energies", arguments)

    def optimize_geometry(
        self,
        structure: str,
        format: str = "xyz",
        model_type: str = "off",
        model_size: str = "medium",
        fmax: float = 0.05,
        steps: int = 200
    ) -> Optional[Dict[str, Any]]:
        """
        Optimize molecular geometry to minimize energy.

        Args:
            structure: Initial molecular structure
            format: Format of structure
            model_type: "off" or "mp"
            model_size: "small", "medium", or "large"
            fmax: Force convergence criterion (eV/Å)
            steps: Maximum optimization steps

        Returns:
            Dictionary with optimized structure, energy, convergence info
        """
        arguments = {
            "atoms_data": structure,
            "format": format,
            "model_type": model_type,
            "model_size": model_size,
            "fmax": fmax,
            "steps": steps
        }

        return self._call_mace_tool("optimize_geometry", arguments)


def smiles_to_xyz(smiles: str) -> str:
    """
    Convert SMILES to XYZ format using RDKit (reuse from gpaw_wrapper fix).

    Args:
        smiles: SMILES string

    Returns:
        XYZ format string
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=False)

        # Optimize with UFF
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except:
            pass  # Use unoptimized if UFF fails

        # Convert to XYZ
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()

        xyz_lines = [f"{num_atoms}", f"Generated from SMILES: {smiles}"]

        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            xyz_lines.append(f"{symbol:2s}  {pos.x:12.8f}  {pos.y:12.8f}  {pos.z:12.8f}")

        return "\n".join(xyz_lines)

    except Exception as e:
        raise ValueError(f"Failed to convert SMILES to XYZ: {e}")


# Convenience functions for common operations

def quick_energy_check(smiles: str, model_type: str = "off") -> Optional[float]:
    """
    Quick energy check for a SMILES string.

    Args:
        smiles: SMILES string
        model_type: "off" or "mp"

    Returns:
        Energy in eV or None on error
    """
    client = MACEClient()
    xyz = smiles_to_xyz(smiles)
    result = client.predict_energy(xyz, format="xyz", model_type=model_type, model_size="small")

    if result and 'energy' in result:
        return result['energy']
    return None


def rapid_screening(smiles_list: List[str], model_type: str = "off") -> List[Dict[str, Any]]:
    """
    Rapidly screen multiple molecules and rank by energy.

    Args:
        smiles_list: List of SMILES strings
        model_type: "off" or "mp"

    Returns:
        List of dicts with 'smiles', 'energy', 'rank' sorted by energy
    """
    client = MACEClient()

    # Convert all SMILES to XYZ
    xyz_structures = []
    valid_smiles = []

    for smiles in smiles_list:
        try:
            xyz = smiles_to_xyz(smiles)
            xyz_structures.append(xyz)
            valid_smiles.append(smiles)
        except Exception as e:
            print(f"Skipping invalid SMILES '{smiles}': {e}")

    if not xyz_structures:
        return []

    # Batch predict
    results = client.batch_predict_energies(
        xyz_structures,
        format="xyz",
        model_type=model_type,
        model_size="medium"
    )

    if not results:
        return []

    # Combine with SMILES and sort
    screening_results = []
    for smiles, result in zip(valid_smiles, results):
        if result and 'energy' in result:
            screening_results.append({
                'smiles': smiles,
                'energy': result['energy'],
                'unit': result.get('unit', 'eV')
            })

    # Sort by energy (lower is better)
    screening_results.sort(key=lambda x: x['energy'])

    # Add rank
    for i, result in enumerate(screening_results, 1):
        result['rank'] = i

    return screening_results


if __name__ == "__main__":
    # Quick test
    print("Testing MACE client (simple JSON-RPC mode)...")

    client = MACEClient()

    # Test 1: Single energy prediction
    water_xyz = """3
Water molecule
O        0.000000    0.000000    0.119262
H        0.000000    0.763239   -0.477047
H        0.000000   -0.763239   -0.477047
"""

    print("\nTest 1: Single energy prediction")
    result = client.predict_energy(water_xyz)
    if result:
        print(f"  ✓ Energy: {result['energy']} {result['unit']}")
    else:
        print("  ✗ Failed")

    # Test 2: SMILES conversion + energy
    print("\nTest 2: SMILES to energy")
    energy = quick_energy_check("N")  # Ammonia
    if energy:
        print(f"  ✓ Ammonia energy: {energy} eV")
    else:
        print("  ✗ Failed")

    # Test 3: Batch screening
    print("\nTest 3: Rapid screening")
    test_smiles = ["N", "NN", "c1ccncc1"]  # NH3, N2H4, Pyridine
    results = rapid_screening(test_smiles)
    if results:
        print(f"  ✓ Screened {len(results)} molecules:")
        for r in results:
            print(f"     Rank {r['rank']}: {r['smiles']:15s} {r['energy']:.2f} eV")
    else:
        print("  ✗ Failed")

    print("\n✓ MACE client test complete!")
