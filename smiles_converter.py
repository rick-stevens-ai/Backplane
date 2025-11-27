#!/usr/bin/env python3
"""
SMILES Converter Module
Uses oss120 (gpt-oss:120b) to convert chemical formulas/descriptions to SMILES strings
and validates them for correctness before passing to MACE and other tools.
"""
import re
from typing import Optional, Dict, Any, Tuple
import logging

logger = logging.getLogger(__name__)


class SmilesConverter:
    """
    Converts chemical formulas and descriptions to SMILES strings using LLM reasoning.
    """

    def __init__(self, agent=None):
        """
        Initialize the SMILES converter.

        Args:
            agent: ComputationalChemistryAgent instance for LLM access
        """
        self.agent = agent

    def formula_to_smiles(
        self,
        formula: str,
        description: Optional[str] = None,
        validate: bool = True
    ) -> Tuple[Optional[str], Dict[str, Any]]:
        """
        Convert chemical formula/description to SMILES using oss120.

        Args:
            formula: Chemical formula (e.g., "NH3", "Fe(N2)(CO)4", "C6H6")
            description: Optional description for context (e.g., "benzene", "iron carbonyl")
            validate: Whether to validate the SMILES (default: True)

        Returns:
            Tuple of (smiles_string, metadata_dict)
            - smiles_string: Valid SMILES or None if conversion failed
            - metadata: Dict with 'valid', 'error', 'llm_response' keys
        """
        if not self.agent:
            logger.error("No agent provided for SMILES conversion")
            return None, {
                'valid': False,
                'error': 'No agent available',
                'llm_response': None
            }

        # Build the prompt for oss120
        prompt = self._build_conversion_prompt(formula, description)

        try:
            # Use agent to call oss120
            result = self.agent.run_agentic_workflow(prompt)

            # Extract SMILES from response
            smiles = self._extract_smiles_from_response(result)

            if not smiles:
                return None, {
                    'valid': False,
                    'error': 'Could not extract SMILES from LLM response',
                    'llm_response': str(result)
                }

            # Validate if requested
            if validate:
                is_valid, validation_error = self.validate_smiles(smiles)
                if not is_valid:
                    logger.warning(f"Generated SMILES '{smiles}' failed validation: {validation_error}")
                    return None, {
                        'valid': False,
                        'error': f'Invalid SMILES: {validation_error}',
                        'smiles_candidate': smiles,
                        'llm_response': str(result)
                    }

            return smiles, {
                'valid': True,
                'formula': formula,
                'description': description,
                'llm_response': str(result)
            }

        except Exception as e:
            logger.error(f"Error converting formula to SMILES: {str(e)}")
            return None, {
                'valid': False,
                'error': str(e),
                'llm_response': None
            }

    def _build_conversion_prompt(self, formula: str, description: Optional[str]) -> str:
        """Build prompt for oss120 to generate SMILES."""

        base_prompt = f"""Convert the following chemical formula to a valid SMILES string:

Formula: {formula}"""

        if description:
            base_prompt += f"\nDescription: {description}"

        base_prompt += """

Requirements:
1. Provide ONLY the SMILES string as your response
2. Ensure the SMILES is chemically valid
3. For simple molecules (H2O, NH3, CO2), use standard SMILES
4. For organic molecules, use standard organic chemistry SMILES notation
5. For metal complexes:
   - If the complex can be reasonably represented in SMILES, provide it
   - If it's a coordination complex that cannot be properly represented in SMILES (e.g., metal clusters, extended structures), respond with "NOT_REPRESENTABLE_IN_SMILES"

Examples:
- Water (H2O) → O
- Ammonia (NH3) → N
- Benzene (C6H6) → c1ccccc1
- Pyridine (C5H5N) → c1ccncc1
- Metal cluster with 10 Ru atoms → NOT_REPRESENTABLE_IN_SMILES

Provide ONLY the SMILES string or "NOT_REPRESENTABLE_IN_SMILES" as your response."""

        return base_prompt

    def _extract_smiles_from_response(self, result: Any) -> Optional[str]:
        """
        Extract SMILES string from LLM response.

        The response structure from agent.run_agentic_workflow() is:
        {'result': {'message': 'text', 'tools_used': []}}
        """
        try:
            if not result:
                return None

            # Navigate the nested structure
            if isinstance(result, dict):
                if 'result' in result:
                    result = result['result']
                if isinstance(result, dict) and 'message' in result:
                    text = result['message']
                else:
                    text = str(result)
            else:
                text = str(result)

            # Check if LLM indicated it's not representable
            if 'NOT_REPRESENTABLE_IN_SMILES' in text.upper():
                logger.info("LLM indicated molecule not representable in SMILES")
                return None

            # Extract SMILES - look for common patterns
            # SMILES typically contains: letters (C, N, O, etc.), numbers,
            # brackets [], parentheses (), bonds (=, #, -, \, /), rings (digits)

            # Try to find a SMILES-like string
            lines = text.strip().split('\n')
            for line in lines:
                line = line.strip()
                # Skip empty lines and common non-SMILES lines
                if not line or line.startswith('#') or line.lower().startswith('smiles'):
                    continue

                # Remove common prefixes
                line = re.sub(r'^(SMILES|Output|Result|Answer):\s*', '', line, flags=re.IGNORECASE)
                line = line.strip()

                # Check if this looks like a SMILES string
                # SMILES typically has only specific characters
                if re.match(r'^[A-Za-z0-9@+\-\[\]()=#$/\\\.%]+$', line):
                    # Verify it has at least one atom symbol
                    if re.search(r'[CNOPSFBIHcnops]', line):
                        return line

            # If no SMILES-like pattern found, return the first non-empty line
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    return line

            return None

        except Exception as e:
            logger.error(f"Error extracting SMILES: {str(e)}")
            return None

    def validate_smiles(self, smiles: str) -> Tuple[bool, Optional[str]]:
        """
        Validate a SMILES string.

        Args:
            smiles: SMILES string to validate

        Returns:
            Tuple of (is_valid, error_message)
        """
        if not smiles or not isinstance(smiles, str):
            return False, "Empty or invalid SMILES string"

        # Try using RDKit if available
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, "RDKit could not parse SMILES"
            return True, None
        except ImportError:
            # RDKit not available, use basic validation
            logger.debug("RDKit not available, using basic validation")
            pass
        except Exception as e:
            return False, f"RDKit error: {str(e)}"

        # Basic validation without RDKit
        # Check for balanced brackets and parentheses
        if smiles.count('[') != smiles.count(']'):
            return False, "Unbalanced square brackets"
        if smiles.count('(') != smiles.count(')'):
            return False, "Unbalanced parentheses"

        # Check for valid characters (SMILES alphabet)
        valid_chars = set('CNOPSFBrClIHcnops0123456789@+\-=[]()#$/\\.%')
        if not all(c in valid_chars for c in smiles):
            invalid = set(smiles) - valid_chars
            return False, f"Invalid characters: {invalid}"

        # Check it has at least one atom
        if not re.search(r'[A-Za-z]', smiles):
            return False, "No atom symbols found"

        # Passed basic validation
        return True, None

    def xyz_to_smiles(
        self,
        xyz_string: str,
        hint: Optional[str] = None,
        validate: bool = True
    ) -> Tuple[Optional[str], Dict[str, Any]]:
        """
        Attempt to convert XYZ structure to SMILES using oss120.

        This is best-effort - many structures (metal clusters, materials)
        cannot be represented as SMILES.

        Args:
            xyz_string: XYZ format structure
            hint: Optional hint about the molecule (e.g., "organic molecule", "benzene")
            validate: Whether to validate the SMILES

        Returns:
            Tuple of (smiles_string, metadata_dict)
        """
        if not self.agent:
            return None, {
                'valid': False,
                'error': 'No agent available',
                'llm_response': None
            }

        # Extract atoms from XYZ
        atoms = self._parse_xyz_atoms(xyz_string)
        if not atoms:
            return None, {
                'valid': False,
                'error': 'Could not parse XYZ structure',
                'llm_response': None
            }

        # Check if it's likely representable in SMILES
        # (only organic elements, no metals)
        organic_elements = {'H', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'}
        elements_present = set(atoms)

        if not elements_present.issubset(organic_elements):
            metal_elements = elements_present - organic_elements
            logger.info(f"XYZ contains non-organic elements {metal_elements}, likely not SMILES-representable")
            return None, {
                'valid': False,
                'error': f'Contains non-organic elements: {metal_elements}',
                'representable': False,
                'llm_response': None
            }

        # Build prompt
        prompt = f"""Convert the following XYZ atomic structure to a SMILES string:

{xyz_string}

"""
        if hint:
            prompt += f"Hint: This is {hint}\n\n"

        prompt += """Provide ONLY the SMILES string. If this structure cannot be represented as a valid SMILES (e.g., it's a 3D cluster or extended material), respond with "NOT_REPRESENTABLE_IN_SMILES"."""

        try:
            result = self.agent.run_agentic_workflow(prompt)
            smiles = self._extract_smiles_from_response(result)

            if not smiles:
                return None, {
                    'valid': False,
                    'error': 'Could not extract SMILES from LLM response',
                    'llm_response': str(result)
                }

            if validate:
                is_valid, validation_error = self.validate_smiles(smiles)
                if not is_valid:
                    return None, {
                        'valid': False,
                        'error': f'Invalid SMILES: {validation_error}',
                        'smiles_candidate': smiles,
                        'llm_response': str(result)
                    }

            return smiles, {
                'valid': True,
                'xyz': xyz_string[:200] + '...' if len(xyz_string) > 200 else xyz_string,
                'llm_response': str(result)
            }

        except Exception as e:
            logger.error(f"Error converting XYZ to SMILES: {str(e)}")
            return None, {
                'valid': False,
                'error': str(e),
                'llm_response': None
            }

    def _parse_xyz_atoms(self, xyz_string: str) -> list:
        """Extract element symbols from XYZ string."""
        try:
            lines = xyz_string.strip().split('\n')
            if len(lines) < 3:
                return []

            atoms = []
            for line in lines[2:]:  # Skip header and comment
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    atoms.append(parts[0])

            return atoms
        except Exception:
            return []


def create_converter(agent) -> SmilesConverter:
    """
    Factory function to create a SMILES converter.

    Args:
        agent: ComputationalChemistryAgent instance

    Returns:
        SmilesConverter instance
    """
    return SmilesConverter(agent=agent)
