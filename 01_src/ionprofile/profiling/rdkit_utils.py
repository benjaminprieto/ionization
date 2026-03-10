"""
rdkit_utils.py - Shared RDKit utilities
=========================================
Low-level functions used by all protonation engines:
    - neutralize_smiles(): strip pre-existing charges
    - get_formal_charge(): count formal charge of a SMILES
    - canonicalize_smiles(): standardize SMILES representation

These are pure helper functions. They don't know about pH or
protonation — that's the engine's job.

Why neutralize?
    SMILES from different sources (SDF, MOL2, Excel) may arrive
    pre-charged from a previous protonation or docking run.
    Dimorphite-DL (and other engines) expect neutral input to
    correctly predict protonation at a target pH. Without
    neutralization, the engine may add charges on top of existing
    ones, producing incorrect results.

Project: ionprofile
"""

import logging
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.debug("RDKit not available")


def neutralize_smiles(smiles: str) -> str:
    """
    Remove pre-existing formal charges from a SMILES string.

    Example:
        "CC[NH+](CC)CC"  ->  "CCN(CC)CC"
        "CC(=O)[O-]"     ->  "CC(=O)O"

    Returns original SMILES if neutralization fails or RDKit
    is not available.
    """
    if not RDKIT_AVAILABLE:
        return smiles
    if pd.isna(smiles) or not smiles:
        return smiles
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return smiles


def get_formal_charge(smiles: str) -> Optional[int]:
    """
    Calculate the formal charge of a SMILES string.

    Example:
        "CCN(CC)CC"       ->  0   (neutral)
        "CC[NH+](CC)CC"   ->  1   (protonated amine)
        "CC(=O)[O-]"      -> -1   (deprotonated acid)

    Returns None if parsing fails or RDKit is not available.
    """
    if not RDKIT_AVAILABLE:
        return None
    if pd.isna(smiles) or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.GetFormalCharge(mol) if mol else None
    except Exception:
        return None


def canonicalize_smiles(smiles: str) -> Optional[str]:
    """
    Canonicalize a SMILES string using RDKit.

    Ensures consistent representation regardless of input format.
    Returns None if the SMILES is invalid.
    """
    if not RDKIT_AVAILABLE:
        return smiles
    if pd.isna(smiles) or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol else None
    except Exception:
        return None


def is_rdkit_available() -> bool:
    """Check if RDKit is installed and importable."""
    return RDKIT_AVAILABLE
