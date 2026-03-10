"""
openbabel_engine.py - OpenBabel protonation engine
=====================================================
pH-dependent protonation using OpenBabel's built-in pH model.

How it works:
    OpenBabel uses ~20-30 SMARTS/pKa rules (phmodel.txt) to add
    hydrogens at a target pH. It modifies the molecule in-place
    and returns the single dominant protonation state.

Best for:    Quick checks, comparison with Dimorphite-DL
Speed:       Milliseconds per molecule
Limitations: Only ~20-30 rules, no enumeration of alternatives,
             known issues with heterocycles and phenols at extreme pH.
             Last updated 2020.

Dependencies: openbabel (install via conda install -c conda-forge openbabel)

Project: ionprofile
"""

import logging
from typing import Optional

import pandas as pd

from ionprofile.profiling.engines.base import BaseEngine

logger = logging.getLogger(__name__)

try:
    from openbabel import pybel, openbabel
    OPENBABEL_AVAILABLE = True
except ImportError:
    try:
        import pybel
        import openbabel
        OPENBABEL_AVAILABLE = True
    except ImportError:
        OPENBABEL_AVAILABLE = False
        logger.debug("OpenBabel not installed")


class OpenBabelEngine(BaseEngine):
    """Protonation engine using OpenBabel's pH model."""

    def __init__(self):
        self._available = OPENBABEL_AVAILABLE
        if self._available:
            # Suppress OpenBabel warnings
            openbabel.obErrorLog.SetOutputLevel(0)

    def calculate_charge_at_ph(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[int]:
        """
        Calculate formal charge at a specific pH using OpenBabel.

        Steps:
            1. Parse SMILES into OpenBabel molecule
            2. Remove existing hydrogens (clean slate)
            3. Add hydrogens at target pH
            4. Read total charge from the molecule
        """
        if not self._available:
            return None
        if pd.isna(smiles) or not smiles:
            return None
        try:
            mol = pybel.readstring("smi", smiles)
            # Remove existing H, then add at target pH
            mol.OBMol.DeleteHydrogens()
            mol.OBMol.AddHydrogens(False, True, ph)
            return mol.OBMol.GetTotalCharge()
        except Exception:
            return None

    def get_protonated_smiles(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[str]:
        """Get protonated SMILES at a specific pH using OpenBabel."""
        if not self._available:
            return None
        if pd.isna(smiles) or not smiles:
            return None
        try:
            mol = pybel.readstring("smi", smiles)
            mol.OBMol.DeleteHydrogens()
            mol.OBMol.AddHydrogens(False, True, ph)
            # Output SMILES without explicit H for cleaner output
            out = pybel.Outputfile("smi", "/dev/null" if __import__('os').name != 'nt' else "NUL")
            # Use write method to get SMILES string
            return mol.write("smi").strip().split("\t")[0]
        except Exception:
            return None

    def is_available(self) -> bool:
        return self._available

    def name(self) -> str:
        return "OpenBabel"