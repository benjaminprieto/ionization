"""
dimorphite_engine.py - Dimorphite-DL protonation engine
=========================================================
Fast empirical protonation using 38 SMARTS substructure patterns.

Best for:    High-throughput screening (thousands of molecules)
Speed:       Milliseconds per molecule
Limitations: Does not predict real pKa values, only enumerates
             protonation states based on empirical rules.

Workflow:
    1. Neutralize input SMILES (via rdkit_utils)
    2. Protonate at target pH (via dimorphite_dl)
    3. Calculate formal charge of result (via rdkit_utils)

Dependencies: dimorphite_dl, rdkit

Project: ionprofile
"""

import logging
from typing import Optional

import pandas as pd

from ionprofile.profiling.engines.base import BaseEngine
from ionprofile.profiling.rdkit_utils import (
    neutralize_smiles,
    get_formal_charge,
    is_rdkit_available,
)

logger = logging.getLogger(__name__)


class DimorphiteEngine(BaseEngine):
    """Protonation engine using Dimorphite-DL."""

    def __init__(self):
        self._dimorphite_available = False
        try:
            from dimorphite_dl import protonate_smiles
            self._protonate_fn = protonate_smiles
            self._dimorphite_available = True
        except ImportError:
            logger.debug("dimorphite_dl not installed")

    def calculate_charge_at_ph(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[int]:
        """
        Calculate formal charge at a specific pH.

        Steps:
            1. Neutralize (strip pre-existing charges)
            2. Protonate at target pH via Dimorphite-DL
            3. Count formal charge via RDKit
        """
        if not self.is_available():
            return None
        if pd.isna(smiles) or not smiles:
            return None
        try:
            neutral = neutralize_smiles(smiles)
            forms = self._protonate_fn(
                neutral, ph_min=ph, ph_max=ph, precision=precision
            )
            return get_formal_charge(forms[0]) if forms else None
        except Exception:
            return None

    def get_protonated_smiles(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[str]:
        """Get the protonated SMILES form at a specific pH."""
        if not self._dimorphite_available:
            return None
        if pd.isna(smiles) or not smiles:
            return None
        try:
            neutral = neutralize_smiles(smiles)
            forms = self._protonate_fn(
                neutral, ph_min=ph, ph_max=ph, precision=precision
            )
            return forms[0] if forms else None
        except Exception:
            return None

    def is_available(self) -> bool:
        """Dimorphite needs both dimorphite_dl AND RDKit."""
        return self._dimorphite_available and is_rdkit_available()

    def name(self) -> str:
        return "Dimorphite-DL"
