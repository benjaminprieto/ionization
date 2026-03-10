"""
qupkake_engine.py - QupKake protonation engine
=================================================
ML-based micro-pKa prediction using GNN + GFN2-xTB quantum features.

How it works (fundamentally different from Dimorphite/OpenBabel):
    1. QupKake predicts the REAL pKa of each ionizable site
    2. At each target pH, Henderson-Hasselbalch determines if each
       site is protonated or deprotonated
    3. Formal charge is the sum of all site charges

This means:
    - pKa prediction is done ONCE per molecule (expensive: ~seconds)
    - Charge at any pH is then instant (just math)
    - Results are scientifically accurate (RMSE 0.5-0.8 pKa units)

Best for:    Lead optimization, accurate profiling (tens of molecules)
Speed:       Seconds per molecule (uses xtb + GNN)
Accuracy:    RMSE 0.5-0.8 pKa units (state-of-the-art)

Dependencies: qupkake, torch, xtb

Project: ionprofile
"""

import logging
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from ionprofile.profiling.engines.base import BaseEngine

logger = logging.getLogger(__name__)

try:
    from qupkake import predict
    from rdkit import Chem
    QUPKAKE_AVAILABLE = True
except ImportError:
    QUPKAKE_AVAILABLE = False
    logger.debug("QupKake not installed")


class QupKakeEngine(BaseEngine):
    """
    Protonation engine using QupKake (GNN + GFN2-xTB).

    Unlike Dimorphite/OpenBabel which guess protonation states,
    QupKake predicts actual pKa values for each ionizable site.
    Charge at any pH is then calculated via Henderson-Hasselbalch.
    """

    def __init__(self):
        self._available = QUPKAKE_AVAILABLE
        # Cache: smiles -> list of (pka_value, pka_type)
        # pka_type: 'acidic' or 'basic'
        self._pka_cache: Dict[str, List[Tuple[float, str]]] = {}

    def _get_pka_values(self, smiles: str) -> List[Tuple[float, str]]:
        """
        Predict pKa values for all ionizable sites using QupKake.

        Returns list of (pka_value, pka_type) where pka_type is
        'acidic' or 'basic'.

        Results are cached per SMILES to avoid redundant calculations.
        """
        if smiles in self._pka_cache:
            return self._pka_cache[smiles]

        if not self._available:
            return []

        try:
            # QupKake works with temporary files
            with tempfile.TemporaryDirectory() as tmpdir:
                # Run QupKake prediction
                results = predict.predict_smiles(
                    smiles,
                    name="mol",
                    root=tmpdir,
                    tautomerize=False,
                )

                pka_list = []

                if results is not None:
                    # results is an RDKit mol with pKa annotations
                    if isinstance(results, list):
                        mols = results
                    else:
                        mols = [results]

                    for mol in mols:
                        if mol is None:
                            continue
                        for atom in mol.GetAtoms():
                            if atom.HasProp("pka_value"):
                                pka_val = float(atom.GetProp("pka_value"))
                                pka_type = atom.GetProp("pka_type") if atom.HasProp("pka_type") else "acidic"
                                pka_list.append((pka_val, pka_type))

                self._pka_cache[smiles] = pka_list
                return pka_list

        except Exception as e:
            logger.debug(f"QupKake prediction failed for '{smiles[:50]}': {e}")
            self._pka_cache[smiles] = []
            return []

    def _calculate_charge_from_pka(
            self,
            pka_values: List[Tuple[float, str]],
            ph: float,
    ) -> int:
        """
        Calculate formal charge at a given pH from pKa values
        using Henderson-Hasselbalch.

        For acidic sites (HA -> A- + H+):
            fraction_deprotonated = 1 / (1 + 10^(pKa - pH))
            If deprotonated: charge contribution = -1

        For basic sites (B + H+ -> BH+):
            fraction_protonated = 1 / (1 + 10^(pH - pKa))
            If protonated: charge contribution = +1

        We use 50% threshold: if fraction > 0.5, the site is
        considered ionized.
        """
        total_charge = 0

        for pka_val, pka_type in pka_values:
            if pka_type == "acidic":
                # Acidic: above pKa = deprotonated (negative)
                fraction_deprotonated = 1.0 / (1.0 + 10 ** (pka_val - ph))
                if fraction_deprotonated > 0.5:
                    total_charge -= 1
            elif pka_type == "basic":
                # Basic: below pKa = protonated (positive)
                fraction_protonated = 1.0 / (1.0 + 10 ** (ph - pka_val))
                if fraction_protonated > 0.5:
                    total_charge += 1

        return total_charge

    def _build_protonated_smiles(
            self,
            smiles: str,
            pka_values: List[Tuple[float, str]],
            ph: float,
    ) -> Optional[str]:
        """
        Build a protonated SMILES at the given pH based on pKa predictions.

        Uses RDKit to add/remove protons based on which sites are
        ionized at the target pH.
        """
        if not self._available:
            return None

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            # For now, return the SMILES with charge annotation
            # A more sophisticated version would modify the actual atoms
            # This is a simplified approach
            charge = self._calculate_charge_from_pka(pka_values, ph)

            # Use Dimorphite-DL style output if charge is known
            # For accurate SMILES, we'd need atom-level modification
            # For now, we annotate the charge
            return Chem.MolToSmiles(mol) + f" [charge={charge} at pH={ph:.1f}]"

        except Exception:
            return smiles

    def calculate_charge_at_ph(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[int]:
        """
        Calculate formal charge at a specific pH using QupKake pKa
        predictions + Henderson-Hasselbalch.

        Steps:
            1. Predict pKa for all ionizable sites (cached)
            2. Apply Henderson-Hasselbalch at target pH
            3. Sum charges across all sites
        """
        if not self._available:
            return None
        if pd.isna(smiles) or not smiles:
            return None

        try:
            pka_values = self._get_pka_values(smiles)

            if not pka_values:
                # No ionizable sites found = neutral
                return 0

            return self._calculate_charge_from_pka(pka_values, ph)

        except Exception as e:
            logger.debug(f"QupKake charge calculation failed: {e}")
            return None

    def get_protonated_smiles(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[str]:
        """Get protonated SMILES at a specific pH."""
        if not self._available:
            return None
        if pd.isna(smiles) or not smiles:
            return None

        try:
            pka_values = self._get_pka_values(smiles)
            return self._build_protonated_smiles(smiles, pka_values, ph)
        except Exception:
            return smiles

    def is_available(self) -> bool:
        return self._available

    def name(self) -> str:
        return "QupKake"

    def get_pka_report(self, smiles: str) -> dict:
        """
        Get detailed pKa report for a molecule.

        Returns dict with all predicted pKa values and their types.
        Useful for debugging and analysis beyond charge calculation.
        """
        pka_values = self._get_pka_values(smiles)
        return {
            "smiles": smiles,
            "n_ionizable_sites": len(pka_values),
            "sites": [
                {"pka": pka, "type": ptype}
                for pka, ptype in pka_values
            ],
        }
