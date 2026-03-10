"""
qupkake_engine.py - QupKake protonation engine
=================================================
ML-based micro-pKa prediction using GNN + GFN2-xTB quantum features.

How it works:
    1. Calls QupKake CLI to predict pKa for each ionizable site
    2. Parses the output SDF: pKa value, type, AND atom index
    3. Henderson-Hasselbalch determines which sites are ionized at each pH
    4. RDKit modifies the actual molecule: adds/removes H, sets charges
    5. Returns real protonated SMILES (e.g. CC[NH+](CC)CC, not annotations)

QupKake SDF output format (molecule-level properties):
    >  <idx>       atom index of ionizable site
    >  <pka_type>  "basic" or "acidic"
    >  <pka>       predicted pKa value

Best for:    Lead optimization (tens of molecules)
Speed:       Seconds to minutes per molecule (xtb + GNN)
Accuracy:    RMSE 0.5-0.8 pKa units (state-of-the-art)

Dependencies: qupkake (CLI), torch, xtb, rdkit

Project: ionprofile
"""

import logging
import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

import pandas as pd

from ionprofile.profiling.engines.base import BaseEngine

logger = logging.getLogger(__name__)

# Check if qupkake CLI is available
try:
    _result = subprocess.run(
        ["qupkake", "--help"],
        capture_output=True, timeout=10,
    )
    QUPKAKE_AVAILABLE = True
except (FileNotFoundError, subprocess.TimeoutExpired):
    QUPKAKE_AVAILABLE = False
    logger.debug("QupKake CLI not found")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolStandardize
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


class QupKakeEngine(BaseEngine):
    """
    Protonation engine using QupKake (GNN + GFN2-xTB).

    Calls QupKake CLI, parses the output SDF for pKa predictions
    (including atom indices), and generates real protonated SMILES
    by modifying the molecule at each pH via RDKit.
    """

    def __init__(self):
        self._available = QUPKAKE_AVAILABLE and RDKIT_AVAILABLE
        # Cache: smiles -> list of (pka_value, pka_type, atom_index)
        self._pka_cache: Dict[str, List[Tuple[float, str, int]]] = {}

    def _predict_pka(self, smiles: str) -> List[Tuple[float, str, int]]:
        """
        Run QupKake CLI and extract pKa predictions with atom indices.

        Returns list of (pka_value, pka_type, atom_index) tuples.
        Results are cached per SMILES.
        """
        if smiles in self._pka_cache:
            return self._pka_cache[smiles]

        if not self._available:
            return []

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_name = "pka_output.sdf"

                cmd = [
                    "qupkake", "smiles", smiles,
                    "-n", "mol",
                    "-o", output_name,
                    "-r", tmpdir,
                ]

                logger.info(f"  QupKake: predicting pKa for "
                            f"'{smiles[:60]}...'")

                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=1800,
                )

                if result.returncode != 0:
                    logger.warning(f"  QupKake CLI failed: "
                                   f"{result.stderr[:200]}")
                    self._pka_cache[smiles] = []
                    return []

                output_path = os.path.join(tmpdir, "output", output_name)
                if not os.path.exists(output_path):
                    logger.warning(f"  QupKake output not found at "
                                   f"{output_path}")
                    self._pka_cache[smiles] = []
                    return []

                pka_list = self._parse_sdf_pka(output_path)
                self._pka_cache[smiles] = pka_list

                if pka_list:
                    logger.info(f"  QupKake: found {len(pka_list)} "
                                f"ionizable sites")
                    for pka_val, pka_type, atom_idx in pka_list:
                        logger.info(f"    atom {atom_idx}: "
                                    f"pKa = {pka_val:.2f} ({pka_type})")
                else:
                    logger.info(f"  QupKake: no ionizable sites found")

                return pka_list

        except subprocess.TimeoutExpired:
            logger.warning(f"  QupKake timed out for '{smiles[:50]}...'")
            self._pka_cache[smiles] = []
            return []
        except Exception as e:
            logger.warning(f"  QupKake failed: {e}")
            self._pka_cache[smiles] = []
            return []

    def _parse_sdf_pka(
            self, sdf_path: str
    ) -> List[Tuple[float, str, int]]:
        """
        Parse pKa values from QupKake output SDF.

        Extracts three properties per entry:
            idx      -> atom index
            pka_type -> "basic" or "acidic"
            pka      -> predicted pKa value
        """
        pka_list = []

        try:
            supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)

            for mol in supplier:
                if mol is None:
                    continue

                props = mol.GetPropsAsDict()

                pka_val = props.get("pka")
                pka_type = props.get("pka_type")
                atom_idx = props.get("idx")

                if pka_val is not None and pka_type is not None:
                    try:
                        idx = int(atom_idx) if atom_idx is not None else -1
                        pka_list.append(
                            (float(pka_val), str(pka_type), idx)
                        )
                    except (ValueError, TypeError):
                        pass

        except Exception as e:
            logger.warning(f"  Failed to parse QupKake SDF: {e}")

        return pka_list

    def _charge_from_pka(
            self,
            pka_values: List[Tuple[float, str, int]],
            ph: float,
    ) -> int:
        """
        Calculate formal charge using Henderson-Hasselbalch.

        Acidic (HA -> A- + H+): pH > pKa → deprotonated → charge -1
        Basic (B + H+ -> BH+): pH < pKa → protonated → charge +1
        """
        total_charge = 0

        for pka_val, pka_type, _ in pka_values:
            if pka_type == "acidic":
                fraction_deprotonated = 1.0 / (1.0 + 10 ** (pka_val - ph))
                if fraction_deprotonated > 0.5:
                    total_charge -= 1
            elif pka_type == "basic":
                fraction_protonated = 1.0 / (1.0 + 10 ** (ph - pka_val))
                if fraction_protonated > 0.5:
                    total_charge += 1

        return total_charge

    def _get_ionized_sites(
            self,
            pka_values: List[Tuple[float, str, int]],
            ph: float,
    ) -> List[Tuple[int, str]]:
        """
        Determine which sites are ionized at a given pH.

        Returns list of (atom_index, action) where action is:
            "deprotonate" -> remove H, set charge -1
            "protonate"   -> add H, set charge +1
        """
        ionized = []

        for pka_val, pka_type, atom_idx in pka_values:
            if atom_idx < 0:
                continue

            if pka_type == "acidic":
                fraction_deprotonated = 1.0 / (1.0 + 10 ** (pka_val - ph))
                if fraction_deprotonated > 0.5:
                    ionized.append((atom_idx, "deprotonate"))
            elif pka_type == "basic":
                fraction_protonated = 1.0 / (1.0 + 10 ** (ph - pka_val))
                if fraction_protonated > 0.5:
                    ionized.append((atom_idx, "protonate"))

        return ionized

    def _build_protonated_smiles(
            self,
            smiles: str,
            pka_values: List[Tuple[float, str, int]],
            ph: float,
    ) -> Optional[str]:
        """
        Build a real protonated SMILES at the given pH.

        Uses RDKit to modify the molecule:
        - Acidic sites above pKa: remove H, set formal charge -1
        - Basic sites below pKa: add H, set formal charge +1
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return smiles

            # Get which sites are ionized at this pH
            ionized_sites = self._get_ionized_sites(pka_values, ph)

            if not ionized_sites:
                return smiles  # No changes needed

            # Work with explicit H so we can add/remove them
            mol = Chem.RWMol(Chem.AddHs(mol))

            for atom_idx, action in ionized_sites:
                if atom_idx >= mol.GetNumAtoms():
                    continue

                atom = mol.GetAtomWithIdx(atom_idx)

                if action == "deprotonate":
                    # Find an H bonded to this atom and remove it
                    h_to_remove = None
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 1:
                            h_to_remove = neighbor.GetIdx()
                            break

                    if h_to_remove is not None:
                        mol.RemoveAtom(h_to_remove)
                        # After removal, atom indices shift
                        # Recalculate the target atom index
                        if atom_idx > h_to_remove:
                            atom_idx -= 1
                        atom = mol.GetAtomWithIdx(atom_idx)

                    atom.SetFormalCharge(-1)

                elif action == "protonate":
                    # Add H to this atom
                    h_idx = mol.AddAtom(Chem.Atom(1))
                    mol.AddBond(atom_idx, h_idx, Chem.BondType.SINGLE)
                    atom.SetFormalCharge(1)

            # Convert back to SMILES
            try:
                Chem.SanitizeMol(mol)
                mol_no_h = Chem.RemoveHs(mol)
                return Chem.MolToSmiles(mol_no_h)
            except Exception:
                # If sanitization fails, try without removing H
                try:
                    return Chem.MolToSmiles(mol)
                except Exception:
                    # Last resort: return original with charge annotation
                    charge = self._charge_from_pka(pka_values, ph)
                    return f"{smiles} [charge={charge}]"

        except Exception as e:
            logger.debug(f"Failed to build protonated SMILES: {e}")
            charge = self._charge_from_pka(pka_values, ph)
            return f"{smiles} [charge={charge}]"

    def calculate_charge_at_ph(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[int]:
        """
        Calculate formal charge at a specific pH.

        Steps:
            1. Predict pKa via QupKake CLI (cached)
            2. Apply Henderson-Hasselbalch at target pH
            3. Sum charges across all ionizable sites
        """
        if not self._available:
            return None
        if pd.isna(smiles) or not smiles:
            return None

        try:
            pka_values = self._predict_pka(smiles)

            if not pka_values:
                return 0

            return self._charge_from_pka(pka_values, ph)

        except Exception as e:
            logger.debug(f"QupKake charge failed: {e}")
            return None

    def get_protonated_smiles(
            self,
            smiles: str,
            ph: float,
            precision: float = 0.5,
    ) -> Optional[str]:
        """
        Get the REAL protonated SMILES at a specific pH.

        Modifies the molecule structure using RDKit based on
        which sites are ionized at the target pH.
        """
        if not self._available:
            return None
        if pd.isna(smiles) or not smiles:
            return None

        try:
            pka_values = self._predict_pka(smiles)

            if not pka_values:
                return smiles

            return self._build_protonated_smiles(smiles, pka_values, ph)

        except Exception:
            return smiles

    def is_available(self) -> bool:
        return self._available

    def name(self) -> str:
        return "QupKake"

    def get_pka_report(self, smiles: str) -> dict:
        """Get detailed pKa report with atom indices."""
        pka_values = self._predict_pka(smiles)
        return {
            "smiles": smiles,
            "n_ionizable_sites": len(pka_values),
            "sites": [
                {"atom_idx": idx, "pka": round(pka, 2), "type": ptype}
                for pka, ptype, idx in pka_values
            ],
        }