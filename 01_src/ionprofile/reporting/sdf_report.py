"""
sdf_report.py - SDF structure output
=======================================
Generates SDF files with protonated molecules at each pH value.
Each pH produces a separate SDF file with 3D coordinates.

Output structure:
    structures/
    ├── ionization_profiling_pH74.sdf   (all molecules at pH 7.4)
    ├── ionization_profiling_pH72.sdf   (all molecules at pH 7.2)
    ├── ionization_profiling_pH70.sdf   (all molecules at pH 7.0)
    └── ...

Each molecule in the SDF includes properties:
    - mol_id
    - pH
    - formal_charge
    - original_smiles

Project: ionprofile
"""

import logging
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def _smiles_to_mol_3d(smiles: str, mol_id: str = "") -> 'Chem.Mol':
    """
    Convert SMILES to RDKit Mol with 3D coordinates.

    Generates a 3D conformer using ETKDG (the best RDKit method).
    Falls back to 2D if 3D generation fails.

    Returns None if SMILES is invalid.
    """
    if not RDKIT_AVAILABLE:
        return None
    if pd.isna(smiles) or not smiles:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Add hydrogens for 3D generation
        mol = Chem.AddHs(mol)

        # Try 3D coordinate generation
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)

        if result == 0:
            # Optimize geometry
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        else:
            # Fallback: 2D coordinates
            AllChem.Compute2DCoords(mol)
            logger.debug(f"  3D failed for '{mol_id}', using 2D")

        # Remove explicit hydrogens for cleaner SDF
        mol = Chem.RemoveHs(mol)

        # Set molecule name
        if mol_id:
            mol.SetProp("_Name", mol_id)

        return mol

    except Exception as e:
        logger.debug(f"  Failed to generate mol for '{mol_id}': {e}")
        return None


def save_sdf(
        df: pd.DataFrame,
        output_dir: Union[str, Path],
        ph_values: List[float],
        output_prefix: str = "ionization_profiling",
) -> Dict[str, str]:
    """
    Save protonated molecules as SDF files (one per pH).

    Args:
        df:            Results DataFrame (must have SMILES_pH{value} columns).
        output_dir:    Directory for SDF files.
        ph_values:     List of pH values.
        output_prefix: Base name for SDF files.

    Returns:
        Dict of pH label -> filepath (e.g. {"pH_7.4": "path/to/file.sdf"})
    """
    if not RDKIT_AVAILABLE:
        logger.warning("RDKit not available - cannot generate SDF files")
        return {}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = {}
    total_written = 0

    for ph in ph_values:
        smiles_col = f"SMILES_pH{int(round(ph * 10))}"
        charge_col = f"Q_pH{int(round(ph * 10))}"

        if smiles_col not in df.columns:
            logger.debug(f"  No SMILES column for pH {ph}, skipping SDF")
            continue

        # Build SDF filename
        ph_label = f"pH{int(round(ph * 10))}"
        sdf_path = output_dir / f"{output_prefix}_{ph_label}.sdf"

        writer = Chem.SDWriter(str(sdf_path))
        n_written = 0

        for _, row in df.iterrows():
            prot_smiles = row.get(smiles_col)
            mol_id = row.get("mol_id", "")
            original_smiles = row.get("smiles", "")

            if pd.isna(prot_smiles) or not prot_smiles:
                continue

            mol = _smiles_to_mol_3d(str(prot_smiles), mol_id)
            if mol is None:
                continue

            # Add properties to SDF
            mol.SetProp("mol_id", str(mol_id))
            mol.SetProp("pH", f"{ph:.1f}")
            mol.SetProp("protonated_smiles", str(prot_smiles))
            mol.SetProp("original_smiles", str(original_smiles))

            charge = row.get(charge_col)
            if pd.notna(charge):
                mol.SetProp("formal_charge", str(int(charge)))

            writer.write(mol)
            n_written += 1

        writer.close()

        if n_written > 0:
            output_files[f"pH_{ph:.1f}"] = str(sdf_path)
            total_written += n_written
            logger.info(f"  Saved SDF: {sdf_path.name} "
                        f"({n_written} molecules)")
        else:
            # Remove empty file
            sdf_path.unlink(missing_ok=True)

    logger.info(f"  Total SDF molecules written: {total_written} "
                f"across {len(output_files)} pH values")

    return output_files
