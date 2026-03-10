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


def _smiles_to_mol_3d(smiles: str, mol_id: str = "", keep_hs: bool = False) -> 'Chem.Mol':
    """
    Convert SMILES to RDKit Mol with 3D coordinates.

    Generates a 3D conformer using ETKDG (the best RDKit method).
    Falls back to 2D if 3D generation fails.

    Args:
        smiles: Input SMILES
        mol_id: Molecule identifier (for logging)
        keep_hs: If True, keep explicit hydrogens in output

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
        try:
            result = AllChem.EmbedMolecule(mol, params)
        except RuntimeError:
            result = -1

        if result == 0:
            # Optimize geometry
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                try:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                except Exception:
                    pass
        else:
            # Fallback: 2D coordinates
            AllChem.Compute2DCoords(mol)
            logger.debug(f"  3D failed for '{mol_id}', using 2D")

        # Remove explicit hydrogens unless caller wants them
        if not keep_hs:
            mol = Chem.RemoveHs(mol)

        # Set molecule name
        if mol_id:
            mol.SetProp("_Name", mol_id)

        return mol

    except Exception as e:
        logger.debug(f"  Failed to generate mol for '{mol_id}': {e}")
        return None


def _sanitize_filename(name: str) -> str:
    """Sanitize a string for use as a filename."""
    import re
    safe = re.sub(r'[^\w\-.]', '_', str(name))
    return safe[:200]  # Limit length


def save_sdf(
        df: pd.DataFrame,
        output_dir: Union[str, Path],
        ph_values: List[float],
        output_prefix: str = "ionization_profiling",
) -> Dict[str, str]:
    """
    Save protonated molecules as SDF files (one multi-mol SDF per pH).

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


def save_sdf_individual(
        df: pd.DataFrame,
        output_dir: Union[str, Path],
        ph_values: List[float],
        output_prefix: str = "ionization_profiling",
) -> Dict[str, str]:
    """
    Save individual SDF files: one molecule per file, organized by pH.

    Output structure:
        output_dir/
        ├── pH72/
        │   ├── 5281613.sdf
        │   ├── HTS1710-00236567-01.sdf
        │   └── ...
        ├── pH63/
        │   ├── 5281613.sdf
        │   └── ...
        └── ...

    Each SDF includes:
        - 3D coordinates (ETKDG)
        - Explicit hydrogens (for antechamber compatibility)
        - Properties: mol_id, pH, formal_charge, protonated_smiles, original_smiles

    Args:
        df:            Results DataFrame with SMILES_pH{value} columns.
        output_dir:    Base directory for structures.
        ph_values:     List of pH values.
        output_prefix: Not used (kept for API consistency).

    Returns:
        Dict of pH label -> directory path
    """
    if not RDKIT_AVAILABLE:
        logger.warning("RDKit not available - cannot generate SDF files")
        return {}

    output_dir = Path(output_dir)
    output_files = {}
    total_written = 0

    logger.info(f"Generating individual SDF files (with explicit H)...")

    for ph in ph_values:
        smiles_col = f"SMILES_pH{int(round(ph * 10))}"
        charge_col = f"Q_pH{int(round(ph * 10))}"
        ph_label = f"pH{int(round(ph * 10))}"

        if smiles_col not in df.columns:
            logger.debug(f"  No SMILES column for pH {ph}, skipping")
            continue

        ph_dir = output_dir / ph_label
        ph_dir.mkdir(parents=True, exist_ok=True)
        n_written = 0

        for _, row in df.iterrows():
            prot_smiles = row.get(smiles_col)
            mol_id = str(row.get("mol_id", ""))
            original_smiles = row.get("smiles", "")

            if pd.isna(prot_smiles) or not prot_smiles:
                continue

            # Generate 3D with explicit H
            mol = _smiles_to_mol_3d(str(prot_smiles), mol_id, keep_hs=True)
            if mol is None:
                logger.debug(f"  Failed to generate 3D for {mol_id} at {ph_label}")
                continue

            # Set properties
            mol.SetProp("_Name", mol_id)
            mol.SetProp("mol_id", mol_id)
            mol.SetProp("pH", f"{ph:.1f}")
            mol.SetProp("protonated_smiles", str(prot_smiles))
            mol.SetProp("original_smiles", str(original_smiles))

            charge = row.get(charge_col)
            if pd.notna(charge):
                mol.SetProp("formal_charge", str(int(charge)))

            # Write individual SDF
            filename = _sanitize_filename(mol_id) + ".sdf"
            sdf_path = ph_dir / filename

            writer = Chem.SDWriter(str(sdf_path))
            writer.SetForceV3000(False)
            writer.write(mol)
            writer.close()

            n_written += 1

        if n_written > 0:
            output_files[ph_label] = str(ph_dir)
            total_written += n_written
            logger.info(f"  {ph_label}/: {n_written} individual SDF files")

    logger.info(f"  Total: {total_written} SDF files "
                f"across {len(output_files)} pH folders")

    return output_files