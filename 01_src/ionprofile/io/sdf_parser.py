"""
sdf_parser.py - SDF/MOL file reader
=====================================
Reads SDF (Structure Data Format) and MOL files using RDKit.
Extracts canonical SMILES from each molecule entry.

SDF files can contain multiple molecules; MOL files contain a single one.
SMILES generated from the molecular graph are more reliable than
SMILES copied from spreadsheets or text files.

Project: ionprofile
"""

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def parse_sdf_file(filepath: str) -> pd.DataFrame:
    """
    Parse SDF or MOL file and extract SMILES.

    Uses RDKit's SDMolSupplier. For each molecule:
    - Extracts canonical SMILES via Chem.MolToSmiles
    - Uses the _Name property as mol_id (or generates one)

    Args:
        filepath: Path to .sdf or .mol file.

    Returns:
        DataFrame with columns: mol_id, smiles

    Raises:
        ImportError: If RDKit is not installed.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required to read SDF/MOL files. "
            "Install with: conda install -c conda-forge rdkit"
        )

    path = Path(filepath)
    logger.info(f"  Parsing SDF: {path.name}")

    # Suppress RDKit warnings during parsing
    from rdkit import RDLogger
    rd_logger = RDLogger.logger()
    original_level = RDLogger.WARNING
    rd_logger.setLevel(RDLogger.ERROR)

    try:
        supplier = Chem.SDMolSupplier(str(filepath), removeHs=True)

        records = []
        n_failed = 0

        for i, mol in enumerate(supplier):
            if mol is None:
                n_failed += 1
                logger.debug(f"  Failed to parse molecule at index {i}")
                continue

            # Extract canonical SMILES from molecular graph
            smiles = Chem.MolToSmiles(mol)

            # Extract name
            mol_id = mol.GetProp("_Name") if mol.HasProp("_Name") else None
            if not mol_id or not mol_id.strip():
                mol_id = f"sdf_{i+1:04d}"
            else:
                mol_id = mol_id.strip()

            records.append({"mol_id": mol_id, "smiles": smiles})

        if n_failed > 0:
            logger.warning(f"  {n_failed} molecules failed to parse in "
                           f"{path.name}")

        logger.info(f"  Extracted {len(records)} molecules from "
                    f"{path.name}")

    finally:
        rd_logger.setLevel(original_level)

    return pd.DataFrame(records)