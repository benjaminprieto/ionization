"""
mol2_parser.py - MOL2 file reader
===================================
Reads Tripos MOL2 files. These are common in docking workflows
(e.g. GOLD, DOCK) and can contain multiple molecules.

Uses RDKit's MolFromMol2Block for parsing.

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


def _split_mol2_blocks(filepath: str) -> list:
    """
    Split a multi-molecule MOL2 file into individual molecule blocks.

    MOL2 files use @<TRIPOS>MOLECULE as the delimiter between entries.
    """
    blocks = []
    current_block = []

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.strip().startswith("@<TRIPOS>MOLECULE"):
                if current_block:
                    blocks.append("".join(current_block))
                current_block = [line]
            else:
                current_block.append(line)

    if current_block:
        blocks.append("".join(current_block))

    return blocks


def _extract_name_from_block(block: str) -> str:
    """Extract molecule name from MOL2 block (line after @<TRIPOS>MOLECULE)."""
    lines = block.strip().split("\n")
    for i, line in enumerate(lines):
        if "@<TRIPOS>MOLECULE" in line and i + 1 < len(lines):
            name = lines[i + 1].strip()
            if name:
                return name
    return ""


def parse_mol2_file(filepath: str) -> pd.DataFrame:
    """
    Parse MOL2 file and extract SMILES via RDKit.

    Args:
        filepath: Path to .mol2 file (single or multi-molecule).

    Returns:
        DataFrame with columns: mol_id, smiles

    Raises:
        ImportError: If RDKit is not installed.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required to read MOL2 files. "
            "Install with: conda install -c conda-forge rdkit"
        )

    path = Path(filepath)
    logger.info(f"  Parsing MOL2: {path.name}")

    # Suppress RDKit warnings
    from rdkit import RDLogger
    rd_logger = RDLogger.logger()
    original_level = rd_logger.level()
    rd_logger.setLevel(RDLogger.ERROR)

    try:
        blocks = _split_mol2_blocks(str(filepath))
        logger.info(f"  Found {len(blocks)} molecule blocks in {path.name}")

        records = []
        n_failed = 0

        for i, block in enumerate(blocks):
            mol = Chem.MolFromMol2Block(block, removeHs=True)

            if mol is None:
                n_failed += 1
                name = _extract_name_from_block(block)
                logger.debug(f"  Failed to parse MOL2 block {i+1}: "
                             f"'{name}'")
                continue

            smiles = Chem.MolToSmiles(mol)

            # Get name from block header
            mol_id = _extract_name_from_block(block)
            if not mol_id:
                mol_id = f"mol2_{i+1:04d}"

            records.append({"mol_id": mol_id, "smiles": smiles})

        if n_failed > 0:
            logger.warning(f"  {n_failed} molecules failed to parse in "
                           f"{path.name}")

        logger.info(f"  Extracted {len(records)} molecules from "
                    f"{path.name}")

    finally:
        rd_logger.setLevel(original_level)

    return pd.DataFrame(records)
