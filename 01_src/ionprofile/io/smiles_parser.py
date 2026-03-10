"""
smiles_parser.py - SMILES file reader (CSV, TSV, TXT, SMI)
============================================================
Reads files containing SMILES strings in tabular or single-column format.
Auto-detects the SMILES and ID columns if not specified.

Supported layouts:
    - CSV/TSV with header: looks for columns named smiles/SMILES etc.
    - Single-column SMILES (one per line, .smi/.smiles)
    - Space/tab-separated: SMILES<sep>Name (standard .smi format)

Project: ionprofile
"""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# Common column names for SMILES
SMILES_COLUMN_NAMES = [
    "smiles", "SMILES", "Smiles",
    "canonical_smiles", "Canonical_SMILES",
    "SMILES_library", "SMILES_canonical",
    "smi", "SMI", "structure", "Structure",
    "isosmiles", "isomericSmiles",
]

# Common column names for molecule IDs
ID_COLUMN_NAMES = [
    "mol_id", "Name", "name", "ID", "id",
    "mol_name", "molecule_id", "compound_id",
    "title", "Title", "compound", "Compound",
    "mol_title",
]


def _detect_smiles_column(columns: list) -> Optional[str]:
    """Find the SMILES column by matching known names."""
    for candidate in SMILES_COLUMN_NAMES:
        if candidate in columns:
            return candidate
    return None


def _detect_id_column(columns: list) -> Optional[str]:
    """Find the ID column by matching known names."""
    for candidate in ID_COLUMN_NAMES:
        if candidate in columns:
            return candidate
    return None


def _detect_separator(filepath: str) -> str:
    """Detect CSV vs TSV by inspecting first lines."""
    with open(filepath, "r", encoding="utf-8") as f:
        first_lines = [f.readline() for _ in range(3)]

    text = "".join(first_lines)
    n_tabs = text.count("\t")
    n_commas = text.count(",")

    if n_tabs > n_commas:
        return "\t"
    return ","


def parse_smiles_file(
        filepath: str,
        smiles_column: Optional[str] = None,
        id_column: Optional[str] = None,
) -> pd.DataFrame:
    """
    Parse a file containing SMILES strings.

    Handles:
        - CSV/TSV with headers
        - .smi format (SMILES<space>name)
        - Single-column (one SMILES per line)

    Args:
        filepath:       Path to input file.
        smiles_column:  Explicit SMILES column name (auto-detect if None).
        id_column:      Explicit ID column name (auto-detect if None).

    Returns:
        DataFrame with columns: mol_id, smiles
    """
    path = Path(filepath)
    ext = path.suffix.lower()

    # Try .smi format first (no header, space-separated)
    if ext in (".smi", ".smiles"):
        return _parse_smi_format(filepath)

    # Read as tabular
    sep = _detect_separator(filepath)

    try:
        df = pd.read_csv(filepath, sep=sep, dtype=str)
    except Exception:
        return _parse_single_column(filepath)

    if len(df.columns) == 1 and df.columns[0] not in SMILES_COLUMN_NAMES:
        return _parse_single_column(filepath)

    # Detect SMILES column
    smi_col = smiles_column or _detect_smiles_column(list(df.columns))
    if smi_col is None:
        if len(df.columns) <= 2:
            smi_col = df.columns[0]
            logger.info(f"  Auto-selected first column as SMILES: "
                        f"'{smi_col}'")
        else:
            raise ValueError(
                f"Cannot detect SMILES column in {path.name}. "
                f"Columns found: {list(df.columns)}. "
                f"Use smiles_column= to specify."
            )

    # Detect ID column
    name_col = id_column or _detect_id_column(list(df.columns))

    # Build output
    result = pd.DataFrame()
    result["smiles"] = df[smi_col].astype(str).str.strip()

    if name_col and name_col in df.columns:
        result["mol_id"] = df[name_col].astype(str).str.strip()
    else:
        result["mol_id"] = [f"mol_{i+1:04d}" for i in range(len(df))]

    logger.info(f"  SMILES column: '{smi_col}', "
                f"ID column: '{name_col or 'auto-generated'}'")

    return result


def _parse_smi_format(filepath: str) -> pd.DataFrame:
    """Parse standard .smi format: SMILES<whitespace>Name per line."""
    records = []
    with open(filepath, "r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            mol_id = parts[1] if len(parts) > 1 else f"mol_{i:04d}"
            records.append({"mol_id": mol_id, "smiles": smiles})

    return pd.DataFrame(records)


def _parse_single_column(filepath: str) -> pd.DataFrame:
    """Parse file with one SMILES per line (no header, no ID)."""
    records = []
    with open(filepath, "r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            records.append({"mol_id": f"mol_{i:04d}", "smiles": line})

    return pd.DataFrame(records)
