"""
reader.py - Universal molecular input reader
==============================================
Auto-detects file format by extension and delegates to the appropriate
parser. Always returns a normalized DataFrame with columns:

    mol_id      : str   - unique molecule identifier
    smiles      : str   - canonical SMILES string
    source_file : str   - path to the originating file
    format      : str   - detected format (smiles_csv, sdf, mol2, pdb)

This is the ONLY entry point for data ingestion. The profiling engine
never needs to know what format the input was in.

Project: ionprofile
"""

import logging
from pathlib import Path
from typing import List, Optional, Union

import pandas as pd

from ionprofile.io.smiles_parser import parse_smiles_file
from ionprofile.io.sdf_parser import parse_sdf_file
from ionprofile.io.mol2_parser import parse_mol2_file
from ionprofile.io.pdb_parser import parse_pdb_file

logger = logging.getLogger(__name__)

# =========================================================================
# FORMAT REGISTRY
# =========================================================================

FORMAT_MAP = {
    ".csv":    ("smiles_csv", parse_smiles_file),
    ".tsv":    ("smiles_csv", parse_smiles_file),
    ".txt":    ("smiles_csv", parse_smiles_file),
    ".smi":    ("smiles_csv", parse_smiles_file),
    ".smiles": ("smiles_csv", parse_smiles_file),
    ".sdf":    ("sdf",        parse_sdf_file),
    ".mol":    ("sdf",        parse_sdf_file),
    ".mol2":   ("mol2",       parse_mol2_file),
    ".pdb":    ("pdb",        parse_pdb_file),
}


def detect_format(filepath: Union[str, Path]) -> str:
    """
    Detect molecular file format from extension.

    Returns:
        Format string: 'smiles_csv', 'sdf', 'mol2', 'pdb'

    Raises:
        ValueError: If extension is not recognized.
    """
    ext = Path(filepath).suffix.lower()
    if ext not in FORMAT_MAP:
        supported = ", ".join(sorted(FORMAT_MAP.keys()))
        raise ValueError(
            f"Unsupported file format: '{ext}'. "
            f"Supported extensions: {supported}"
        )
    return FORMAT_MAP[ext][0]


def read_molecules(
        filepath: Union[str, Path],
        format_hint: Optional[str] = None,
        smiles_column: Optional[str] = None,
        id_column: Optional[str] = None,
) -> pd.DataFrame:
    """
    Read molecules from any supported format.

    Args:
        filepath:       Path to input file.
        format_hint:    Override auto-detection ('smiles_csv', 'sdf',
                        'mol2', 'pdb').
        smiles_column:  For CSV/TSV: name of column containing SMILES.
                        Auto-detected if None.
        id_column:      For CSV/TSV: name of column for molecule IDs.
                        Auto-detected if None.

    Returns:
        DataFrame with columns: mol_id, smiles, source_file, format
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    # Detect or validate format
    if format_hint:
        fmt = format_hint
        parser_fn = None
        for ext_info in FORMAT_MAP.values():
            if ext_info[0] == fmt:
                parser_fn = ext_info[1]
                break
        if parser_fn is None:
            raise ValueError(f"Unknown format_hint: '{fmt}'")
    else:
        ext = filepath.suffix.lower()
        if ext not in FORMAT_MAP:
            supported = ", ".join(sorted(FORMAT_MAP.keys()))
            raise ValueError(
                f"Cannot auto-detect format for '{ext}'. "
                f"Supported: {supported}. "
                f"Use format_hint= to override."
            )
        fmt, parser_fn = FORMAT_MAP[ext]

    logger.info(f"Reading {filepath.name} as '{fmt}'...")

    # Dispatch to parser
    kwargs = {}
    if fmt == "smiles_csv":
        if smiles_column:
            kwargs["smiles_column"] = smiles_column
        if id_column:
            kwargs["id_column"] = id_column

    df = parser_fn(str(filepath), **kwargs)

    # Normalize output schema
    df["source_file"] = str(filepath)
    df["format"] = fmt

    # Ensure required columns
    for col in ["mol_id", "smiles"]:
        if col not in df.columns:
            raise ValueError(
                f"Parser for '{fmt}' did not produce required column: {col}"
            )

    # Drop rows with no SMILES
    n_before = len(df)
    df = df.dropna(subset=["smiles"])
    df = df[df["smiles"].str.strip() != ""]
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        logger.warning(f"  Dropped {n_dropped} entries with no SMILES")

    logger.info(f"  Loaded {len(df)} molecules from {filepath.name}")
    return df.reset_index(drop=True)


def read_directory(
        dirpath: Union[str, Path],
        extensions: Optional[List[str]] = None,
        recursive: bool = False,
        **kwargs,
) -> pd.DataFrame:
    """
    Read all molecular files from a directory.

    Args:
        dirpath:    Path to directory containing molecular files.
        extensions: Limit to these extensions (e.g. ['.sdf', '.csv']).
                    Default: all supported formats.
        recursive:  Search subdirectories.
        **kwargs:   Passed to read_molecules() for each file.

    Returns:
        Combined DataFrame with molecules from all files.
    """
    dirpath = Path(dirpath)
    if not dirpath.is_dir():
        raise NotADirectoryError(f"Not a directory: {dirpath}")

    if extensions is None:
        extensions = list(FORMAT_MAP.keys())

    # Collect files
    files = []
    pattern_fn = dirpath.rglob if recursive else dirpath.glob
    for ext in extensions:
        ext_clean = ext if ext.startswith(".") else f".{ext}"
        files.extend(pattern_fn(f"*{ext_clean}"))

    files = sorted(set(files))

    if not files:
        logger.warning(f"No molecular files found in {dirpath}")
        return pd.DataFrame(columns=["mol_id", "smiles", "source_file",
                                      "format"])

    logger.info(f"Found {len(files)} molecular files in {dirpath}")

    # Read and combine
    frames = []
    for f in files:
        try:
            df = read_molecules(f, **kwargs)
            frames.append(df)
        except Exception as e:
            logger.warning(f"  Skipped {f.name}: {e}")

    if not frames:
        return pd.DataFrame(columns=["mol_id", "smiles", "source_file",
                                      "format"])

    combined = pd.concat(frames, ignore_index=True)

    # Check for duplicate mol_ids
    dups = combined["mol_id"].duplicated()
    if dups.any():
        n_dups = dups.sum()
        logger.warning(f"  {n_dups} duplicate mol_id values found; "
                       f"appending source suffix")
        combined.loc[dups, "mol_id"] = (
            combined.loc[dups, "mol_id"] + "_" +
            combined.loc[dups, "source_file"].apply(
                lambda x: Path(x).stem
            )
        )

    logger.info(f"  Total: {len(combined)} molecules from {len(files)} files")
    return combined
