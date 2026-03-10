"""
csv_report.py - CSV output
============================
Saves profiling results as a flat CSV table.
The simplest output — one row per molecule, one column per pH.

Project: ionprofile
"""

import logging
from pathlib import Path
from typing import Union

import pandas as pd

logger = logging.getLogger(__name__)


def save_csv(
        df: pd.DataFrame,
        output_path: Union[str, Path],
) -> str:
    """
    Save ionization profiling results as CSV.

    Args:
        df:          Results DataFrame.
        output_path: Output file path.

    Returns:
        Absolute path to saved file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_path, index=False, encoding="utf-8")
    logger.info(f"Saved CSV: {output_path} ({len(df)} rows)")

    return str(output_path.resolve())
