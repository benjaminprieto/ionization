"""
json_report.py - Structured JSON output
=========================================
Saves profiling results as structured JSON with:
    - Metadata (run config, timestamps)
    - Per-molecule ionization data
    - Aggregate statistics

Useful for programmatic consumption by downstream pipelines.

Project: ionprofile
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _to_serializable(obj: Any) -> Any:
    """Recursively convert numpy/special types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_to_serializable(v) for v in obj]
    elif isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj)
    elif isinstance(obj, (np.bool_,)):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, set):
        return sorted(obj)
    elif pd.isna(obj):
        return None
    return obj


def save_json(
        df: pd.DataFrame,
        output_path: Union[str, Path],
        ph_values: List[float] = None,
        config: Dict[str, Any] = None,
) -> str:
    """
    Save ionization profiling results as structured JSON.

    Args:
        df:          Results DataFrame.
        output_path: Output file path (.json).
        ph_values:   List of pH values (for statistics).
        config:      Run configuration metadata.

    Returns:
        Absolute path to saved file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if config is None:
        config = {}
    if ph_values is None:
        ph_values = []

    data = {
        "ionprofile": {
            "version": "1.0.0",
            "timestamp": datetime.now().isoformat(),
            "configuration": config,
        },
        "statistics": _build_statistics(df, ph_values),
        "molecules": _build_molecule_records(df, ph_values),
    }

    data = _to_serializable(data)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    logger.info(f"Saved JSON: {output_path}")
    return str(output_path.resolve())


def _build_statistics(
        df: pd.DataFrame,
        ph_values: List[float],
) -> Dict[str, Any]:
    """Build aggregate statistics."""
    stats = {"n_molecules": len(df)}

    if not ph_values:
        return stats

    # Charge distributions per pH
    distributions = {}
    for ph in ph_values:
        col = f"Q_pH{int(round(ph * 10))}"
        if col in df.columns:
            series = df[col]
            distributions[f"pH_{ph:.1f}"] = {
                "negative": int((series < 0).sum()),
                "neutral": int((series == 0).sum()),
                "positive": int((series > 0).sum()),
                "na": int(series.isna().sum()),
            }
    stats["charge_distributions"] = distributions

    # Key populations
    if len(ph_values) >= 2:
        col_max = f"Q_pH{int(round(ph_values[0] * 10))}"
        col_min = f"Q_pH{int(round(ph_values[-1] * 10))}"
        if col_max in df.columns and col_min in df.columns:
            stats["key_populations"] = {
                f"neutral_at_pH_{ph_values[0]}": int(
                    (df[col_max] == 0).sum()),
                f"positive_at_pH_{ph_values[-1]}": int(
                    (df[col_min] > 0).sum()),
                "neutral_to_positive": int(
                    ((df[col_max] == 0) & (df[col_min] > 0)).sum()),
            }

    # Transitions
    if "N_Transitions" in df.columns:
        trans = df["N_Transitions"].dropna()
        if len(trans) > 0:
            stats["transitions"] = {
                "no_transition": int((trans == 0).sum()),
                "single_transition": int((trans == 1).sum()),
                "multiple_transitions": int((trans > 1).sum()),
            }

    return stats


def _build_molecule_records(
        df: pd.DataFrame,
        ph_values: List[float],
) -> List[Dict[str, Any]]:
    """Build per-molecule records."""
    records = []

    charge_cols = [f"Q_pH{int(round(ph * 10))}" for ph in ph_values]
    charge_cols = [c for c in charge_cols if c in df.columns]

    for _, row in df.iterrows():
        rec = {
            "mol_id": row.get("mol_id", ""),
            "smiles": row.get("smiles", ""),
        }

        if "source_file" in df.columns:
            rec["source_file"] = row.get("source_file", "")
        if "format" in df.columns:
            rec["format"] = row.get("format", "")

        # Charges per pH
        charges = {}
        for col in charge_cols:
            val = row.get(col)
            ph_label = col.replace("Q_pH", "")
            charges[f"pH_{int(ph_label)/10:.1f}"] = (
                int(val) if pd.notna(val) else None
            )
        rec["charges"] = charges

        # Transitions
        if "N_Transitions" in df.columns:
            val = row.get("N_Transitions")
            rec["n_transitions"] = int(val) if pd.notna(val) else None
        if "First_Transition_pH" in df.columns:
            val = row.get("First_Transition_pH")
            rec["first_transition_ph"] = float(val) if pd.notna(val) else None

        records.append(rec)

    return records
