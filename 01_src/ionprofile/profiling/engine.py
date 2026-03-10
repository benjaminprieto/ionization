"""
engine.py - Ionization profiling orchestrator
===============================================
Main entry point for running ionization profiles. Coordinates:
    1. Input reading (via ionprofile.io)
    2. pH gradient generation
    3. Charge calculation + protonated SMILES capture (via engines)
    4. Transition metrics
    5. Output generation: tables (CSV, Excel, JSON) + structures (SDF)

This is the function you import when using ionprofile as a library:
    from ionprofile import run_profiling

Project: ionprofile
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from ionprofile.io.reader import read_molecules, read_directory
from ionprofile.profiling.ionizer import get_engine, check_dependencies
from ionprofile.reporting.csv_report import save_csv
from ionprofile.reporting.excel_report import save_excel
from ionprofile.reporting.json_report import save_json
from ionprofile.reporting.sdf_report import save_sdf

logger = logging.getLogger(__name__)


# =========================================================================
# pH GRADIENT
# =========================================================================

def generate_ph_values(
        ph_max: float,
        ph_min: float,
        ph_step: float,
) -> List[float]:
    """
    Generate descending pH values from ph_max to ph_min.

    Example:
        generate_ph_values(7.4, 6.0, 0.2)
        -> [7.4, 7.2, 7.0, 6.8, 6.6, 6.4, 6.2, 6.0]
    """
    values = []
    ph = ph_max
    while ph >= ph_min - (ph_step * 0.01):
        values.append(round(ph, 2))
        ph -= ph_step
    return values


def _ph_col_name(ph: float) -> str:
    """Charge column name: Q_pH74, Q_pH62, etc."""
    return f"Q_pH{int(round(ph * 10))}"


def _ph_smiles_col_name(ph: float) -> str:
    """Protonated SMILES column name: SMILES_pH74, SMILES_pH62, etc."""
    return f"SMILES_pH{int(round(ph * 10))}"


# =========================================================================
# IONIZATION CALCULATION
# =========================================================================

def calculate_ionization_profile(
        df: pd.DataFrame,
        ph_values: List[float],
        engine_name: str = "dimorphite",
        precision: float = 0.5,
) -> pd.DataFrame:
    """
    Calculate formal charge AND protonated SMILES at each pH.

    For each molecule x each pH, captures:
        - Q_pH{value}: formal charge (int)
        - SMILES_pH{value}: protonated SMILES (str)

    Also computes transition metrics:
        - N_Transitions: how many pH steps show a charge change
        - First_Transition_pH: pH where charge first differs from ph_max
    """
    engine = get_engine(engine_name)

    if not engine.is_available():
        logger.warning(f"Engine '{engine.name()}' not available. "
                       f"Charge columns will be null.")
        for ph in ph_values:
            df[_ph_col_name(ph)] = None
            df[_ph_smiles_col_name(ph)] = None
        df["N_Transitions"] = None
        df["First_Transition_pH"] = None
        return df

    # Silence RDKit warnings
    try:
        from rdkit import RDLogger
        RDLogger.logger().setLevel(RDLogger.ERROR)
    except ImportError:
        pass

    n_molecules = len(df)
    logger.info(f"Calculating ionization for {n_molecules} molecules "
                f"across {len(ph_values)} pH values...")
    logger.info(f"  Engine:   {engine.name()}")
    logger.info(f"  pH range: {ph_values[0]} -> {ph_values[-1]}")

    # Calculate charge AND protonated SMILES at each pH
    charge_data = {ph: [] for ph in ph_values}
    smiles_data = {ph: [] for ph in ph_values}

    for idx, row in df.iterrows():
        smiles = row.get("smiles", "")

        for ph in ph_values:
            q = engine.calculate_charge_at_ph(smiles, ph, precision)
            prot_smi = engine.get_protonated_smiles(smiles, ph, precision)
            charge_data[ph].append(q)
            smiles_data[ph].append(prot_smi)

        if (idx + 1) % 100 == 0 or (idx + 1) == n_molecules:
            logger.info(f"  Progress: {idx + 1}/{n_molecules}")

    # Add columns to DataFrame
    df = df.copy()
    for ph in ph_values:
        df[_ph_col_name(ph)] = charge_data[ph]
        df[_ph_smiles_col_name(ph)] = smiles_data[ph]

    # Add transition metrics
    _add_transition_metrics(df, ph_values)

    # Restore RDKit logger
    try:
        from rdkit import RDLogger
        RDLogger.logger().setLevel(RDLogger.WARNING)
    except ImportError:
        pass

    # Log summary
    col_max = _ph_col_name(ph_values[0])
    if col_max in df.columns:
        n_calc = df[col_max].notna().sum()
        n_neutral = (df[col_max] == 0).sum()
        logger.info(f"  Charges calculated: {n_calc}/{n_molecules}")
        logger.info(f"  Neutral at pH {ph_values[0]}: {n_neutral} "
                     f"({100 * n_neutral / max(n_molecules, 1):.1f}%)")

    return df


def _add_transition_metrics(
        df: pd.DataFrame,
        ph_values: List[float],
) -> None:
    """
    Add transition metrics (in-place).

    N_Transitions:      how many pH steps show a charge change.
    First_Transition_pH: pH where charge first differs from ph_max.
    """
    col_names = [_ph_col_name(ph) for ph in ph_values]

    n_transitions = []
    first_transition = []

    for _, row in df.iterrows():
        charges = [row.get(c) for c in col_names]

        transitions = 0
        first_ph = None

        if charges[0] is not None:
            ref_charge = charges[0]
            prev_charge = charges[0]

            for i in range(1, len(charges)):
                if charges[i] is not None and charges[i] != prev_charge:
                    transitions += 1
                    if first_ph is None and charges[i] != ref_charge:
                        first_ph = ph_values[i]
                    prev_charge = charges[i]

        n_transitions.append(
            transitions if charges[0] is not None else None
        )
        first_transition.append(first_ph)

    df["N_Transitions"] = n_transitions
    df["First_Transition_pH"] = first_transition


# =========================================================================
# MAIN PIPELINE
# =========================================================================

def run_profiling(
        input_path: Union[str, Path],
        output_dir: Union[str, Path],
        ph_max: float = 7.4,
        ph_min: float = 6.0,
        ph_step: float = 0.1,
        precision: float = 0.5,
        engine: str = "dimorphite",
        output_prefix: str = "ionization_profiling",
        output_formats: Optional[List[str]] = None,
        run_id: Optional[str] = None,
        smiles_column: Optional[str] = None,
        id_column: Optional[str] = None,
        format_hint: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run the complete ionization profiling pipeline.

    Args:
        input_path:     Path to input file or directory.
        output_dir:     Base output directory.
        ph_max:         Starting pH (default 7.4).
        ph_min:         Ending pH (default 6.0).
        ph_step:        Step between pH values (default 0.1).
        precision:      Engine precision parameter.
        engine:         Protonation engine: "dimorphite" (more coming).
        output_prefix:  Base name for output files.
        output_formats: List: ['csv', 'excel', 'json', 'sdf'].
                        Default: ['csv', 'json'].
        run_id:         Custom run identifier (default: timestamp).
        smiles_column:  Override SMILES column name (CSV input).
        id_column:      Override molecule ID column name (CSV input).
        format_hint:    Override input format auto-detection.

    Returns:
        Dict with: run_id, n_molecules, ph_values, output_dir,
                   output_files, dataframe, dependencies
    """
    # --- Setup ---
    if output_formats is None:
        output_formats = ["csv", "json"]

    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    input_path = Path(input_path)
    output_path = Path(output_dir) / run_id
    output_path.mkdir(parents=True, exist_ok=True)

    deps = check_dependencies(engine)

    # --- Header ---
    logger.info("=" * 60)
    logger.info("IONPROFILE v1.0.0")
    logger.info("=" * 60)
    logger.info(f"Timestamp:    {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Run ID:       {run_id}")
    logger.info(f"Input:        {input_path}")
    logger.info(f"Engine:       {engine}")
    logger.info(f"pH gradient:  {ph_max} -> {ph_min} (step {ph_step})")
    logger.info(f"Precision:    {precision}")
    logger.info(f"Output:       {output_path}")
    logger.info(f"Formats:      {output_formats}")
    logger.info(f"Dependencies: {deps}")
    logger.info("-" * 60)

    # --- Read input ---
    if input_path.is_dir():
        df = read_directory(input_path)
    else:
        kwargs = {}
        if smiles_column:
            kwargs["smiles_column"] = smiles_column
        if id_column:
            kwargs["id_column"] = id_column
        if format_hint:
            kwargs["format_hint"] = format_hint
        df = read_molecules(input_path, **kwargs)

    n_molecules = len(df)
    if n_molecules == 0:
        logger.warning("No molecules loaded. Check input file/directory.")
        return {
            "run_id": run_id,
            "n_molecules": 0,
            "ph_values": [],
            "output_dir": str(output_path),
            "output_files": {},
            "dataframe": df,
            "dependencies": deps,
        }

    logger.info(f"Loaded {n_molecules} molecules")

    # --- Generate pH gradient ---
    ph_values = generate_ph_values(ph_max, ph_min, ph_step)
    logger.info(f"pH values ({len(ph_values)}): {ph_values}")

    # --- Calculate ionization ---
    df = calculate_ionization_profile(
        df, ph_values,
        engine_name=engine,
        precision=precision,
    )

    # --- Generate outputs ---
    output_files = {}

    # Columns for table outputs (CSV, Excel, JSON)
    base_cols = ["mol_id", "smiles"]
    meta_cols = [c for c in ["source_file", "format"] if c in df.columns]
    charge_cols = [_ph_col_name(ph) for ph in ph_values
                   if _ph_col_name(ph) in df.columns]
    smiles_cols = [_ph_smiles_col_name(ph) for ph in ph_values
                   if _ph_smiles_col_name(ph) in df.columns]
    metric_cols = [c for c in ["N_Transitions", "First_Transition_pH"]
                   if c in df.columns]
    output_cols = (base_cols + meta_cols + charge_cols
                   + smiles_cols + metric_cols)
    df_out = df[output_cols].copy()

    if "csv" in output_formats:
        csv_path = output_path / f"{output_prefix}.csv"
        save_csv(df_out, csv_path)
        output_files["csv"] = str(csv_path)

    if "excel" in output_formats:
        xlsx_path = output_path / f"{output_prefix}.xlsx"
        save_excel(df_out, xlsx_path, ph_values)
        output_files["excel"] = str(xlsx_path)

    if "json" in output_formats:
        json_path = output_path / f"{output_prefix}.json"
        save_json(df_out, json_path, ph_values, {
            "run_id": run_id,
            "engine": engine,
            "ph_max": ph_max,
            "ph_min": ph_min,
            "ph_step": ph_step,
            "precision": precision,
            "input_path": str(input_path),
        })
        output_files["json"] = str(json_path)

    if "sdf" in output_formats:
        sdf_dir = output_path / "structures"
        sdf_files = save_sdf(df_out, sdf_dir, ph_values, output_prefix)
        output_files["sdf"] = sdf_files

    # --- Summary ---
    logger.info("")
    logger.info("=" * 60)
    logger.info("RESULTS")
    logger.info("=" * 60)
    logger.info(f"Molecules:    {n_molecules}")
    logger.info(f"Engine:       {engine}")
    logger.info(f"pH values:    {len(ph_values)}")
    for fmt, fpath in output_files.items():
        if isinstance(fpath, dict):
            for sub_key, sub_path in fpath.items():
                logger.info(f"  SDF {sub_key}: {sub_path}")
        else:
            logger.info(f"  {fmt.upper():10s}: {fpath}")
    logger.info("=" * 60)
    logger.info("IONPROFILE COMPLETE")
    logger.info("=" * 60)

    return {
        "run_id": run_id,
        "n_molecules": n_molecules,
        "ph_values": ph_values,
        "output_dir": str(output_path),
        "output_files": output_files,
        "dataframe": df_out,
        "dependencies": deps,
    }