#!/usr/bin/env python3
"""
run_profiling.py - ionprofile CLI
===================================
Command-line interface for ionization profiling.

Two modes:
    1. Config mode (recommended):
       python 02_scripts/run_profiling.py \
           --config 03_configs/profiling.yaml \
           --data-dir 04_data/molecules/

    2. Direct mode:
       python 02_scripts/run_profiling.py \
           --input molecules.sdf \
           --ph-max 7.4 --ph-min 6.0

Project: ionprofile
Version: 1.0.0
"""

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

import yaml

# Add 01_src to path for standalone execution
sys.path.insert(0, str(Path(__file__).parent.parent / '01_src'))

from ionprofile.profiling.engine import run_profiling
from ionprofile.profiling.ionizer import list_engines

__version__ = "1.0.0"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s'
)
logger = logging.getLogger(__name__)


def load_yaml(path: str) -> dict:
    """Load YAML configuration file."""
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    """Add file handler to root logger."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(str(log_path), mode='w', encoding='utf-8')
    handler.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    handler.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s'
    ))
    logging.getLogger().addHandler(handler)
    logger.info(f"Log file: {log_path}")


def main():
    parser = argparse.ArgumentParser(
        description=f'ionprofile v{__version__} -- pH-dependent '
                    f'ionization profiling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Config mode (recommended)
  python 02_scripts/run_profiling.py \\
      --config 03_configs/profiling.yaml \\
      --data-dir 04_data/molecules/

  # Single file, direct mode
  python 02_scripts/run_profiling.py \\
      --input my_compounds.sdf \\
      --ph-max 7.4 --ph-min 6.0 --ph-step 0.1

  # Custom run ID and output formats
  python 02_scripts/run_profiling.py \\
      --config 03_configs/profiling.yaml \\
      --data-dir 04_data/molecules/ \\
      --run-id my_experiment_01 \\
      --formats csv excel json sdf

  # CSV with explicit columns
  python 02_scripts/run_profiling.py \\
      --input compounds.csv \\
      --smiles-column canonical_smiles \\
      --id-column compound_id

  # List available engines
  python 02_scripts/run_profiling.py --list-engines
        '''
    )

    # Config mode
    parser.add_argument('--config', '-c', type=str,
                        help='YAML config file (03_configs/profiling.yaml)')
    parser.add_argument('--data-dir', '-d', type=str,
                        help='Directory with molecular files')

    # Direct mode
    parser.add_argument('--input', '-i', type=str, default=None,
                        help='Input file (SDF, MOL2, CSV, SMILES, PDB)')

    # Overrides
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output directory override')
    parser.add_argument('--run-id', type=str, default=None,
                        help='Custom run identifier (default: timestamp)')
    parser.add_argument('--engine', type=str, default=None,
                        help='Protonation engine (default: dimorphite)')
    parser.add_argument('--ph-max', type=float, default=None,
                        help='Override max pH (e.g. 7.4)')
    parser.add_argument('--ph-min', type=float, default=None,
                        help='Override min pH (e.g. 6.0)')
    parser.add_argument('--ph-step', type=float, default=None,
                        help='Override pH step (e.g. 0.1)')
    parser.add_argument('--precision', type=float, default=None,
                        help='Override engine precision')
    parser.add_argument('--formats', type=str, nargs='+', default=None,
                        choices=['csv', 'excel', 'json', 'sdf'],
                        help='Output formats (default from config)')
    parser.add_argument('--smiles-column', type=str, default=None,
                        help='SMILES column name for CSV input')
    parser.add_argument('--id-column', type=str, default=None,
                        help='Molecule ID column name for CSV input')
    parser.add_argument('--format-hint', type=str, default=None,
                        choices=['smiles_csv', 'sdf', 'mol2', 'pdb'],
                        help='Override input format auto-detection')
    parser.add_argument('--log-level', type=str, default=None,
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='Log level')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable DEBUG logging')
    parser.add_argument('--list-engines', action='store_true',
                        help='List available protonation engines and exit')
    parser.add_argument('--version', action='version',
                        version=f'ionprofile {__version__}')

    args = parser.parse_args()

    # --- List engines ---
    if args.list_engines:
        engines = list_engines()
        print(f"\nionprofile v{__version__} - Available engines:\n")
        for key, info in engines.items():
            status = "READY" if info["available"] else "NOT INSTALLED"
            print(f"  {key:15s}  {info['name']:20s}  [{status}]")
        print()
        return 0

    # --- Log level ---
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.log_level:
        logging.getLogger().setLevel(
            getattr(logging, args.log_level, logging.INFO)
        )

    # =========================================================================
    # DEFAULTS
    # =========================================================================
    engine = "dimorphite"
    ph_max = 7.4
    ph_min = 6.0
    ph_step = 0.1
    precision = 0.5
    output_prefix = "ionization_profiling"
    output_formats = ["csv", "json"]
    output_dir = "05_results"
    log_level = "INFO"
    run_id = args.run_id

    # =========================================================================
    # LOAD CONFIG
    # =========================================================================
    if args.config:
        config = load_yaml(args.config)
        params = config.get('parameters', {})
        outputs = config.get('outputs', {})

        engine = params.get('engine', engine)
        ph_max = params.get('ph_max', ph_max)
        ph_min = params.get('ph_min', ph_min)
        ph_step = params.get('ph_step', ph_step)
        precision = params.get('precision', precision)
        log_level = params.get('log_level', log_level)

        output_prefix = outputs.get('output_prefix', output_prefix)
        output_formats = outputs.get('formats', output_formats)
        output_dir = outputs.get('output_dir', output_dir)

    # --- Input source ---
    if args.data_dir:
        input_path = args.data_dir
    elif args.input:
        input_path = args.input
    else:
        parser.error(
            "Requires either --data-dir (config mode) or --input "
            "(direct mode). Use --list-engines to see available engines."
        )

    # --- CLI overrides ---
    if args.engine is not None:
        engine = args.engine
    if args.ph_max is not None:
        ph_max = args.ph_max
    if args.ph_min is not None:
        ph_min = args.ph_min
    if args.ph_step is not None:
        ph_step = args.ph_step
    if args.precision is not None:
        precision = args.precision
    if args.formats is not None:
        output_formats = args.formats
    if args.output is not None:
        output_dir = args.output
    if args.log_level:
        log_level = args.log_level

    # =========================================================================
    # VALIDATE
    # =========================================================================
    input_p = Path(input_path)
    if not input_p.exists():
        logger.error(f"Input not found: {input_path}")
        return 1

    if ph_max <= ph_min:
        logger.error(f"ph_max ({ph_max}) must be > ph_min ({ph_min})")
        return 1

    if ph_step <= 0:
        logger.error(f"ph_step ({ph_step}) must be > 0")
        return 1

    # =========================================================================
    # SETUP LOG FILE
    # =========================================================================
    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    log_dir = Path(output_dir) / run_id
    log_dir.mkdir(parents=True, exist_ok=True)
    setup_log_file(log_dir / 'ionprofile.log', log_level)

    # =========================================================================
    # RUN
    # =========================================================================
    try:
        result = run_profiling(
            input_path=input_path,
            output_dir=output_dir,
            ph_max=ph_max,
            ph_min=ph_min,
            ph_step=ph_step,
            precision=precision,
            engine=engine,
            output_prefix=output_prefix,
            output_formats=output_formats,
            run_id=run_id,
            smiles_column=args.smiles_column,
            id_column=args.id_column,
            format_hint=args.format_hint,
        )

        logger.info(f"\nAll outputs in: {result['output_dir']}")
        return 0

    except Exception as e:
        logger.error(f"Failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
