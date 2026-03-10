#!/usr/bin/env python3
"""
run_profiling.py - ionprofile CLI
===================================
Simple command-line interface for ionization profiling.

Usage:
    # Simplest: just point to your molecules
    python 02_scripts/run_profiling.py 04_data/molecules/

    # With engine choice
    python 02_scripts/run_profiling.py 04_data/molecules/ -e openbabel

    # Single file
    python 02_scripts/run_profiling.py my_compounds.sdf

    # Custom pH range
    python 02_scripts/run_profiling.py 04_data/molecules/ --ph 7.4 6.0

    # Full control (single line — works on Windows and Linux)
    python 02_scripts/run_profiling.py 04_data/molecules/ -e qupkake --ph 7.4 6.0 --step 0.2 --name my_run

    # Full control with YAML config
    python 02_scripts/run_profiling.py --config 03_configs/profiling.yaml --data-dir 04_data/molecules/ --run-id test_ob --engine openbabel --formats csv excel json sdf

 nohup python 02_scripts/run_profiling.py --data-dir 04_data/molecules/ --engine qupkake --run-id overnight_open --formats csv excel json sdf > overnight_qupkake.log 2>&1 &
Project: ionprofile
Version: 1.1.0
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

__version__ = "1.1.0"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s'
)
logger = logging.getLogger(__name__)


def setup_log_file(log_path: Path):
    """Add file handler to root logger."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(str(log_path), mode='w', encoding='utf-8')
    handler.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s'
    ))
    logging.getLogger().addHandler(handler)
    logger.info(f"Log file: {log_path}")


def main():
    parser = argparse.ArgumentParser(
        description=f'ionprofile v{__version__} — pH-dependent '
                    f'ionization profiling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s 04_data/molecules/                                          # all files, default engine
  %(prog)s 04_data/molecules/ -e openbabel                             # use OpenBabel
  %(prog)s 04_data/molecules/ -e qupkake --name my_run                 # use QupKake (slow, accurate)
  %(prog)s my_compounds.sdf                                            # single SDF file
  %(prog)s compounds.csv --ph 7.4 5.0                                  # custom pH range
  %(prog)s 04_data/molecules/ --name golgi_screen                      # custom run name
  %(prog)s 04_data/molecules/ --formats csv excel                      # choose outputs
  %(prog)s --config 03_configs/profiling.yaml                          # full YAML config
  %(prog)s --engines                                                   # list available engines

Full control (single line):
  %(prog)s 04_data/molecules/ -e qupkake --ph 7.4 6.0 --step 0.2 --name my_run --formats csv excel json sdf

With YAML config:
  %(prog)s --config 03_configs/profiling.yaml --data-dir 04_data/molecules/ --run-id test_ob --engine openbabel
        '''
    )

    # === THE ONLY REQUIRED ARGUMENT ===
    parser.add_argument(
        'input', nargs='?', default=None,
        help='Input file or directory with molecules '
             '(SDF, MOL2, CSV, SMILES, PDB)'
    )

    # === COMMON OPTIONS (short flags) ===
    parser.add_argument(
        '-e', '--engine', type=str, default='dimorphite',
        choices=['dimorphite', 'openbabel', 'qupkake'],
        help='Protonation engine (default: dimorphite)'
    )
    parser.add_argument(
        '--ph', type=float, nargs=2, metavar=('MAX', 'MIN'),
        default=[7.4, 6.0],
        help='pH range as MAX MIN (default: 7.4 6.0)'
    )
    parser.add_argument(
        '--step', type=float, default=0.1,
        help='pH step size (default: 0.1)'
    )
    parser.add_argument(
        '--name', '-n', type=str, default=None,
        help='Run name for output folder (default: auto timestamp)'
    )
    parser.add_argument(
        '--formats', '-f', type=str, nargs='+',
        default=['csv', 'excel', 'json', 'sdf'],
        choices=['csv', 'excel', 'json', 'sdf'],
        help='Output formats (default: all)'
    )
    parser.add_argument(
        '--output', '-o', type=str, default='05_results',
        help='Output base directory (default: 05_results)'
    )

    # === LESS COMMON OPTIONS ===
    parser.add_argument(
        '--config', '-c', type=str, default=None,
        help='YAML config file (overrides all other options)'
    )
    parser.add_argument(
        '--data-dir', '-d', type=str, default=None,
        help='Directory with molecular files (alias for positional input)'
    )
    parser.add_argument(
        '--run-id', type=str, default=None,
        help='Alias for --name (backwards compatibility)'
    )
    parser.add_argument(
        '--precision', type=float, default=0.5,
        help='Engine pKa precision (default: 0.5)'
    )
    parser.add_argument(
        '--smiles-col', type=str, default=None,
        help='SMILES column name for CSV input'
    )
    parser.add_argument(
        '--id-col', type=str, default=None,
        help='Molecule ID column name for CSV input'
    )

    # === UTILITY ===
    parser.add_argument(
        '--engines', action='store_true',
        help='List available engines and exit'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Show debug output'
    )
    parser.add_argument(
        '--version', action='version',
        version=f'ionprofile {__version__}'
    )

    args = parser.parse_args()

    # --- List engines ---
    if args.engines:
        engines = list_engines()
        print(f"\nionprofile v{__version__} — Available engines:\n")
        for key, info in engines.items():
            status = "✓ READY" if info["available"] else "✗ NOT INSTALLED"
            print(f"  {key:15s}  {info['name']:20s}  [{status}]")
        print()
        return 0

    # --- Verbose ---
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    # Defaults
    engine = args.engine
    ph_max, ph_min = args.ph[0], args.ph[1]
    ph_step = args.step
    precision = args.precision
    output_formats = args.formats
    output_dir = args.output
    output_prefix = "ionization_profiling"
    run_id = args.name or args.run_id
    input_path = args.input or args.data_dir

    # Config mode overrides everything
    if args.config:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        params = config.get('parameters', {})
        outputs = config.get('outputs', {})

        engine = params.get('engine', engine)
        ph_max = params.get('ph_max', ph_max)
        ph_min = params.get('ph_min', ph_min)
        ph_step = params.get('ph_step', ph_step)
        precision = params.get('precision', precision)
        output_prefix = outputs.get('output_prefix', output_prefix)
        output_formats = outputs.get('formats', output_formats)
        output_dir = outputs.get('output_dir', output_dir)

        # Config can specify input too
        if input_path is None:
            input_path = config.get('input', {}).get('path')

    # --- Validate input ---
    if input_path is None:
        parser.print_help()
        print("\nError: provide an input file/directory or use --config")
        return 1

    input_p = Path(input_path)
    if not input_p.exists():
        logger.error(f"Input not found: {input_path}")
        return 1

    if ph_max <= ph_min:
        logger.error(f"pH max ({ph_max}) must be > pH min ({ph_min})")
        return 1

    # --- Run ID ---
    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    # --- Log file ---
    log_dir = Path(output_dir) / run_id
    log_dir.mkdir(parents=True, exist_ok=True)
    setup_log_file(log_dir / 'ionprofile.log')

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
            smiles_column=args.smiles_col,
            id_column=args.id_col,
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