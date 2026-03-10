# ionprofile

**pH-dependent ionization profiling for drug discovery molecules.**

Calculates formal charge and protonation states across a pH gradient using multiple engines. Generates structured outputs (CSV, Excel, JSON, SDF) ready for downstream docking and ADMET workflows.

---

## Installation

```bash
# From GitHub
pip install git+https://github.com/benjaminprieto/ionization.git

# Or clone and install in editable mode
git clone https://github.com/benjaminprieto/ionization.git
cd ionization
pip install -e ".[all]"
```

### Dependencies

**Required:** pandas, numpy, pyyaml, openpyxl

**Chemistry (install via conda):**
```bash
conda install -c conda-forge rdkit openbabel
pip install dimorphite-dl qupkake
```

## Quick Start

```bash
# Simplest — point to your molecules
python 02_scripts/run_profiling.py 04_data/molecules/

# Choose engine
python 02_scripts/run_profiling.py 04_data/molecules/ -e openbabel

# Custom pH range
python 02_scripts/run_profiling.py 04_data/molecules/ --ph 7.4 5.0

# Single file
python 02_scripts/run_profiling.py my_compounds.sdf

# Full control (single line, works on Windows and Linux)
python 02_scripts/run_profiling.py 04_data/molecules/ -e qupkake --ph 7.4 6.0 --step 0.2 --name my_run --formats csv excel json sdf

# List available engines
python 02_scripts/run_profiling.py --engines
```

## Engines

| Engine | Speed | Accuracy | Best for |
|--------|-------|----------|----------|
| `dimorphite` | ms/mol | Empirical rules | Virtual screening |
| `openbabel` | ms/mol | ~30 SMARTS rules | Quick estimates |
| `qupkake` | min/mol | RMSE 0.5-0.8 pKa | Lead optimization |

```bash
python 02_scripts/run_profiling.py --engines   # check what's installed
```

## Outputs

Each run generates a folder in `05_results/{run_id}/` with:

| Format | Description |
|--------|-------------|
| CSV | Flat table: mol_id, smiles, Q_pH74, SMILES_pH74, ..., N_Transitions |
| Excel | Formatted workbook with color-coded charges |
| JSON | Structured metadata + per-molecule records |
| SDF | Individual 3D structures per molecule per pH (with explicit H) |

### SDF Output Structure

```
05_results/{run_id}/structures/
├── pH74/
│   ├── HTS1710-00236567-01.sdf
│   ├── HTS1710-00277847-01.sdf
│   └── ...
├── pH72/
│   └── ...
└── pH60/
    └── ...
```

Each SDF includes 3D coordinates, explicit hydrogens (antechamber-ready), and properties: `mol_id`, `pH`, `formal_charge`, `protonated_smiles`.

## Use as Library

```python
from ionprofile import run_profiling

result = run_profiling(
    input_path="04_data/molecules/",
    output_dir="05_results",
    ph_max=7.4,
    ph_min=6.0,
    engine="openbabel",
    output_formats=["csv", "excel", "json", "sdf"],
    run_id="my_experiment",
)

df = result["dataframe"]
print(f"Processed {result['n_molecules']} molecules")
```

### Use from another project (e.g. molecular_docking)

In `environment.yaml`:
```yaml
dependencies:
  - pip:
    - git+https://github.com/benjaminprieto/ionization.git
```

Then:
```python
from ionprofile import run_profiling

result = run_profiling(
    input_path="path/to/molecules/",
    output_dir="results/",
    ph_max=7.2,
    ph_min=6.2,
    engine="dimorphite",
)

# Use structures directly for docking
sdf_dir = result["output_files"]["sdf"]["pH72"]
# → "results/{run_id}/structures/pH72/"
```

## Project Structure

```
ionization/
├── 01_src/ionprofile/
│   ├── __init__.py
│   ├── io/                     # Format readers (SDF, MOL2, CSV, PDB)
│   │   ├── reader.py           # Auto-detect format + unified reader
│   │   ├── smiles_parser.py
│   │   ├── sdf_parser.py
│   │   ├── mol2_parser.py
│   │   └── pdb_parser.py
│   ├── profiling/
│   │   ├── engine.py           # Main orchestrator (run_profiling)
│   │   ├── ionizer.py          # Engine dispatcher
│   │   ├── rdkit_utils.py      # Neutralization, charge calculation
│   │   └── engines/
│   │       ├── base.py         # Abstract engine interface
│   │       ├── dimorphite_engine.py
│   │       ├── openbabel_engine.py
│   │       └── qupkake_engine.py
│   └── reporting/
│       ├── csv_report.py
│       ├── excel_report.py
│       ├── json_report.py
│       └── sdf_report.py       # Multi-mol + individual SDF output
├── 02_scripts/
│   └── run_profiling.py        # CLI entry point
├── 03_configs/
│   └── profiling.yaml          # Default configuration
├── 04_data/
│   ├── molecules/              # Input: SDF, MOL2, CSV, SMILES
│   └── proteins/               # PDB, mmCIF (future)
├── 05_results/                  # Output organized by run_id
├── 06_tests/
├── environment.yaml
├── pyproject.toml
└── README.md
```

## Key Features

- **Auto-format detection** — reads SDF, MOL2, CSV, SMILES, PDB
- **Pre-ionization neutralization** — handles charged SMILES from docking SDF files
- **Individual SDF per molecule** — one file per molecule per pH, with explicit H for antechamber
- **Multi-engine** — swap engines with `-e` flag, no code changes
- **Config-driven or CLI** — YAML for reproducibility, flags for quick runs
- **Installable** — `pip install` from GitHub, import as library

## Adding a New Engine

1. Create `engines/my_engine.py` inheriting from `BaseEngine`
2. Implement `calculate_charge_at_ph()` and `get_protonated_smiles()`
3. Add one line to `ionizer.py`:
```python
from ionprofile.profiling.engines.my_engine import MyEngine
ENGINE_REGISTRY["my_engine"] = MyEngine
```

## License

MIT
