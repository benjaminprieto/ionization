"""
pdb_parser.py - PDB ligand extractor
======================================
Extracts non-standard residues (ligands) from PDB files and converts
them to SMILES via RDKit.

Filters out common solvents, ions, and buffer molecules that are not
drug-like.

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

# Common non-drug residues to skip
EXCLUDE_RESIDUES = {
    "HOH", "WAT", "DOD",           # Water
    "SO4", "PO4", "CL", "NA",      # Ions
    "GOL", "EDO", "PEG", "PGE",    # Solvents / buffers
    "ACE", "NME", "NH2",           # Capping groups
    "MES", "TRS", "HED", "EPE",    # Buffer molecules
    "DMS", "ACT", "BME", "IMD",    # Small additives
    "FMT", "NO3", "SCN", "CIT",    # Crystallization additives
    "UNX", "UNL",                   # Unknown atoms
}


def parse_pdb_file(filepath: str) -> pd.DataFrame:
    """
    Extract ligands from a PDB file and convert to SMILES.

    Reads HETATM records, groups by residue name + chain + sequence
    number, and attempts to convert each unique ligand to SMILES.

    Args:
        filepath: Path to .pdb file.

    Returns:
        DataFrame with columns: mol_id, smiles

    Raises:
        ImportError: If RDKit is not installed.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required to read PDB files. "
            "Install with: conda install -c conda-forge rdkit"
        )

    path = Path(filepath)
    logger.info(f"  Parsing PDB ligands: {path.name}")

    # Suppress RDKit warnings
    from rdkit import RDLogger
    rd_logger = RDLogger.logger()
    original_level = rd_logger.level()
    rd_logger.setLevel(RDLogger.ERROR)

    try:
        # Parse PDB to get unique HETATM residue names
        ligand_residues = set()

        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("HETATM"):
                    resname = line[17:20].strip()
                    if resname and resname not in EXCLUDE_RESIDUES:
                        ligand_residues.add(resname)

        if not ligand_residues:
            logger.info(f"  No ligands found in {path.name}")
            return pd.DataFrame(columns=["mol_id", "smiles"])

        logger.info(f"  Found {len(ligand_residues)} unique ligand "
                    f"residues: {sorted(ligand_residues)}")

        # Use RDKit to parse the full PDB and extract ligands
        mol = Chem.MolFromPDBFile(str(filepath), removeHs=True,
                                   sanitize=False)

        records = []
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                logger.debug("  Full PDB sanitization failed; "
                             "trying per-residue extraction")

            # Group atoms by residue
            residue_atoms = {}
            for atom in mol.GetAtoms():
                info = atom.GetPDBResidueInfo()
                if info is None:
                    continue
                resname = info.GetResidueName().strip()
                if resname in EXCLUDE_RESIDUES:
                    continue
                if resname not in ligand_residues:
                    continue
                chain = info.GetChainId().strip()
                resnum = info.GetResidueNumber()
                key = f"{resname}_{chain}_{resnum}"
                if key not in residue_atoms:
                    residue_atoms[key] = []
                residue_atoms[key].append(atom.GetIdx())

            # Extract each ligand as a fragment
            for lig_id, atom_indices in residue_atoms.items():
                try:
                    emol = Chem.RWMol(mol)
                    atoms_to_remove = sorted(
                        set(range(mol.GetNumAtoms())) -
                        set(atom_indices),
                        reverse=True
                    )
                    for idx in atoms_to_remove:
                        emol.RemoveAtom(idx)
                    try:
                        Chem.SanitizeMol(emol)
                        smiles = Chem.MolToSmiles(emol)
                        if smiles:
                            records.append({
                                "mol_id": lig_id,
                                "smiles": smiles
                            })
                    except Exception:
                        logger.debug(f"  Failed to sanitize: {lig_id}")
                except Exception:
                    logger.debug(f"  Failed to extract: {lig_id}")

        if not records:
            logger.warning(f"  Could not extract SMILES from any ligand "
                           f"in {path.name}")

        logger.info(f"  Extracted {len(records)} ligands from {path.name}")

    finally:
        rd_logger.setLevel(original_level)

    return pd.DataFrame(records)
