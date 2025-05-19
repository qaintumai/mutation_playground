# utils/chem_utils.py

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw

def canonicalize_smiles(smiles):
    """
    Convert SMILES to canonical form.
    Returns None if invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None
