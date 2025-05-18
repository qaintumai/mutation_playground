# utils/chem_utils.py

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw

def draw_molecule_with_labels(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    for atom in mol.GetAtoms():
        atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))
    return Draw.MolToImage(mol)

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
