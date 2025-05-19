# utils/chem_utils.py

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw

# Conditional import for Draw module
try:
    from rdkit.Chem import Draw
    CAN_DRAW = True
except ImportError:
    CAN_DRAW = False

def draw_molecule_with_labels(smiles):
    """
    Draws a molecule with atom indices as labels.
    Returns PIL image if available, else None.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    # Add atom map numbers (indices)
    for atom in mol.GetAtoms():
        atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))

    # Only try to draw if Draw is available
    if CAN_DRAW:
        try:
            return Draw.MolToImage(mol, size=(300, 300))
        except:
            return None
    else:
        return None

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
