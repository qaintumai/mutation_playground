# data.py
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski

def calculate_physicochemical_properties(smiles):
    """
    Calculate detailed physicochemical properties of a molecule from its SMILES.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    props = {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "H-Bond Donors": Lipinski.NumHDonors(mol),
        "H-Bond Acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Formal Charge": Chem.GetFormalCharge(mol),
        "Number of Rings": rdMolDescriptors.CalcNumRings(mol),
        "Heavy Atom Count": rdMolDescriptors.CalcNumHeavyAtoms(mol),
        "Fraction Csp³": rdMolDescriptors.CalcFractionCSP3(mol),
    }

    return props