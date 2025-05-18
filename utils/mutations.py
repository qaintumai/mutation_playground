# utils/mutations.py
from rdkit import Chem
from rdkit.Chem import AllChem

def apply_smarts_reaction(smiles, smarts):
    rxn = AllChem.ReactionFromSmarts(smarts)
    if not rxn:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    products = rxn.RunReactants((mol,))
    if products and products[0]:
        product = products[0][0]
        try:
            Chem.SanitizeMol(product)
            return Chem.MolToSmiles(product)
        except:
            return None
    else:
        return None

def mutate_add_methyl(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    editable_mol = Chem.RWMol(mol)
    methyl_idx = editable_mol.AddAtom(Chem.Atom(6))
    editable_mol.AddBond(0, methyl_idx, order=Chem.BondType.SINGLE)
    mutated_mol = editable_mol.GetMol()
    try:
        Chem.SanitizeMol(mutated_mol)
        return Chem.MolToSmiles(mutated_mol)
    except:
        return None

def mutate_aromatic_methyl(smiles): return apply_smarts_reaction(smiles, '[cH:1]>>[c:1]C')
def mutate_aromatic_fluoro(smiles): return apply_smarts_reaction(smiles, '[cH:1]>>[c:1]F')
def mutate_aromatic_chloro(smiles): return apply_smarts_reaction(smiles, '[cH:1]>>[c:1]Cl')
def mutate_aromatic_methoxy(smiles): return apply_smarts_reaction(smiles, '[cH:1]>>[c:1]OC')
def mutate_bioisostere_cyan(smiles): return apply_smarts_reaction(smiles, '[C:1]=[O:2]>>[C:1][C:3]#[N:4]')

def get_possible_mutations(smiles):
    from rdkit import Chem
    from rdkit.Chem import Lipinski, rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return [], {}

    possible = []
    reasons = {}

    has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    has_carbonyl = any(
        atom.GetSymbol() == 'C' and any(
            bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetEndAtom().GetSymbol() == 'O'
            for bond in atom.GetBonds()
        )
        for atom in mol.GetAtoms()
    )

    # Fluorination
    if has_aromatic:
        possible.append("Fluorination")
        reasons["Fluorination"] = "Molecule contains aromatic ring(s)"
    else:
        reasons["Fluorination"] = "No aromatic H to replace"

    # Chlorination
    if has_aromatic:
        possible.append("Chlorination")
        reasons["Chlorination"] = "Molecule contains aromatic ring(s)"
    else:
        reasons["Chlorination"] = "No aromatic H to replace"

    # Methylation
    if has_aromatic:
        possible.append("Aromatic Methyl")
        reasons["Aromatic Methyl"] = "Molecule contains aromatic ring(s)"
    else:
        reasons["Aromatic Methyl"] = "No aromatic H to replace"

    # Methoxylation
    if has_aromatic:
        possible.append("Methoxylation")
        reasons["Methoxylation"] = "Molecule contains aromatic ring(s)"
    else:
        reasons["Methoxylation"] = "No aromatic H to replace"

    # Bioisostere (CN)
    if has_carbonyl:
        possible.append("Bioisostere (CN)")
        reasons["Bioisostere (CN)"] = "Molecule contains carbonyl group"
    else:
        reasons["Bioisostere (CN)"] = "No carbonyl group to replace"

    # General Methyl Addition
    if mol.GetNumAtoms() >= 1:
        possible.append("Add Methyl Group")
        reasons["Add Methyl Group"] = "Molecule has at least one carbon"
    else:
        reasons["Add Methyl Group"] = "Too small to add methyl group"

    return possible, reasons