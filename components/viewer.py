# components/viewer.py
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def generate_3d_html(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return ""
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        pass
    block = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=600, height=400)
    view.addModel(block, 'mol')
    view.setStyle({'stick': {}})
    view.zoomTo()
    return view._make_html()