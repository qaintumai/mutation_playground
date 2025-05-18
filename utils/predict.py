# utils/predict.py
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem
import math
import numpy as np
from data import calculate_physicochemical_properties

binding_model = joblib.load("models/binding_affinity_model.pkl")
toxicity_model = joblib.load("models/toxicity_classifier.pkl")

def featurize_molecule(mol):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
    return np.array(fp).reshape(1, -1)

def predict_toxicity(mol):
    from rdkit.Chem import Descriptors
    features = featurize_molecule(mol)
    tox_class = toxicity_model.predict(features)[0]
    proba = toxicity_model.predict_proba(features)[0]
    return {
        "class": "Toxic" if tox_class == 1 else "Non-Toxic",
        "confidence": max(proba)
    }

def predict_binding_affinity(mol):
    features = featurize_molecule(mol)
    delta_g = binding_model.predict(features)[0]
    return delta_g

def interpret_binding_affinity(delta_g):
    R_kcal = 1.987e-3
    T = 298.15
    try:
        Kd = math.exp(delta_g / (R_kcal * T))
        pKd = -math.log10(Kd)
    except:
        return {"pKd": "<0", "interpretation": "No meaningful binding"}
    if Kd > 1:
        interpretation = "Very weak or no binding"
    elif Kd > 0.001:
        interpretation = "Weak binding"
    elif Kd > 1e-6:
        interpretation = "Moderate binding"
    elif Kd > 1e-9:
        interpretation = "Strong binding"
    else:
        interpretation = "Very strong binding"
    return {
        "pKd": f"{pKd:.2f}",
        "interpretation": interpretation
    }

def check_lipinski(mol):
    from rdkit.Chem import Descriptors
    molwt = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    violations = sum([
        molwt > 500,
        logp > 5,
        h_donors > 10,
        h_acceptors > 20
    ])
    return {"Pass": violations <= 1}