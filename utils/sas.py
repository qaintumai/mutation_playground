# utils/sascorer.py

from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
import logging
import os
import sys

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Try to load the RDKit SA_Score module
try:
    from rdkit.Chem.SA_Score.sascorer import calculateScore as rdkit_calculate_score
    HAS_RD_SASCORE = True
except ImportError:
    logger.warning("RDKit SA_Score not available. Using fallback scoring.")
    HAS_RD_SASCORE = False

def calculate_sas(mol):
    """
    Calculate Ertl's Synthetic Accessibility Score (SAS).
    Falls back to a simple heuristic if RDKit's SA_Score is unavailable.

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        float: Synthetic Accessibility Score (lower = more synthetic)
    """
    if mol is None:
        logger.error("Invalid molecule input for SAS calculation")
        return None

    if HAS_RD_SASCORE:
        try:
            return rdkit_calculate_score(mol)
        except Exception as e:
            logger.warning(f"Failed to compute SAS with RDKit: {str(e)}")

    # --- FALLBACK METHOD ---
    logger.info("Using fallback SAS method based on basic molecular properties")
    try:
        ring_count = sum(1 for x in mol.GetRingInfo().AtomRings() if len(x) > 6)
        sp3_ratio = rdMolDescriptors.CalcFractionCSP3(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        logp = Descriptors.MolLogP(mol)

        # Very rough approximation of SAS using known trends
        sas = (
            1.0 +
            0.2 * ring_count +       # Penalize large rings
            0.5 * sp3_ratio +         # SP3-rich = harder to synthesize
            0.1 * rot_bonds +         # More flexibility = harder
            0.3 * max(logp - 4, 0)     # Penalize high lipophilicity
        )
        return round(sas, 2)
    except Exception as e:
        logger.error(f"Fallback SAS failed: {str(e)}")
        return None