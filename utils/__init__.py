# utils/__init__.py

# Import mutation functions from mutations.py
from .mutations import (
    apply_smarts_reaction,
    mutate_add_methyl,
    mutate_aromatic_methyl,
    mutate_aromatic_fluoro,
    mutate_aromatic_chloro,
    mutate_aromatic_methoxy,
    mutate_bioisostere_cyan,
    get_possible_mutations
)

# Import prediction functions from predict.py
from .predict import (
    predict_binding_affinity,
    predict_toxicity,
    check_lipinski,
    interpret_binding_affinity
)

# Import chem utils functions from chem_utils.py
from .chem_utils import (
    draw_molecule_with_labels,
    canonicalize_smiles,
)

from .pdf_utils import generate_pdf_report

from .sas import calculate_sas

# Import the function to get SwissADEMT from admet.py
from .admet import get_swissadmet_data, parse_swissadmet

# Optional: Expose __all__ for clarity
__all__ = [
    # Mutation functions
    'apply_smarts_reaction',
    'mutate_add_methyl',
    'mutate_aromatic_methyl',
    'mutate_aromatic_fluoro',
    'mutate_aromatic_chloro',
    'mutate_aromatic_methoxy',
    'mutate_bioisostere_cyan',
    'get_possible_mutations',

    # Prediction functions
    'predict_binding_affinity',
    'predict_toxicity',
    'check_lipinski',
    'interpret_binding_affinity',

    # Chem_utils
    'draw_molecule_with_labels'
    'canonicalize_smiles',

    # pdf utils
    'generate_pdf_report',

    # SAS
    'calculate_sas'
]