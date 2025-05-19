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

# Import prediction functions
from .predict import (
    predict_binding_affinity,
    predict_toxicity,
    check_lipinski,
    interpret_binding_affinity
)

# Import core chem utilities
from .chem_utils import canonicalize_smiles

# Import PDF report generator
try:
    from .pdf_utils import generate_pdf_report
    _PDF_AVAILABLE = True
except ImportError:
    generate_pdf_report = None
    _PDF_AVAILABLE = False

# Import SAS calculator
try:
    from .sas import calculate_sas
    _SAS_AVAILABLE = True
except ImportError:
    calculate_sas = None
    _SAS_AVAILABLE = False

# Import ADMEt data fetcher
try:
    from .admet import get_swissadmet_data, parse_swissadmet
    _ADMET_AVAILABLE = True
except ImportError:
    get_swissadmet_data = None
    parse_swissadmet = None
    _ADMET_AVAILABLE = False

# --- Optional: Expose __all__ dynamically ---
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
    'interpret_binding_affinity'

    'canonicalize_smiles',  # Always available
    'calculate_sas',        # Fallback or real
    'get_swissadmet_data',  # SwissADMET support
    'parse_swissadmet'
]

if _PDF_AVAILABLE:
    __all__.append('generate_pdf_report')