# molecule_playground/models/__init__.py

import pickle
import os

# Get the directory of this __init__.py file
MODEL_DIR = os.path.dirname(__file__)

# Load models
with open(os.path.join(MODEL_DIR, 'binding_affinity_model.pkl'), 'rb') as f:
    binding_affinity_model = pickle.load(f)

with open(os.path.join(MODEL_DIR, 'toxicity_classifier.pkl'), 'rb') as f:
    toxicity_model = pickle.load(f)

# Optional: Expose them explicitly
__all__ = ['binding_affinity_model', 'toxicity_classifier']