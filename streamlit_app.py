# --- Molecular Mutation Playground ---
# AI-powered molecule mutation playground with predictive models

import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem

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

# Import utilities
from utils import (
    get_possible_mutations,
    mutate_add_methyl,
    mutate_aromatic_methyl,
    mutate_aromatic_fluoro,
    mutate_aromatic_chloro,
    mutate_aromatic_methoxy,
    mutate_bioisostere_cyan,
    predict_binding_affinity,
    predict_toxicity,
    interpret_binding_affinity,
    check_lipinski,
    canonicalize_smiles,
    calculate_sas,
    generate_pdf_report
)
from components.viewer import generate_3d_html
from data import calculate_physicochemical_properties

# --- Streamlit Page Setup ---
st.set_page_config(page_title="🧬 Molecular Mutation Playground", layout="wide")
st.title("🧬 AI-Powered Molecular Mutation Playground")

# Initialize session state
if "mutation_history" not in st.session_state:
    st.session_state.mutation_history = []

if "last_result" not in st.session_state:
    st.session_state.last_result = None

# Input SMILES
raw_smiles = st.text_input("Enter a base molecule (e.g., c1ccccc1):")
smiles_input = canonicalize_smiles(raw_smiles)

if not smiles_input:
    st.stop()

# Get possible mutations
possible_mutations, reasons = get_possible_mutations(smiles_input)

# Define mutation explanations
mutation_explanations = {
    "Fluorination": "Adds fluorine to aromatic rings. Most effective for improving metabolic stability.",
    "Chlorination": "Replaces aromatic H with chlorine. Useful for modulating lipophilicity and binding.",
    "Aromatic Methyl": "Adds methyl group to aromatic rings. Known as 'magic methyl' effect.",
    "Methoxylation": "Adds methoxy group (-OCH₃) to aromatic rings. Improves solubility and binding.",
    "Bioisostere (CN)": "Replaces carbonyl groups with nitrile. Often improves metabolic stability.",
    "Add Methyl Group": "Adds methyl to any carbon. Enhances potency or blocks metabolism."
}

# Map mutation names to functions
mutation_functions = {
    "Add Methyl Group": mutate_add_methyl,
    "Aromatic Methyl": mutate_aromatic_methyl,
    "Fluorination": mutate_aromatic_fluoro,
    "Chlorination": mutate_aromatic_chloro,
    "Methoxylation": mutate_aromatic_methoxy,
    "Bioisostere (CN)": mutate_bioisostere_cyan
}

# Display Mutations in Rows
st.markdown("### 🧪 Choose a Mutation to Apply")

for i, (mutation_name, mutator_func) in enumerate(mutation_functions.items()):
    col1, col2 = st.columns([3, 7])

    with col1:
        mutated_smi = mutator_func(smiles_input)
        mol_valid = False
        if mutated_smi:
            mol = Chem.MolFromSmiles(mutated_smi)
            if mol:
                mol_valid = True

        explanation = mutation_explanations.get(mutation_name, "")
        reason = reasons.get(mutation_name, "No reason available.")

        if mol_valid:
            if st.button(f"🧪 {mutation_name}", key=f"mut_{i}"):
                mutated_smi = mutator_func(smiles_input)

                # Validate molecule
                mol = Chem.MolFromSmiles(mutated_smi)
                if mol:
                    try:
                        Chem.SanitizeMol(mol)
                    except Exception as e:
                        st.warning(f"⚠️ Molecule failed sanitization: {str(e)}")
                        mol = None

                if not mol:
                    st.error("❌ Cannot run predictions on invalid molecule.")
                    st.stop()

                # Run predictions
                tox_result = predict_toxicity(mol)
                delta_g = predict_binding_affinity(mol)
                interpretation = interpret_binding_affinity(delta_g)
                lipinski_pass = check_lipinski(mol)["Pass"]
                sas_score = calculate_sas(mol)

                # Save to history
                st.session_state.mutation_history.append({
                    "Original": smiles_input,
                    "Mutated": mutated_smi,
                    "Mutation Type": mutation_name,
                    "Binding Affinity (ΔG)": round(delta_g, 2),
                    "pKd (~pIC50)": interpretation["pKd"],
                    "Toxicity Prediction": tox_result["class"],
                    "Lipinski Rule Pass": lipinski_pass,
                    "SAS Score": round(sas_score, 2) if sas_score else None,
                })

                # Update last result
                st.session_state.last_result = {
                    "mutated_smi": mutated_smi,
                    "tox_result": tox_result,
                    "delta_g": delta_g,
                    "interpretation": interpretation,
                    "lipinski_pass": lipinski_pass,
                    "sas_score": sas_score,
                }

        else:
            st.button(f"❌ {mutation_name}", disabled=True, help=reason)

    with col2:
        st.markdown(f"**{mutation_name}:** {explanation}")
        if not mol_valid:
            st.caption(f"⚠️ Why it's unavailable: {reason}")

# --- Mutation Result Section (Full Width) ---
st.markdown("---")
if st.session_state.last_result:
    result = st.session_state.last_result
    mutated_smi = result["mutated_smi"]
    st.subheader("🧬 Mutation Result")

    col1, col2 = st.columns([3, 3])

    with col1:
        st.markdown("##### Original Molecule")
        # Generate full HTML content including wrapper and SMILES
        html_content = f"""
        <div style="border:1px solid #ccc; padding:10px; border-radius:8px;">
            {generate_3d_html(smiles_input)}
            <p><strong>SMILES:</strong> <code>{smiles_input}</code></p>
        </div>
        """

        # Render it safely
        st.components.v1.html(html_content, height=500)

    with col2:
        st.markdown("##### Mutated Molecule")
        # Generate full HTML content including wrapper and SMILES
        html_content = f"""
        <div style="border:1px solid #ccc; padding:10px; border-radius:8px;">
            {generate_3d_html(mutated_smi)}
            <p><strong>SMILES:</strong> <code>{mutated_smi}</code></p>
        </div>
        """

        # Render it safely
        st.components.v1.html(html_content, height=500)


    # --- Predicted Properties Section ---
    st.markdown("### 📊 Predicted Properties")

    property_data = {
        "Property": [
            "Binding Affinity (ΔG)",
            "pKd (~pIC50)",
            "Toxicity Prediction",
            "Synthetic Accessibility Score",
            "Lipinski Rule Pass"
        ],
        "Value": [
            f"{result['delta_g']:.2f} kcal/mol" if result['delta_g'] is not None else "N/A",
            result['interpretation'].get('pKd', 'N/A'),
            result['tox_result'].get('class', 'Unknown'),
            f"{result['sas_score']:.2f}" if result['sas_score'] is not None else "N/A",
            "Yes ✅" if result.get('lipinski_pass', False) else "No ❌"
        ]
    }

    property_df = pd.DataFrame(property_data)
    st.table(property_df)

    # --- Physicochemical Properties Section ---
    st.markdown("### 🧪 Physicochemical Properties")
    props = calculate_physicochemical_properties(mutated_smi)
    prop_df = pd.DataFrame([props]).T.rename(columns={0: "Value"})
    st.dataframe(prop_df, use_container_width=True)

else:
    st.info("Click a mutation button to see results here.")

# --- Mutation History + Filtering ---
if st.session_state.mutation_history:
    st.subheader("📜 Mutation History")
    df = pd.DataFrame(st.session_state.mutation_history)
    st.dataframe(df)

    # Filter promising candidates
    filtered = df[df["Toxicity Prediction"] == "Non-Toxic"]
    filtered = filtered.copy()
    filtered["pKd (~pIC50)"] = pd.to_numeric(filtered["pKd (~pIC50)"], errors='coerce')
    filtered = filtered[filtered["pKd (~pIC50)"] > 5]
    filtered = filtered[filtered["Lipinski Rule Pass"]]
    filtered = filtered[filtered["SAS Score"].notna()]

    if not filtered.empty:
        st.success(f"✅ Found {len(filtered)} promising candidate(s):")
        st.dataframe(filtered)

        @st.cache_data
        def convert_df(dataframe):
            return dataframe.to_csv(index=False).encode('utf-8')

        csv = convert_df(filtered)
        st.download_button(
            label="📥 Download Filtered Candidates",
            data=csv,
            file_name='filtered_mutants.csv',
            mime='text/csv'
        )