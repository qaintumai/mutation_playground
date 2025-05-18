# mutation_playground 🧬

🧪 AI-powered molecule mutation playground with predictive models
Predict binding affinity (ΔG), toxicity, synthetic accessibility score (SAS), Lipinski Rule of Five, and more.

---

## 🔍 Overview

This tool allows chemists and computational biologists to:
- Mutate molecules interactively
- Predict drug-likeness properties
- Filter promising candidates
- Visualize mutations in 2D/3D
- Export results for lab testing

Great for early-stage lead optimization or teaching medicinal chemistry concepts.

---

## 🚀 Features

| Feature | Description |
|--------|-------------|
| 🧪 Interactive Mutation Engine | Apply methyl, fluorine, chloro, methoxy, bioisostere changes |
| 💊 Property Prediction | Binding affinity, toxicity, Lipinski, SAScore |
| 🧮 Physicochemical Properties | MW, LogP, H-bond counts, TPSA |
| 🖼️ Molecule Viewer | 2D & 3D rendering
| 📦 Streamlit UI | Easy-to-use interface for scientists
| 🧠 Modular Design | Extendable for docking, ADME, or generative models

---

## 🛠 Requirements

- Python 3.9+
- RDKit
- Streamlit
- Numpy, Pandas, Scikit-learn

Install dependencies via:

```bash
pip install -r requirements.txt
```

## 🚀 Run Locally
```bash
streamlit run app.py
```

## 🧪 Testing
We use pytest for unit tests:
```bash
python -m pytest tests/
```

## 🤝 Contributing

PRs welcome! Especially for:

* Real docking integration (e.g., Smina)
* Better mutation logic
* SwissADMET or Tox21 integration
* PDF report generator
* Performance improvements