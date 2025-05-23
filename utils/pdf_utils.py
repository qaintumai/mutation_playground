# utils/pdf_utils.py
from fpdf import FPDF
import os

def generate_pdf_report(result, output_path="mutation_report.pdf"):
    # Create PDF with Unicode support
    pdf = FPDF()
    pdf.add_page()
    pdf.set_auto_page_break(auto=True, margin=15)

    # Register a Unicode-compatible font
    font_path = os.path.join(os.path.dirname(__file__), "fonts", "DejaVuSans-Bold.ttf")
    if not os.path.exists(font_path):
        raise FileNotFoundError(f"Font file not found: {font_path}")

    pdf.add_font("DejaVu", fname=font_path, uni=True)
    pdf.set_font("DejaVu", size=12)

    # Title with emoji
    pdf.set_font("DejaVu", 'B', 16)
    pdf.cell(0, 10, "🧬 Mutation Result Report", ln=True, align='C')
    pdf.ln(10)

    # Mutation Details
    pdf.set_font("DejaVu", 'B', 12)
    pdf.cell(0, 10, "🔬 Mutation Details", ln=True)
    pdf.set_font("DejaVu", size=12)
    pdf.cell(0, 8, f"Original SMILES: {result['Original']}", ln=True)
    pdf.cell(0, 8, f"Mutated SMILES: {result['mutated_smi']}", ln=True)
    pdf.ln(5)

    # Predicted Properties
    pdf.set_font("DejaVu", 'B', 12)
    pdf.cell(0, 10, "📊 Predicted Properties", ln=True)
    pdf.set_font("DejaVu", size=12)

    pdf.cell(0, 8, f"Binding Affinity (ΔG): {result['delta_g']:.2f} kcal/mol", ln=True)
    pdf.cell(0, 8, f"pKd (~pIC50): {result['interpretation']['pKd']}", ln=True)
    pdf.cell(0, 8, f"Toxicity Prediction: {result['tox_result'].get('class', 'Unknown')}", ln=True)
    sas = f"{result['sas_score']:.2f}" if result.get('sas_score') else "N/A"
    pdf.cell(0, 8, f"Synthetic Accessibility Score: {sas}", ln=True)
    lipinski = "Yes ✅" if result.get('lipinski_pass', False) else "No ❌"
    pdf.cell(0, 8, f"Lipinski Rule Pass: {lipinski}", ln=True)
    pdf.ln(5)

    # Physicochemical Properties
    props = result.get("physprop", {})
    if props:
        pdf.set_font("DejaVu", 'B', 12)
        pdf.cell(0, 10, "🧪 Physicochemical Properties", ln=True)
        pdf.set_font("DejaVu", size=12)
        for key, value in props.items():
            pdf.cell(0, 8, f"{key}: {value}", ln=True)
        pdf.ln(5)

    # Footer
    pdf.set_font("DejaVu", size=10, style="I")
    pdf.cell(0, 10, "Generated by Molecular Mutation Playground", align='C', ln=True)

    # Output
    pdf.output(output_path)
    return output_path