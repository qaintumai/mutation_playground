# utils/admet.py

import requests
import time

def get_swissadmet_data(smiles):
    """
    Submit SMILES to SwissADMET and return parsed results.
    Returns None if failed.
    """
    url = "https://swissadmet.ch/predict "
    payload = {"smiles": smiles}

    try:
        response = requests.post(url, json=payload)
        if response.status_code == 200:
            data = response.json()
            # Wait a bit to avoid rate limiting
            time.sleep(1)
            return data
        else:
            print(f"Error fetching SwissADMET data: {response.status_code}")
            return None
    except Exception as e:
        print(f"Exception during SwissADMET request: {str(e)}")
        return None

def parse_swissadmet(data):
    if not data or "descriptors" not in data:
        return {}

    d = data["descriptors"]

    admet_data = {
        "Solubility (log mol/L)": d.get("Solubility"),
        "Permeability (Papp)": d.get("Papp"),
        "BBB Class": d.get("bbb_class"),
        "CYP Inhibition": ", ".join(d.get("CYP_inhibition", ["None"])),
        "T1/2 (h)": d.get("T1_2"),
        "Clearance (ml/min/kg)": d.get("Clearance"),
        "F(%)": d.get("F20_absorption"),
        "Hepatic Stability": d.get("HepaticStabilityClass"),
        "Renal Clearance": d.get("RenalClearanceClass"),
        "Plasma Protein Binding": d.get("PPBClass"),
    }

    return admet_data