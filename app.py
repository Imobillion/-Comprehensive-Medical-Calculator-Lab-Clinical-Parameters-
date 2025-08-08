# app.py
import streamlit as st
import pandas as pd
import numpy as np
from io import BytesIO, StringIO
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet

st.set_page_config(page_title="Comprehensive Medical Calculator", layout="wide")
st.title("🩺 Comprehensive Medical Calculator (Lab & Clinical Parameters)")

st.markdown("""
**Use:** Enter single patient data or upload CSV for batch.  
**Output:** Computed lab parameters + explanations and downloadable PDF report per patient.
""")

# ---------- Helper calculation functions ----------
def kg_from_lb(lb): return lb * 0.45359237
def cm_from_inch(inch): return inch * 2.54

# General body
def calc_bmi(weight_kg, height_cm):
    h_m = height_cm / 100.0
    return weight_kg / (h_m*h_m) if h_m > 0 else np.nan

def calc_bsa_mosteller(weight_kg, height_cm):
    return np.sqrt((height_cm * weight_kg) / 3600.0)

def calc_map(sbp, dbp):
    return (sbp + 2*dbp) / 3.0

# Hematology
def calc_mcv(hct_percent, rbc_million_per_uL):
    # hct in %, rbc in million/µL
    return (hct_percent * 10) / rbc_million_per_uL if rbc_million_per_uL else np.nan

def calc_mch(hb_g_dL, rbc_million_per_uL):
    return (hb_g_dL * 10) / rbc_million_per_uL if rbc_million_per_uL else np.nan

def calc_mchc(hb_g_dL, hct_percent):
    return (hb_g_dL * 100) / hct_percent if hct_percent else np.nan

# Lipid profile
def calc_lipids(tc, hdl, tg):
    # Friedewald if TG < 400
    vldl = tg / 5.0 if pd.notna(tg) else np.nan
    ldl = tc - hdl - vldl if pd.notna(tc) and pd.notna(hdl) and pd.notna(vldl) else np.nan
    ratios = {"TC/HDL": tc/hdl if pd.notna(tc) and pd.notna(hdl) and hdl!=0 else np.nan}
    return {"TotalChol": tc, "HDL": hdl, "Triglycerides": tg, "VLDL": vldl, "LDL": ldl, **ratios}

# Renal
def cockcroft_gault(age, weight_kg, scr_mg_dl, sex):
    if scr_mg_dl <= 0 or pd.isna(scr_mg_dl): return np.nan
    crcl = ((140 - age) * weight_kg) / (72 * scr_mg_dl)
    if str(sex).lower().startswith("f"): crcl *= 0.85
    return crcl

def egfr_ckd_epi(creatinine, age, sex, race_black=False):
    # very simplified CKD-EPI approximation (not clinical exactness), we include simple formula placeholder
    # For production use implement validated CKD-EPI equations with units and race consideration
    return np.nan

# Liver
def bilirubin_fraction(total, direct):
    indirect = total - direct if pd.notna(total) and pd.notna(direct) else np.nan
    return indirect

# Simple ESR note (ESR measured in lab; no validated simple formula from CBC)
def interpret_esr(esr_value, age, sex):
    # normal adult ESR rule-of-thumb: male: age/2, female: (age+10)/2
    normal = (age/2) if sex.lower().startswith("m") else ((age + 10)/2)
    return {"ESR": esr_value, "ExpectedUpper": normal, "Raised": bool(esr_value > normal)}

# Thyroid/hormone placeholders: we echo input and provide simple flags
def interpret_tft(tsh, t3, t4):
    out = {}
    out["TSH"] = tsh
    out["T3"] = t3
    out["T4"] = t4
    # very simplistic interpretation
    if pd.notna(tsh):
        if tsh > 4.0: out["TFT_Interpretation"] = "Suggestive of hypothyroid (confirm clinically)"
        elif tsh < 0.4: out["TFT_Interpretation"] = "Suggestive of hyperthyroid (confirm clinically)"
        else: out["TFT_Interpretation"] = "TSH within common reference range"
    return out

# Glucose interpretation
def interpret_glucose(fasting, random):
    out = {}
    if pd.notna(fasting):
        if fasting >= 126: out["Fasting"] = "Diabetic range (>=126 mg/dL)"
        elif fasting >= 100: out["Fasting"] = "Impaired fasting glucose (100-125 mg/dL)"
        else: out["Fasting"] = "Normal"
    if pd.notna(random):
        if random >= 200: out["Random"] = "Diabetic range (>=200 mg/dL)"
        elif random >= 140: out["Random"] = "Impaired (140-199 mg/dL)"
        else: out["Random"] = "Normal"
    return out

# ---------- PDF report generation ----------
def create_pdf_report(patient_info, results_dict):
    """
    patient_info: dict of patient metadata (Name, Age, Sex, ID)
    results_dict: ordered dict or dict of {parameter: (value, interpretation)}
    returns bytes of PDF
    """
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
    styles = getSampleStyleSheet()
    story = []

    title = Paragraph("Medical Test Report", styles["Title"])
    story.append(title)
    story.append(Spacer(1, 8))

    # patient info table
    info_data = [["Name", patient_info.get("Name", "")],
                 ["Age", str(patient_info.get("Age", ""))],
                 ["Sex", patient_info.get("Sex", "")],
                 ["ID", patient_info.get("ID", "")]]
    t = Table(info_data, hAlign="LEFT")
    t.setStyle(TableStyle([("BACKGROUND", (0,0),(0,-1), colors.lightgrey),
                           ("GRID", (0,0), (-1,-1), 0.25, colors.black)]))
    story.append(t)
    story.append(Spacer(1, 12))

    # results table
    table_data = [["Parameter", "Value", "Interpretation / Note"]]
    for param, item in results_dict.items():
        val = item.get("value", "")
        interp = item.get("note", "")
        table_data.append([param, str(val), interp])

    table = Table(table_data, colWidths=[160, 110, 240])
    table.setStyle(TableStyle([
        ("GRID", (0,0), (-1,-1), 0.25, colors.grey),
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#f0f0f0")),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ]))
    story.append(table)
    story.append(Spacer(1, 12))
    footer = Paragraph("This report is for informational purposes only. Not a substitute for clinical lab report or physician advice.", styles["Normal"])
    story.append(footer)

    doc.build(story)
    pdf = buffer.getvalue()
    buffer.close()
    return pdf

# ---------- UI: Input mode ----------
st.sidebar.header("Input / Upload")
mode = st.sidebar.radio("Mode", ["Single patient", "Batch CSV"])

def run_calculation_row(row):
    # row: dict-like pandas Series
    # standardize column keys lowercase
    r = {k.lower(): row.get(k) for k in row.index}
    # patient meta
    name = r.get("name", "")
    age = float(r.get("age", np.nan)) if pd.notna(r.get("age", np.nan)) else np.nan
    sex = str(r.get("sex", "Male"))
    # weight/height units handling
    weight = r.get("weight", np.nan)
    if pd.isna(weight): weight_kg = np.nan
    else:
        unit = str(r.get("weight_unit", "kg")).lower()
        weight_kg = float(weight) if unit=="kg" else kg_from_lb(float(weight))
    height = r.get("height", np.nan)
    if pd.isna(height): height_cm = np.nan
    else:
        hunit = str(r.get("height_unit", "cm")).lower()
        height_cm = float(height) if hunit=="cm" else cm_from_inch(float(height))

    # basic vitals
    sbp = r.get("systolic", np.nan); dbp = r.get("diastolic", np.nan)
    sbp = float(sbp) if pd.notna(sbp) else np.nan
    dbp = float(dbp) if pd.notna(dbp) else np.nan

    # CBC inputs
    hb = float(r.get("hb", r.get("hemoglobin", np.nan))) if pd.notna(r.get("hb", r.get("hemoglobin", np.nan))) else np.nan
    hct = float(r.get("hct", r.get("hematocrit", np.nan))) if pd.notna(r.get("hct", r.get("hematocrit", np.nan))) else np.nan
    rbc = float(r.get("rbc", np.nan)) if pd.notna(r.get("rbc", np.nan)) else np.nan
    # lipids
    tc = r.get("totalcholesterol", r.get("tc", np.nan)); hdl = r.get("hdl", np.nan); tg = r.get("triglycerides", r.get("tg", np.nan))
    tc = float(tc) if pd.notna(tc) else np.nan
    hdl = float(hdl) if pd.notna(hdl) else np.nan
    tg = float(tg) if pd.notna(tg) else np.nan
    # renal
    scr = r.get("serumcreatinine", r.get("scr", np.nan)); urea = r.get("urea", np.nan)
    scr = float(scr) if pd.notna(scr) else np.nan
    urea = float(urea) if pd.notna(urea) else np.nan
    # liver
    alt = r.get("alt", np.nan); ast = r.get("ast", np.nan); alp = r.get("alp", np.nan); tbil = r.get("tbil", np.nan); dbil = r.get("dbil", np.nan)
    alt = float(alt) if pd.notna(alt) else np.nan
    ast = float(ast) if pd.notna(ast) else np.nan
    alp = float(alp) if pd.notna(alp) else np.nan
    tbil = float(tbil) if pd.notna(tbil) else np.nan
    dbil = float(dbil) if pd.notna(dbil) else np.nan

    # thyroid
    tsh = r.get("tsh", np.nan); t3 = r.get("t3", np.nan); t4 = r.get("t4", np.nan)
    tsh = float(tsh) if pd.notna(tsh) else np.nan
    t3 = float(t3) if pd.notna(t3) else np.nan
    t4 = float(t4) if pd.notna(t4) else np.nan

    # glucose
    fasting = r.get("fasting_glucose", r.get("fasting", np.nan)); random_g = r.get("random_glucose", r.get("random", np.nan))
    fasting = float(fasting) if pd.notna(fasting) else np.nan
    random_g = float(random_g) if pd.notna(random_g) else np.nan

    results = {}

    # General body
    results["BMI"] = {"value": round(calc_bmi(weight_kg, height_cm),2) if pd.notna(weight_kg) and pd.notna(height_cm) else "", "note": "Body Mass Index"}
    results["BSA (Mosteller)"] = {"value": round(calc_bsa_mosteller(weight_kg, height_cm),3) if pd.notna(weight_kg) and pd.notna(height_cm) else "", "note": "Body Surface Area (m^2)"}
    results["MAP (mmHg)"] = {"value": round(calc_map(sbp, dbp),1) if pd.notna(sbp) and pd.notna(dbp) else "", "note": "Mean Arterial Pressure"}

    # Hematology
    results["Hemoglobin (g/dL)"] = {"value": hb if pd.notna(hb) else "", "note": "Hemoglobin"}
    results["Hematocrit (%)"] = {"value": hct if pd.notna(hct) else "", "note": "Hematocrit (packed cell volume)"}
    results["RBC (million/uL)"] = {"value": rbc if pd.notna(rbc) else "", "note": "Red Blood Cell count"}
    results["MCV (fL)"] = {"value": round(calc_mcv(hct, rbc),1) if pd.notna(hct) and pd.notna(rbc) else "", "note": "MCV = (Hct% * 10) / RBC(millions/µL)"}
    results["MCH (pg)"] = {"value": round(calc_mch(hb, rbc),1) if pd.notna(hb) and pd.notna(rbc) else "", "note": "MCH = (Hb * 10) / RBC"}
    results["MCHC (g/dL)"] = {"value": round(calc_mchc(hb, hct),1) if pd.notna(hb) and pd.notna(hct) else "", "note": "MCHC = (Hb*100)/Hct"}

    # ESR
    esr_in = r.get("esr", np.nan)
    if pd.notna(esr_in):
        esr_dict = interpret_esr(float(esr_in), age if pd.notna(age) else 40, sex)
        results["ESR (mm/hr)"] = {"value": esr_dict["ESR"], "note": f"Expected upper ~ {round(esr_dict['ExpectedUpper'],1)}; Raised: {esr_dict['Raised']}"}
    else:
        results["ESR (mm/hr)"] = {"value": "", "note": ""}

    # Lipids
    lip = calc_lipids(tc, hdl, tg)
    results["Total Cholesterol (mg/dL)"] = {"value": lip.get("TotalChol", ""), "note": ""}
    results["HDL (mg/dL)"] = {"value": lip.get("HDL",""), "note": "Higher HDL is protective"}
    results["Triglycerides (mg/dL)"] = {"value": lip.get("Triglycerides",""), "note": ""}
    results["LDL (mg/dL)"] = {"value": round(lip.get("LDL", np.nan),1) if pd.notna(lip.get("LDL", np.nan)) else "", "note": "Friedewald (if TG<400): LDL = TC - HDL - TG/5"}
    results["TC/HDL ratio"] = {"value": round(lip.get("TC/HDL", np.nan),2) if pd.notna(lip.get("TC/HDL", np.nan)) else "", "note": "Cardiovascular risk marker"}

    # Glucose
    ginterp = interpret_glucose(fasting, random_g)
    results["Fasting Glucose (mg/dL)"] = {"value": fasting if pd.notna(fasting) else "", "note": ginterp.get("Fasting","")}
    results["Random Glucose (mg/dL)"] = {"value": random_g if pd.notna(random_g) else "", "note": ginterp.get("Random","")}

    # Renal
    results["Serum Creatinine (mg/dL)"] = {"value": scr if pd.notna(scr) else "", "note": ""}
    results["Urea (mg/dL)"] = {"value": urea if pd.notna(urea) else "", "note": ""}
    crcl = cockcroft_gault(age if pd.notna(age) else 40, weight_kg if pd.notna(weight_kg) else 70, scr if pd.notna(scr) else np.nan, sex)
    results["Creatinine Clearance (mL/min, Cockcroft-Gault)"] = {"value": round(crcl,1) if pd.notna(crcl) else "", "note": "Use cautiously in extremes of muscle mass"}

    # Liver
    results["ALT (U/L)"] = {"value": alt if pd.notna(alt) else "", "note": ""}
    results["AST (U/L)"] = {"value": ast if pd.notna(ast) else "", "note": ""}
    results["ALP (U/L)"] = {"value": alp if pd.notna(alp) else "", "note": ""}
    results["Total Bilirubin (mg/dL)"] = {"value": tbil if pd.notna(tbil) else "", "note": ""}
    results["Direct Bilirubin (mg/dL)"] = {"value": dbil if pd.notna(dbil) else "", "note": ""}
    ind = bilirubin_fraction(tbil, dbil) if pd.notna(tbil) and pd.notna(dbil) else ""
    results["Indirect Bilirubin (mg/dL)"] = {"value": ind if pd.notna(ind) else "", "note": ""}

    # Thyroid
    tft = interpret_tft(tsh, t3, t4)
    results["TSH (µIU/mL)"] = {"value": tsh if pd.notna(tsh) else "", "note": tft.get("TFT_Interpretation","")}
    results["T3"] = {"value": t3 if pd.notna(t3) else "", "note": ""}
    results["T4"] = {"value": t4 if pd.notna(t4) else "", "note": ""}

    # Parasitology placeholder (user may input stool results or malaria RDT)
    stool = r.get("stool_ova", "")
    results["Stool Ova/Cysts"] = {"value": stool if stool!="" else "", "note": "Enter 'Positive' or 'Negative' or organism name if known"}
    malaria = r.get("malaria_rdt", "")
    results["Malaria RDT"] = {"value": malaria if malaria!="" else "", "note": "Positive/Negative"}

    return name, age, sex, results

# ---------- Single patient UI ----------
if mode == "Single patient":
    st.header("Single patient entry")
    with st.form("patient_form"):
        col1, col2, col3 = st.columns(3)
        with col1:
            name = st.text_input("Name")
            age = st.number_input("Age (years)", min_value=0, value=40)
            sex = st.selectbox("Sex", ["Male", "Female"])
            patient_id = st.text_input("Patient ID (optional)")
        with col2:
            weight = st.number_input("Weight", min_value=0.0, format="%.2f", value=70.0)
            weight_unit = st.selectbox("Weight unit", ["kg","lb"])
            height = st.number_input("Height", min_value=0.0, format="%.1f", value=170.0)
            height_unit = st.selectbox("Height unit", ["cm","inch"])
        with col3:
            hb = st.number_input("Hemoglobin (g/dL)", value=13.5, format="%.2f")
            hct = st.number_input("Hematocrit (%)", value=40.0, format="%.2f")
            rbc = st.number_input("RBC (million/µL)", value=4.5, format="%.2f")
        st.markdown("---")
        col4, col5, col6 = st.columns(3)
        with col4:
            sbp = st.number_input("Systolic BP (mmHg)", value=120)
            dbp = st.number_input("Diastolic BP (mmHg)", value=80)
            esr = st.number_input("ESR (mm/hr)", value=10)
        with col5:
            tc = st.number_input("Total Cholesterol (mg/dL)", value=180)
            hdl = st.number_input("HDL (mg/dL)", value=45)
            tg = st.number_input("Triglycerides (mg/dL)", value=120)
        with col6:
            scr = st.number_input("Serum creatinine (mg/dL)", value=1.0, format="%.2f")
            urea = st.number_input("Urea (mg/dL)", value=30.0, format="%.2f")
            alt = st.number_input("ALT (U/L)", value=25.0, format="%.1f")
        st.markdown("---")
        col7, col8 = st.columns(2)
        with col7:
            tsh = st.number_input("TSH (µIU/mL)", value=2.0, format="%.3f")
            t3 = st.number_input("T3", value=1.0, format="%.2f")
            t4 = st.number_input("T4", value=8.0, format="%.2f")
        with col8:
            fasting = st.number_input("Fasting Glucose (mg/dL)", value=90.0, format="%.1f")
            random_g = st.number_input("Random Glucose (mg/dL)", value=120.0, format="%.1f")
            stool_ova = st.text_input("Stool Ova/Cysts result (optional)")
        submitted = st.form_submit_button("Calculate")

    if submitted:
        # Build a Series-like row for reuse
        row = pd.Series({
            "name": name, "age": age, "sex": sex,
            "weight": weight, "weight_unit": weight_unit,
            "height": height, "height_unit": height_unit,
            "hb": hb, "hct": hct, "rbc": rbc,
            "systolic": sbp, "diastolic": dbp, "esr": esr,
            "totalcholesterol": tc, "hdl": hdl, "triglycerides": tg,
            "serumcreatinine": scr, "urea": urea, "alt": alt,
            "tsh": tsh, "t3": t3, "t4": t4,
            "fasting_glucose": fasting, "random_glucose": random_g,
            "stool_ova": stool_ova
        })
        name, age, sex, results = run_calculation_row(row)
        st.success("Calculations complete. See results below.")
        # show results in two columns
        left, right = st.columns(2)
        with left:
            st.subheader("Parameter — Value")
            for k,v in results.items():
                st.write(f"**{k}**: {v.get('value','')} — {v.get('note','')}")
        with right:
            st.subheader("Interpretation & Notes")
            for k,v in results.items():
                # keep short interpretive notes if any
                note = v.get("note","")
                st.write(f"**{k}** — {note}")

        # PDF report
        patient_info = {"Name": name, "Age": age, "Sex": sex, "ID": patient_id}
        pdf_bytes = create_pdf_report(patient_info, results)
        st.download_button("Download PDF Report", data=pdf_bytes, file_name=f"{name or 'patient'}_report.pdf", mime="application/pdf")

# ---------- Batch CSV mode ----------
else:
    st.header("Batch processing (CSV)")
    st.markdown("Upload CSV with columns: name, age, sex, weight, weight_unit, height, height_unit, hb, hct, rbc, systolic, diastolic, esr, totalcholesterol, hdl, triglycerides, serumcreatinine, urea, alt, tsh, t3, t4, fasting_glucose, random_glucose, stool_ova")
    uploaded = st.file_uploader("Upload CSV", type=["csv"])
    if uploaded:
        df_in = pd.read_csv(uploaded)
        st.info(f"Loaded {len(df_in)} rows")
        results_list = []
        pdfs = []
        for idx, row in df_in.iterrows():
            name, age, sex, results = run_calculation_row(row)
            # flatten results for CSV
            flat = {"name": name, "age": age, "sex": sex}
            for k,v in results.items():
                flat_key = k.replace(" ", "_").replace("/","_").replace("(","").replace(")","")
                flat[flat_key] = v.get("value","")
            results_list.append(flat)
            # create pdf per patient and store
            patient_info = {"Name": name, "Age": age, "Sex": sex, "ID": idx}
            pdf_bytes = create_pdf_report(patient_info, results)
            pdfs.append((f"{name or 'patient'}_{idx}.pdf", pdf_bytes))

        out_df = pd.DataFrame(results_list)
        st.subheader("Results preview")
        st.dataframe(out_df.head(), use_container_width=True)
        csv = out_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download results CSV", data=csv, file_name="batch_results.csv", mime="text/csv")

        # bundle PDFs into a zip for download (optional)
        import zipfile, os, tempfile
        tmp = BytesIO()
        with zipfile.ZipFile(tmp, "w") as z:
            for fname, b in pdfs:
                z.writestr(fname, b)
        tmp.seek(0)
        st.download_button("Download all PDF reports (zip)", data=tmp.getvalue(), file_name="patient_reports.zip", mime="application/zip")

st.markdown("---")
st.caption("Disclaimer: Calculations are for informational use only. Confirm with clinical lab and practitioner.")
