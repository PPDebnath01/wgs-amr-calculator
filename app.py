import math
from io import BytesIO
from datetime import date

import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="WGS-AMR Sample Size Calculator",
    page_icon="🧬",
    layout="wide",
)

APP_VERSION = "1.0"
APP_DATE = "2026-05-30"


def ceil_safe(x: float) -> int:
    if x is None or not math.isfinite(x):
        return 0
    return int(math.ceil(x))


def detection_outputs(p: float, se: float, c: float, m: float, rho: float, f: float) -> dict:
    p_eff = p * se
    n0 = math.log(1 - c) / math.log(1 - p) if 0 < p < 1 and 0 < c < 1 else float("nan")
    n_se = math.log(1 - c) / math.log(1 - p_eff) if 0 < p_eff < 1 and 0 < c < 1 else float("nan")
    deff = 1 + (m - 1) * rho
    n_cluster = n_se * deff
    n_final = n_cluster / (1 - f) if f < 1 else float("nan")
    n_final_ceil = ceil_safe(n_final)
    effective_usable_n = n_final_ceil * (1 - f) / deff if deff > 0 else float("nan")
    p_detect = 1 - (1 - p_eff) ** effective_usable_n if 0 <= p_eff < 1 else float("nan")
    return {
        "Baseline detection n": ceil_safe(n0),
        "Effective prevalence": p_eff,
        "Sensitivity-adjusted n": ceil_safe(n_se),
        "Design effect (DEFF)": deff,
        "Cluster-adjusted n": ceil_safe(n_cluster),
        "Final integrated n": n_final_ceil,
        "Expected achieved detection probability": p_detect,
    }


def prevalence_outputs(p_prev: float, z: float, d: float, f_prev: float) -> dict:
    n_raw = (z ** 2) * p_prev * (1 - p_prev) / (d ** 2) if d > 0 else float("nan")
    n_final = n_raw / (1 - f_prev) if f_prev < 1 else float("nan")
    return {"Prevalence estimation sample size": ceil_safe(n_final)}


def association_outputs(p0: float, odds_ratio: float, z_alpha: float, z_beta: float) -> dict:
    p1 = (odds_ratio * p0) / (1 - p0 + odds_ratio * p0) if 0 < p0 < 1 and odds_ratio > 0 else float("nan")
    denominator = (p1 - p0) ** 2
    n_group = ((z_alpha + z_beta) ** 2) * (p0 * (1 - p0) + p1 * (1 - p1)) / denominator if denominator > 0 else float("nan")
    return {
        "Case gene frequency": p1,
        "Sample size per group": ceil_safe(n_group),
        "Total association sample size": 2 * ceil_safe(n_group),
    }


def outputs_dataframe(results: dict) -> pd.DataFrame:
    rows = []
    for key, value in results.items():
        if isinstance(value, float):
            display = f"{value:.4f}" if abs(value) < 10 else f"{value:.2f}"
        else:
            display = value
        rows.append({"Output": key, "Value": display})
    return pd.DataFrame(rows)


def make_excel_download(detection, prevalence, association, scenarios):
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
        outputs_dataframe(detection).to_excel(writer, index=False, sheet_name="Detection outputs")
        outputs_dataframe(prevalence).to_excel(writer, index=False, sheet_name="Prevalence outputs")
        outputs_dataframe(association).to_excel(writer, index=False, sheet_name="Association outputs")
        scenarios.to_excel(writer, index=False, sheet_name="Scenario analysis")
    return buffer.getvalue()


st.title("WGS-based AMR Surveillance Sample Size Calculator")
st.caption(f"Version {APP_VERSION} · validated against the Excel calculator equations · {APP_DATE}")

with st.sidebar:
    st.header("Detection inputs")
    p = st.number_input("Expected prevalence, p", min_value=0.001, max_value=0.5, value=0.05, step=0.001, format="%.3f")
    se = st.number_input("Genomic detection sensitivity, Se", min_value=0.50, max_value=1.00, value=0.80, step=0.01, format="%.2f")
    c = st.number_input("Target detection confidence, C", min_value=0.80, max_value=0.99, value=0.95, step=0.01, format="%.2f")
    m = st.number_input("Average cluster size, m", min_value=1.0, max_value=50.0, value=8.0, step=1.0)
    rho = st.number_input("Intracluster correlation coefficient, ρ", min_value=0.0, max_value=0.5, value=0.10, step=0.01, format="%.2f")
    f = st.number_input("Sequencing failure rate, f", min_value=0.0, max_value=0.5, value=0.10, step=0.01, format="%.2f")

results = detection_outputs(p, se, c, m, rho, f)

k1, k2, k3, k4 = st.columns(4)
k1.metric("Baseline n", results["Baseline detection n"])
k2.metric("Se-adjusted n", results["Sensitivity-adjusted n"])
k3.metric("DEFF", f"{results['Design effect (DEFF)']:.2f}")
k4.metric("Final integrated n", results["Final integrated n"])

st.progress(min(max(results["Expected achieved detection probability"], 0), 1), text=f"Expected achieved detection probability: {results['Expected achieved detection probability']:.1%}")

tab1, tab2, tab3, tab4, tab5 = st.tabs(["Detection calculator", "Prevalence estimation", "Association", "Scenario analysis", "User guide"])

with tab1:
    st.subheader("Detection-based WGS/metagenomic AMR surveillance")
    st.dataframe(outputs_dataframe(results), use_container_width=True, hide_index=True)
    st.markdown(
        """
        **Equation notes**

        Baseline detection: `n = ln(1−C) / ln(1−p)`  
        Effective prevalence: `p_eff = p × Se`  
        Sensitivity-adjusted detection: `n_Se = ln(1−C) / ln(1−p × Se)`  
        Design effect: `DEFF = 1 + (m−1)ρ`  
        Final integrated sample size: `n_final = n_Se × DEFF × 1/(1−f)`
        """
    )

with tab2:
    st.subheader("Precision-based prevalence estimation")
    c1, c2, c3, c4 = st.columns(4)
    p_prev = c1.number_input("Expected prevalence", min_value=0.001, max_value=0.5, value=0.10, step=0.001, format="%.3f")
    z = c2.selectbox("Z-score", options=[1.64, 1.96, 2.58], index=1)
    d = c3.number_input("Margin of error, d", min_value=0.01, max_value=0.10, value=0.05, step=0.01, format="%.2f")
    f_prev = c4.number_input("Failure rate", min_value=0.0, max_value=0.5, value=0.10, step=0.01, format="%.2f")
    prev_results = prevalence_outputs(p_prev, z, d, f_prev)
    st.metric("Prevalence estimation sample size", prev_results["Prevalence estimation sample size"])
    st.markdown("Equation: `n = Z²p(1−p)/d²`, inflated by `1/(1−f)` for unusable genomes.")

with tab3:
    st.subheader("Genotype–phenotype association")
    c1, c2, c3, c4 = st.columns(4)
    p0 = c1.number_input("Control gene frequency, p0", min_value=0.01, max_value=0.80, value=0.20, step=0.01, format="%.2f")
    odds_ratio = c2.number_input("Assumed odds ratio", min_value=1.01, max_value=20.0, value=2.0, step=0.1, format="%.2f")
    z_alpha = c3.number_input("Z alpha/2", min_value=1.0, max_value=3.5, value=1.96, step=0.01, format="%.2f")
    z_beta = c4.selectbox("Z beta / power", options=[0.84, 1.28], index=0, help="0.84 ≈ 80% power; 1.28 ≈ 90% power")
    assoc_results = association_outputs(p0, odds_ratio, z_alpha, z_beta)
    st.dataframe(outputs_dataframe(assoc_results), use_container_width=True, hide_index=True)
    st.markdown("Case gene frequency: `p1 = OR × p0 / (1 − p0 + OR × p0)`. Sample size is calculated per group and doubled for total sample size.")

with tab4:
    st.subheader("Scenario analysis")
    base_scenarios = [
        ["Tilapia aquaculture base case", 0.05, 0.80, 0.95, 8, 0.10, 0.10],
        ["Lower prevalence", 0.02, 0.80, 0.95, 8, 0.10, 0.10],
        ["Higher sensitivity", 0.05, 0.95, 0.95, 8, 0.10, 0.10],
        ["Increased clustering", 0.05, 0.80, 0.95, 10, 0.15, 0.10],
        ["Lower failure rate", 0.05, 0.80, 0.95, 8, 0.10, 0.05],
        ["Hospital effluent metagenomics", 0.10, 0.70, 0.95, 4, 0.20, 0.15],
        ["Current sidebar inputs", p, se, c, m, rho, f],
    ]
    scenario_rows = []
    for name, sp, sse, sc, sm, srho, sf in base_scenarios:
        out = detection_outputs(sp, sse, sc, sm, srho, sf)
        scenario_rows.append({
            "Scenario": name,
            "p": sp,
            "Se": sse,
            "Confidence": sc,
            "m": sm,
            "ICC/rho": srho,
            "Failure": sf,
            "Baseline n": out["Baseline detection n"],
            "Se-adjusted n": out["Sensitivity-adjusted n"],
            "DEFF": round(out["Design effect (DEFF)"], 3),
            "Final n": out["Final integrated n"],
            "Achieved probability": round(out["Expected achieved detection probability"], 4),
        })
    scenarios = pd.DataFrame(scenario_rows)
    edited_scenarios = st.data_editor(scenarios, use_container_width=True, hide_index=True, disabled=["Baseline n", "Se-adjusted n", "DEFF", "Final n", "Achieved probability"])
    st.bar_chart(edited_scenarios.set_index("Scenario")[["Baseline n", "Final n"]])

    st.markdown("**One-way sensitivity by prevalence**")
    sensitivity_rows = []
    for sp in [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]:
        out = detection_outputs(sp, se, c, m, rho, f)
        sensitivity_rows.append({"Prevalence": sp, "Final n": out["Final integrated n"]})
    sensitivity_df = pd.DataFrame(sensitivity_rows)
    st.dataframe(sensitivity_df, use_container_width=True, hide_index=True)
    st.line_chart(sensitivity_df.set_index("Prevalence"))

with tab5:
    st.subheader("User guide")
    st.markdown(
        """
        1. Enter expected prevalence, genomic detection sensitivity, target confidence, average cluster size, ICC/rho, and sequencing failure rate.
        2. Use **Final integrated n** for detection-based WGS or metagenomic AMR surveillance planning.
        3. Use the prevalence tab when the surveillance objective is estimating AMR prevalence with specified precision.
        4. Use the association tab when estimating sample size for genotype–phenotype association analysis.
        5. Use scenario analysis to compare conventional and integrated framework estimates.

        **Recommended parameter ranges preserved from the Excel calculator**

        - p: 0.001–0.5
        - Se: 0.5–1.0
        - C: 0.80–0.99
        - m: 1–50
        - ρ: 0–0.5
        - f: 0–0.5

        This tool supports planning and transparent assumptions. Parameter values should be justified using literature, pilot data, surveillance records, or expert elicitation.
        """
    )

st.divider()

current_detection = outputs_dataframe(results)
csv = current_detection.to_csv(index=False).encode("utf-8")
excel_bytes = make_excel_download(results, prev_results if "prev_results" in locals() else prevalence_outputs(0.1, 1.96, 0.05, 0.1), assoc_results if "assoc_results" in locals() else association_outputs(0.2, 2.0, 1.96, 0.84), scenarios if "scenarios" in locals() else pd.DataFrame())

col1, col2 = st.columns(2)
col1.download_button("Download current detection results as CSV", data=csv, file_name="wgs_amr_detection_results.csv", mime="text/csv")
col2.download_button("Download full results workbook", data=excel_bytes, file_name="wgs_amr_calculator_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

st.markdown("---")
st.subheader("Citation")

st.info("""
Debnath PP et al. WGS-AMR Surveillance Sample Size Calculator v1.0 (2026).

This application implements a systems-level framework for genomic surveillance sample size determination incorporating prevalence, genomic observability, clustering, and sequencing failure.
""")
