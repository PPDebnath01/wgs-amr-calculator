# WGS-AMR Surveillance Sample Size Calculator

This Streamlit app converts the uploaded Excel calculator into an interactive web app for WGS/metagenomic AMR surveillance sample size planning.

## Files

- `app.py` — Streamlit app
- `requirements.txt` — Python dependencies

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Main calculator equation

```text
n_final = [ln(1-C) / ln(1-p×Se)] × [1+(m−1)ρ] × 1/(1−f)
```

## Preserved calculator components

- Detection calculator
- Prevalence estimation calculator
- Genotype–phenotype association calculator
- Scenario analysis
- One-way prevalence sensitivity analysis
- User guide and equation notes
- CSV and Excel result downloads

## Suggested deployment

You can deploy this app using Streamlit Community Cloud, Hugging Face Spaces, or an institutional server.
