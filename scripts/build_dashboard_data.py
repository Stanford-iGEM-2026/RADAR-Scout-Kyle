"""Assemble per-disease dashboard data from the Modal-scored outputs.

Reads outputs/manifest.json + each <name>_racs.parquet / _umap.parquet / _meta.json
and writes compact JSON into dashboard/public/data/ for the multi-disease dashboard:
  - diseases.json                 : the disease/cohort switcher manifest
  - <name>_racs.json              : ranked gene table (top N, all metrics + DSS)
  - <name>_umap.json              : subsampled UMAP points
  - cross_cohort.json             : keloid cross-cohort consensus (if available)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("outputs")
DASH = Path("dashboard/public/data")
DASH.mkdir(parents=True, exist_ok=True)
TOP_N = 600


def add_dss(df):
    if "DSS" not in df.columns and {"log2FC", "mean_P"}.issubset(df.columns):
        df["DSS"] = df["log2FC"].clip(lower=0) * np.log10(df["mean_P"].clip(lower=0) + 1)
    return df


def emit_disease(name, disease, cell_type, cohort="CELLxGENE"):
    p = OUT / f"{name}_racs.parquet"
    if not p.exists():
        print(f"  [skip {name}] no racs parquet")
        return None
    df = add_dss(pd.read_parquet(p))
    df = df[df["RACS"].notna()].sort_values("RACS", ascending=False)
    fl = df.select_dtypes("float").columns
    df[fl] = df[fl].round(4)
    df.head(TOP_N).to_json(DASH / f"{name}_racs.json", orient="records")

    has_umap = (OUT / f"{name}_umap.parquet").exists()
    if has_umap:
        ud = pd.read_parquet(OUT / f"{name}_umap.parquet")
        f2 = ud.select_dtypes("float").columns
        ud[f2] = ud[f2].round(3)
        ud.to_json(DASH / f"{name}_umap.json", orient="records")

    meta = {}
    if (OUT / f"{name}_meta.json").exists():
        meta = json.load(open(OUT / f"{name}_meta.json"))
    return {"key": name, "disease": disease, "cell_type": cell_type, "cohort": cohort,
            "n_genes": int(len(df)), "n_donors": meta.get("n_donors", {}),
            "reference_pop": meta.get("reference_pop"), "has_umap": has_umap,
            "subpop": meta.get("subpop", {})}


def main():
    diseases = []
    seen = set()
    if (OUT / "manifest.json").exists():
        for d in json.load(open(OUT / "manifest.json"))["diseases"]:
            ent = emit_disease(d["out_name"], d["disease"], d.get("cell_type"))
            if ent:
                diseases.append(ent)
                seen.add(d["out_name"])
    # extra cohorts not in the census manifest (GEO + the Deng MFB sensor ranking)
    for name, dis, ct, coh in [
        ("keloid_geo", "keloid", "skin fibroblast", "GEO GSE163973"),
        ("keloid_deng_mfb", "keloid", "pathological MFB (Deng C3)", "Deng MFB vs normal scar"),
    ]:
        if name not in seen and (OUT / f"{name}_racs.parquet").exists():
            ent = emit_disease(name, dis, ct, cohort=coh)
            if ent:
                diseases.append(ent)

    json.dump({"diseases": diseases}, open(DASH / "diseases.json", "w"), indent=2)

    # cross-cohort consensus (from cross_cohort.py output)
    cc = OUT / "cross_cohort_DSS.csv"
    if cc.exists():
        c = pd.read_csv(cc).head(40).round(1)
        c.to_json(DASH / "cross_cohort.json", orient="records")

    print("diseases:", [f"{d['disease']}:{d['key']}(umap={d['has_umap']})" for d in diseases])


if __name__ == "__main__":
    main()
