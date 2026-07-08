"""Modal app: pull single-cell data from the CELLxGENE Census and score targets.

Heavy single-cell work runs here (per the project's Modal-first policy), not
locally. Two entrypoints:

  * ``census_disease_scan`` — CHEAP smoke test. Reads only the Census's
    precomputed summary-cell-counts table (tiny) to report which skin/scar/
    fibrosis diseases and how many cells are actually available. Run this first;
    it tells us whether keloid is in the Census or whether we fall back to a GEO
    loader / a related fibrosis for the first vertical.

  * ``build_and_score`` — pulls per-cell expression for the pathogenic (P),
    healthy (H), bystander (B), and related-disease (R) populations, normalizes
    to CP10k, and runs the donor-aware RACS scorer, writing a ranked table to a
    Modal Volume.

Usage:
    modal run modal_app/census_pull.py                 # smoke test (scan)
    modal run modal_app/census_pull.py::build_and_score --disease keloid ...
"""

from __future__ import annotations

import json

import modal

app = modal.App("radar-scout")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        # must be recent enough to read the current Census SOMA encoding (>=1.1.0);
        # pulls a compatible tiledbsoma. Older pins fail with "Unsupported SOMA
        # object encoding version".
        "cellxgene-census==1.18.0",
        "scanpy",
        "pyarrow",
    )
    # ship the scoring package so RACS runs on Modal
    .add_local_dir("radar_scout", remote_path="/root/radar_scout")
)

vol = modal.Volume.from_name("radar-scout-data", create_if_missing=True)
DATA = "/data"

CENSUS_VERSION = "2025-11-08"  # pinned LTS for reproducibility (was "stable" on 2026-07-08)

# Terms we care about for the keloid / fibrosis first vertical. Extended via CLI.
SCAN_PATTERNS = ["keloid", "scar", "fibro", "sclerosis", "dermat", "skin", "wound", "cheloid"]


@app.function(image=image, timeout=900)
def census_disease_scan(patterns: list[str] | None = None) -> dict:
    """Report available diseases (cell counts) matching patterns. Cheap.

    Schema-robust: the summary_cell_counts column names vary across Census
    versions, so we detect the organism/category/label/count columns dynamically
    and echo the real column list.
    """
    import cellxgene_census

    patterns = [p.lower() for p in (patterns or SCAN_PATTERNS)]
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        sccc = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()

    print("SCCC COLUMNS:", list(sccc.columns), "| SHAPE", sccc.shape)
    lc = {c.lower(): c for c in sccc.columns}

    def col(*cands):
        for c in cands:
            if c in lc:
                return lc[c]
        return None

    org_c = col("organism")
    cat_c = col("category")
    lab_c = col("label", "ontology_term_label", "ontology_term_id")
    cnt_c = col("total_cell_count", "unique_cell_count", "n_cells", "count")

    # probe the actual vocabulary (varies across Census versions)
    print("CATEGORIES:", sorted(sccc[cat_c].astype(str).unique())[:30] if cat_c else None)
    print("ORGANISMS:", sorted(sccc[org_c].astype(str).unique())[:10] if org_c else None)

    diseases = sccc[sccc[cat_c].astype(str).str.lower() == "disease"].copy() if cat_c else sccc.copy()
    if org_c is not None and diseases.shape[0]:
        org_vals = diseases[org_c].astype(str)
        mask = org_vals.str.contains("sapiens", case=False, na=False) | org_vals.str.contains("9606", na=False)
        if mask.any():
            diseases = diseases[mask].copy()
    diseases["_label_l"] = diseases[lab_c].astype(str).str.lower()
    hit = diseases[diseases["_label_l"].apply(lambda s: any(p in s for p in patterns))]
    have_cnt = cnt_c is not None and cnt_c in hit.columns
    if have_cnt:
        hit = hit.sort_values(cnt_c, ascending=False)

    out = {
        "census_version": CENSUS_VERSION,
        "columns": list(sccc.columns),
        "count_col": cnt_c,
        "n_disease_terms": int(diseases.shape[0]),
        "matched_diseases": [
            {"disease": r[lab_c], "cells": int(r[cnt_c]) if have_cnt else None}
            for _, r in hit.iterrows()
        ],
    }
    print(json.dumps(out, indent=2))
    return out


def _fetch_pop(census, obs_filter, label, tissue):
    """Fetch one population as AnnData, tagged with population label."""
    import cellxgene_census

    adata = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        obs_value_filter=obs_filter,
        obs_column_names=["donor_id", "cell_type", "disease", "tissue_general", "assay"],
    )
    adata.obs["radar_pop"] = label
    return adata


@app.function(image=image, volumes={DATA: vol}, timeout=3600, memory=65536, cpu=8.0)
def build_and_score(
    disease: str = "keloid",
    pathogenic_cell_types: list[str] | None = None,
    related_diseases: list[str] | None = None,
    tissue: str = "skin",
    min_detect_frac_in_P: float = 0.10,
    out_name: str = "keloid_fibroblast",
) -> dict:
    """Assemble P/H/B/R populations, normalize, and run RACS. Writes to the Volume."""
    import cellxgene_census
    import numpy as np
    import scanpy as sc

    from radar_scout.scoring import score_matrix

    pathogenic_cell_types = pathogenic_cell_types or ["fibroblast"]
    related_diseases = related_diseases or []
    ct_in = "[" + ", ".join(f"'{c}'" for c in pathogenic_cell_types) + "]"

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        parts = []
        # P: pathogenic cell type in the target disease
        parts.append(_fetch_pop(
            census, f"disease == '{disease}' and cell_type in {ct_in} and tissue_general == '{tissue}'",
            "P", tissue))
        # H: same cell type, healthy tissue
        parts.append(_fetch_pop(
            census, f"disease == 'normal' and cell_type in {ct_in} and tissue_general == '{tissue}'",
            "H", tissue))
        # B: bystander cell types in the diseased tissue
        parts.append(_fetch_pop(
            census, f"disease == '{disease}' and cell_type not in {ct_in} and tissue_general == '{tissue}'",
            "B", tissue))
        # R: pathogenic cell type in related diseases
        if related_diseases:
            rd_in = "[" + ", ".join(f"'{d}'" for d in related_diseases) + "]"
            parts.append(_fetch_pop(
                census, f"disease in {rd_in} and cell_type in {ct_in} and tissue_general == '{tissue}'",
                "R", tissue))

    import anndata as ad
    adata = ad.concat([p for p in parts if p.n_obs > 0], join="inner", merge="same")
    if adata.n_obs == 0:
        raise RuntimeError(f"No cells found for disease='{disease}' — run census_disease_scan first.")

    # normalize to CP10k (linear) — the units the Hill model expects
    sc.pp.normalize_total(adata, target_sum=1e4)
    expr = adata.X
    expr = np.asarray(expr.todense()) if hasattr(expr, "todense") else np.asarray(expr)

    pop = adata.obs["radar_pop"].to_numpy()
    donor = adata.obs["donor_id"].astype(str).to_numpy()

    # candidate genes: detected in >= min_detect_frac of pathogenic cells
    is_p = pop == "P"
    detect = (expr[is_p] > 0).mean(axis=0)
    keep = np.where(detect >= min_detect_frac_in_P)[0]
    genes = adata.var["feature_name"].to_numpy() if "feature_name" in adata.var else adata.var_names.to_numpy()

    df = score_matrix(expr[:, keep], genes[keep], donor, pop, pos_label="P")

    counts = {p: int((pop == p).sum()) for p in ["P", "H", "B", "R"]}
    n_donor = {p: int(len(set(donor[pop == p]))) for p in ["P", "H", "B", "R"]}

    out_path = f"{DATA}/{out_name}_racs.parquet"
    df.to_parquet(out_path)
    meta = {"disease": disease, "tissue": tissue, "cell_counts": counts,
            "n_donors": n_donor, "n_candidate_genes": int(len(keep)),
            "census_version": CENSUS_VERSION, "out": out_path,
            "top20": df.head(20).to_dict(orient="records")}
    with open(f"{DATA}/{out_name}_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=float)
    vol.commit()
    print(json.dumps({k: meta[k] for k in ("cell_counts", "n_donors", "n_candidate_genes")}, indent=2))
    return meta


@app.local_entrypoint()
def main():
    """Default: run the cheap disease scan and print what's available."""
    res = census_disease_scan.remote()
    print("\n=== RADAR-Scout Census scan ===")
    print("summary_cell_counts columns:", res.get("columns"))
    for d in res["matched_diseases"]:
        cells = d.get("cells")
        print(f"  {d['disease']:45s} {cells:>12,} cells" if cells is not None else f"  {d['disease']}")
    if not res["matched_diseases"]:
        print("  (no matching diseases in Census — fall back to a GEO loader)")
