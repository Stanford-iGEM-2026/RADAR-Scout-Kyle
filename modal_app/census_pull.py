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


@app.function(image=image, timeout=900)
def population_probe(disease: str = "keloid", cell_types: str = "fibroblast",
                     tissue: str = "skin", related: str = "") -> dict:
    """Metadata-only donor/cell/dataset counts for P and H (and related diseases).

    Cheap (reads obs columns only, no expression). Answers the critical question:
    does this disease have enough biological donors for donor-aware statistics?

    ``cell_types`` and ``related`` are comma-separated (Modal CLI can't parse list
    annotations). e.g. --cell-types "fibroblast,myofibroblast" --related "cystic fibrosis".
    """
    import cellxgene_census

    cts = [c.strip() for c in cell_types.split(",") if c.strip()]
    related_list = [r.strip() for r in related.split(",") if r.strip()]
    ct_in = "[" + ", ".join(f"'{c}'" for c in cts) + "]"
    cols = ["donor_id", "cell_type", "disease", "tissue_general", "assay", "dataset_id"]

    def summarize(df, name):
        return {
            "population": name, "n_cells": int(len(df)),
            "n_donors": int(df["donor_id"].nunique()) if len(df) else 0,
            "n_datasets": int(df["dataset_id"].nunique()) if len(df) else 0,
            "assays": sorted(df["assay"].astype(str).unique())[:6] if len(df) else [],
            "cells_per_donor": (df.groupby("donor_id").size().describe()[["min", "50%", "max"]]
                                .astype(int).to_dict() if len(df) else {}),
        }

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        obs = census["census_data"]["homo_sapiens"].obs

        def q(vf):
            return obs.read(value_filter=vf, column_names=cols).concat().to_pandas()

        out = {"disease": disease, "cell_types": cts, "tissue": tissue, "populations": []}
        out["populations"].append(summarize(
            q(f"disease == '{disease}' and cell_type in {ct_in} and tissue_general == '{tissue}'"), "P (pathogenic)"))
        out["populations"].append(summarize(
            q(f"disease == 'normal' and cell_type in {ct_in} and tissue_general == '{tissue}'"), "H (healthy)"))
        for rd in related_list:
            out["populations"].append(summarize(
                q(f"disease == '{rd}' and cell_type in {ct_in}"), f"R: {rd}"))

    print(json.dumps(out, indent=2))
    return out


@app.function(image=image, timeout=900)
def disease_composition(disease: str = "keloid", tissue: str = "") -> dict:
    """Dump the ACTUAL cell_type / tissue labels used for a disease (no expression).

    Reveals the real annotation vocabulary so we can build correct filters (the
    generic 'fibroblast' CL term often isn't what datasets use)."""
    import cellxgene_census

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        obs = census["census_data"]["homo_sapiens"].obs
        vf = f"disease == '{disease}'"
        if tissue:
            vf += f" and tissue_general == '{tissue}'"
        df = obs.read(
            value_filter=vf,
            column_names=["donor_id", "cell_type", "tissue_general", "tissue", "assay", "dataset_id"],
        ).concat().to_pandas()

    out = {
        "disease": disease,
        "n_cells": int(len(df)),
        "n_donors": int(df["donor_id"].nunique()) if len(df) else 0,
        "n_datasets": int(df["dataset_id"].nunique()) if len(df) else 0,
        "tissue_general": {k: int(v) for k, v in df["tissue_general"].value_counts().head(10).items()},
        "cell_types": {k: int(v) for k, v in df["cell_type"].value_counts().head(30).items()},
    }
    print(json.dumps(out, indent=2))
    return out


@app.function(image=image, timeout=900)
def celltype_across_disease(cell_type: str = "skin fibroblast", tissue: str = "") -> dict:
    """Per-disease cell/donor/dataset counts for a given cell type. Gives P, H, and
    the candidate R (related-disease) set in one shot."""
    import cellxgene_census

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        obs = census["census_data"]["homo_sapiens"].obs
        vf = f"cell_type == '{cell_type}'"
        if tissue:
            vf += f" and tissue_general == '{tissue}'"
        df = obs.read(
            value_filter=vf,
            column_names=["donor_id", "disease", "tissue_general", "dataset_id"],
        ).concat().to_pandas()

    g = (df.groupby("disease")
         .agg(n_cells=("donor_id", "size"), n_donors=("donor_id", "nunique"),
              n_datasets=("dataset_id", "nunique"))
         .sort_values("n_cells", ascending=False))
    g = g[g["n_cells"] > 0]
    out = {"cell_type": cell_type, "total_cells": int(len(df)),
           "by_disease": [{"disease": d, "n_cells": int(r.n_cells),
                           "n_donors": int(r.n_donors), "n_datasets": int(r.n_datasets)}
                          for d, r in g.head(40).iterrows()]}
    print(json.dumps(out, indent=2))
    return out


def _fetch_pop(census, obs_filter, label, max_cells, rng):
    """Fetch one population as AnnData, tag it, subsample to max_cells.

    Robust: returns None on filter error (e.g. unsupported 'not in') or empty result.
    """
    import cellxgene_census
    import numpy as np

    try:
        adata = cellxgene_census.get_anndata(
            census, organism="Homo sapiens", obs_value_filter=obs_filter,
            obs_column_names=["donor_id", "cell_type", "disease", "tissue_general", "assay"],
        )
    except Exception as e:
        print(f"[skip {label}] {type(e).__name__}: {e}")
        return None
    if adata.n_obs == 0:
        print(f"[empty {label}]")
        return None
    if adata.n_obs > max_cells:
        idx = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
        adata = adata[idx].copy()
    adata.obs["radar_pop"] = label
    print(f"[{label}] {adata.n_obs} cells, {adata.obs['donor_id'].nunique()} donors")
    return adata


@app.function(image=image, volumes={DATA: vol}, timeout=3600, memory=49152, cpu=8.0)
def build_and_score(
    disease: str = "keloid",
    pathogenic_cell_types: str = "skin fibroblast",
    related_diseases: str = "",
    tissue: str = "skin of body",
    min_detect_frac: float = 0.10,
    max_cells_per_pop: int = 40000,
    exclude_technical: bool = True,
    out_name: str = "keloid_fibroblast",
) -> dict:
    """Assemble P/H/B/R populations, normalize, and run RACS. Writes to the Volume.

    ``pathogenic_cell_types`` and ``related_diseases`` are comma-separated strings.
    Related diseases are collapsed into one off-target group 'R' for scoring.
    """
    import cellxgene_census
    import numpy as np
    import scanpy as sc
    import anndata as ad

    from radar_scout.scoring import score_matrix
    from radar_scout.genesets import filter_technical, reasons

    pct = [c.strip() for c in pathogenic_cell_types.split(",") if c.strip()]
    related = [r.strip() for r in related_diseases.split(",") if r.strip()]
    ct_in = "[" + ", ".join(f"'{c}'" for c in pct) + "]"
    tis = f" and tissue_general == '{tissue}'" if tissue else ""
    rng = np.random.default_rng(0)

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        parts = []
        p = _fetch_pop(census, f"disease == '{disease}' and cell_type in {ct_in}{tis}",
                       "P", max_cells_per_pop, rng)
        if p is None:
            raise RuntimeError(f"No pathogenic cells: disease='{disease}', cell_type in {pct}")
        parts.append(p)
        parts.append(_fetch_pop(census, f"disease == 'normal' and cell_type in {ct_in}{tis}",
                                "H", max_cells_per_pop, rng))
        parts.append(_fetch_pop(census, f"disease == '{disease}' and cell_type not in {ct_in}{tis}",
                                "B", max_cells_per_pop, rng))
        for rd in related:  # each related disease -> its own tag, collapsed to R below
            parts.append(_fetch_pop(census, f"disease == '{rd}' and cell_type in {ct_in}",
                                    f"R::{rd}", max_cells_per_pop, rng))

    parts = [x for x in parts if x is not None]
    adata = ad.concat(parts, join="inner", merge="same")

    # normalize to CP10k (linear) — the units the Hill model expects
    sc.pp.normalize_total(adata, target_sum=1e4)

    pop_detail = adata.obs["radar_pop"].astype(str).to_numpy()
    pop = np.array(["R" if s.startswith("R::") else s for s in pop_detail])  # collapse related
    donor = adata.obs["donor_id"].astype(str).to_numpy()

    # candidate genes: detected in >= min_detect_frac of pathogenic cells (sparse-safe)
    Xp = adata[pop == "P"].X
    detect = np.asarray((Xp > 0).mean(axis=0)).ravel()
    keep = np.where(detect >= min_detect_frac)[0]

    sub = adata[:, keep]
    X = sub.X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)  # only candidate genes
    genes = (sub.var["feature_name"] if "feature_name" in sub.var else sub.var_names).to_numpy().astype(str)

    n_before = len(genes)
    if exclude_technical:
        mask = filter_technical(genes)
        removed = reasons(genes[~mask])
        print(f"[filter] dropped {int((~mask).sum())}/{n_before} technical genes "
              f"(sex/ribosomal/mito/ncRNA/IEG); e.g. {list(removed)[:8]}")
        X = X[:, mask]
        genes = genes[mask]

    df = score_matrix(X, genes, donor, pop, pos_label="P")

    ucounts = sorted(set(pop_detail))
    counts = {u: int((pop_detail == u).sum()) for u in ucounts}
    ndon = {u: int(len(set(donor[pop_detail == u]))) for u in ucounts}

    out_path = f"{DATA}/{out_name}_racs.parquet"
    df.to_parquet(out_path)
    meta = {"disease": disease, "tissue": tissue, "pathogenic_cell_types": pct,
            "related_diseases": related, "cell_counts": counts, "n_donors": ndon,
            "n_candidate_genes": int(len(genes)), "n_candidate_prefilter": int(n_before),
            "exclude_technical": bool(exclude_technical), "census_version": CENSUS_VERSION,
            "out": out_path, "top25": df.head(25).to_dict(orient="records")}
    with open(f"{DATA}/{out_name}_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=float)
    vol.commit()

    print(json.dumps({"cell_counts": counts, "n_donors": ndon,
                      "n_candidate_genes": int(len(keep))}, indent=2))
    print("TOP 15 RADAR targets:")
    for r in df.head(15).to_dict("records"):
        print(f"  {str(r['gene'])[:16]:16s} RACS={r['RACS']:.3f}  Sep={r['Sep']:.3f}  "
              f"Feas={r['Feas']:.3f}  Repro={r.get('Repro', float('nan')):.3f}  Off={r['OffMax']:.3f}")
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
