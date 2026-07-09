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
        "leidenalg",
        "igraph",
        "harmonypy",   # batch correction
        "scikit-image",  # scrublet dep
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


def _save_paga(lg, pop, path):
    """Save the PAGA graph (cluster connectivities + UMAP centroids + population mix)."""
    import json

    import numpy as np
    import pandas as pd

    conn = lg.uns["paga"]["connectivities"]
    conn = conn.toarray() if hasattr(conn, "toarray") else np.asarray(conn)
    clusters = lg.obs["leiden"].to_numpy()
    cats = list(lg.obs["leiden"].cat.categories)
    xy = lg.obsm["X_umap"]
    nodes = [{"cluster": str(c), "x": float(xy[clusters == c, 0].mean()),
              "y": float(xy[clusters == c, 1].mean()), "n": int((clusters == c).sum()),
              "pop": {k: round(float(v), 2) for k, v in
                      pd.Series(pop[clusters == c]).value_counts(normalize=True).items()}}
             for c in cats]
    edges = [{"a": str(cats[i]), "b": str(cats[j]), "w": round(float(conn[i, j]), 3)}
             for i in range(len(cats)) for j in range(i + 1, len(cats)) if conn[i, j] > 0.05]
    json.dump({"nodes": nodes, "edges": edges}, open(path, "w"))


@app.function(image=image, volumes={DATA: vol}, timeout=7200, memory=65536, cpu=8.0)
def build_and_score(
    disease: str = "melanoma",
    pathogenic_cell_types: str = "malignant cell",
    related_diseases: str = "",
    tissue: str = "",
    min_detect_frac: float = 0.10,
    max_cells_per_pop: int = 40000,
    exclude_technical: bool = True,
    compute_umap: bool = True,
    subcluster: bool = False,
    subcluster_resolution: float = 1.0,
    out_name: str = "",
) -> dict:
    """Disease-agnostic RADAR target prioritization on the CELLxGENE Census.

    Pulls the whole diseased tissue and splits it into pathogenic P (cell types
    matching ``pathogenic_cell_types``) vs bystander B (all other cell types ->
    cell-type specificity), plus healthy counterpart H (same cell types in 'normal')
    and related-disease R. Computes the full RACS metric table, a UMAP, and
    donor-level DE; writes everything to the Volume and updates a dashboard manifest.
    """
    import cellxgene_census
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import anndata as ad

    from radar_scout.scoring import score_matrix
    from radar_scout.genesets import filter_technical
    from radar_scout.de import pseudobulk_de

    pct = [c.strip().lower() for c in pathogenic_cell_types.split(",") if c.strip()]
    related = [r.strip() for r in related_diseases.split(",") if r.strip()]
    out_name = out_name or disease.replace(" ", "_")
    tis = f" and tissue_general == '{tissue}'" if tissue else ""
    rng = np.random.default_rng(0)

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        # resolve the target cell-type labels actually present for this disease (cheap obs read),
        # so the pathogenic population is fetched directly and is well-sampled even when it is
        # a minority of the tissue (e.g. fibroblasts in lung).
        obs = census["census_data"]["homo_sapiens"].obs
        avail = obs.read(value_filter=f"disease == '{disease}'{tis}",
                         column_names=["cell_type"]).concat().to_pandas()["cell_type"].astype(str)
        uniq = sorted(avail.unique())
        # exact match, else the requested term is contained in a more specific label
        # (e.g. "fibroblast" -> "alveolar fibroblast"). Only the query-in-candidate
        # direction — the reverse wrongly matched "T cell" inside "malignan-t cell".
        resolved = [c for c in uniq if c.lower() in pct or any(t in c.lower() for t in pct)]
        if not resolved:
            raise RuntimeError(f"No cell types matched {pct}; available: {uniq[:25]}")
        lab_in = "[" + ", ".join(f"'{c}'" for c in resolved) + "]"
        p_labels = resolved
        print(f"[resolve] pathogenic cell types: {resolved}")

        # P: dedicated fetch of the target cell type in the disease (well-sampled)
        p = _fetch_pop(census, f"disease == '{disease}' and cell_type in {lab_in}{tis}",
                       "P", max_cells_per_pop, rng)
        if p is None:
            raise RuntimeError(f"No pathogenic cells for {disease} / {resolved}")
        parts = [p]
        # B: bystander cell types in the diseased tissue (cell-type specificity + UMAP context)
        b = _fetch_pop(census, f"disease == '{disease}'{tis}", "B", max_cells_per_pop, rng)
        if b is not None:
            b = b[~b.obs["cell_type"].astype(str).isin(resolved)].copy()
            if b.n_obs > 0:
                parts.append(b)
        # H: healthy counterpart (same cell types, disease == normal) — may be empty (cancer)
        parts.append(_fetch_pop(census, f"disease == 'normal' and cell_type in {lab_in}{tis}",
                                "H", max_cells_per_pop, rng))
        # R: related diseases (same cell types)
        for rd in related:
            parts.append(_fetch_pop(census, f"disease == '{rd}' and cell_type in {lab_in}",
                                    f"R::{rd}", max_cells_per_pop, rng))

    parts = [x for x in parts if x is not None and x.n_obs > 0]
    adata = ad.concat(parts, join="inner", merge="same")
    adata.obs_names_make_unique()

    # QC: doublet removal (per donor) on raw counts (spec: data harmonization)
    try:
        sc.pp.scrublet(adata, batch_key="donor_id")
        db = adata.obs.get("predicted_doublet")
        if db is not None:
            n_db = int(db.sum())
            adata = adata[~db.to_numpy()].copy()
            print(f"[qc] removed {n_db} predicted doublets")
    except Exception as e:
        print(f"[qc] scrublet skipped: {type(e).__name__}: {e}")

    adata.layers["counts"] = adata.X.copy()         # raw counts (carried for DE)
    sc.pp.normalize_total(adata, target_sum=1e4)    # linear CP10k for scoring

    pop_detail = adata.obs["radar_pop"].astype(str).to_numpy()
    pop = np.array(["R" if s.startswith("R::") else s for s in pop_detail])
    donor = adata.obs["donor_id"].astype(str).to_numpy()
    genes_all = (adata.var["feature_name"] if "feature_name" in adata.var
                 else adata.var_names).to_numpy().astype(str)

    # spec Task 4: identify the disease-enriched (pathogenic) subpopulation within the
    # target cell type and restrict P to it. This surfaces cell-STATE markers that are
    # diluted across the whole cell type (e.g. keloid mesenchymal fibroblasts:
    # POSTN/ASPN/ADAM12), mirroring how melanoma's malignant cells are already a state.
    subpop_info = {}
    if subcluster:
        mask_pt = np.isin(pop, ["P", "H"])
        subad = adata[mask_pt].copy()
        sc.pp.log1p(subad)
        sc.pp.highly_variable_genes(subad, n_top_genes=2000)
        sc.pp.pca(subad, n_comps=30)
        sc.pp.neighbors(subad, n_neighbors=15)
        try:
            sc.tl.leiden(subad, resolution=subcluster_resolution, flavor="igraph",
                         n_iterations=2, directed=False)
        except (TypeError, ImportError):
            sc.tl.leiden(subad, resolution=subcluster_resolution)
        clusters = subad.obs["leiden"].to_numpy()
        is_dis = (pop[mask_pt] == "P").astype(float)
        baseline = float(is_dis.mean())
        grp = pd.Series(is_dis).groupby(clusters)
        dfrac, csize = grp.mean(), grp.size()
        # a genuine pathogenic state must be enriched above the disease baseline AND sizeable
        elig = dfrac[(csize >= 30) & (dfrac > max(1.5 * baseline, baseline + 0.05))]
        idx_pt = np.where(mask_pt)[0]
        if len(elig):
            path_cluster = elig.idxmax()
            n_after = int(((pop[idx_pt] == "P") & (clusters == path_cluster)).sum())
            if n_after >= 20:
                reassign = idx_pt[(pop[idx_pt] == "P") & (clusters != path_cluster)]
                pop = pop.copy()
                pop[reassign] = "B"  # disease cells in non-pathogenic states -> bystander
                pop_detail = pop.copy()
                subpop_info = {"path_cluster": str(path_cluster), "disease_frac": float(dfrac[path_cluster]),
                               "baseline_frac": baseline, "n_clusters": int(csize.size), "n_P_after": n_after}
                print(f"[subpop] cluster {path_cluster}: disease_frac={dfrac[path_cluster]:.2f} "
                      f"(baseline {baseline:.2f}), P={n_after} state cells")
            else:
                print(f"[subpop] enriched cluster too small (n={n_after}); keeping whole-cell-type P")
        else:
            print(f"[subpop] no disease-enriched cluster (baseline {baseline:.2f}); keeping whole-cell-type P")

    # candidate genes: detected in >= min_detect_frac of pathogenic cells
    Xp = adata[pop == "P"].X
    detect = np.asarray((Xp > 0).mean(axis=0)).ravel()
    keep = np.where(detect >= min_detect_frac)[0]
    n_before = len(keep)
    if exclude_technical:
        m = filter_technical(genes_all[keep])
        print(f"[filter] dropped {int((~m).sum())}/{n_before} technical genes")
        keep = keep[m]

    # a single candidate-gene subset drives BOTH scoring and DE (aligned shapes)
    sub = adata[:, keep]
    genes = (sub.var["feature_name"] if "feature_name" in sub.var else sub.var_names).to_numpy().astype(str)
    Xc = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
    df = score_matrix(Xc, genes, donor, pop, pos_label="P")

    # donor-level DE: P vs primary reference (healthy if present, else bystander)
    ref = "H" if (pop == "H").any() else "B"
    de = None
    try:
        rc = sub.layers["counts"]
        rc = np.asarray(rc.todense()) if hasattr(rc, "todense") else np.asarray(rc)
        de = pseudobulk_de(rc, donor, pop, genes, pos_label="P", neg_labels=[ref])
        de.to_parquet(f"{DATA}/{out_name}_de.parquet")
    except Exception as e:
        print(f"[de] skipped: {type(e).__name__}: {e}")

    # UMAP for the dashboard (subsampled to <=8k cells), batch-corrected + PAGA
    umap_saved = False
    if compute_umap:
        try:
            lg = adata.copy()
            sc.pp.log1p(lg)
            sc.pp.highly_variable_genes(lg, n_top_genes=2000)
            sc.pp.pca(lg, n_comps=30)
            rep = "X_pca"
            try:  # batch correction across donors (spec: data harmonization)
                sc.external.pp.harmony_integrate(lg, "donor_id")
                rep = "X_pca_harmony"
            except Exception as e:
                print(f"[umap] harmony skipped: {type(e).__name__}: {e}")
            sc.pp.neighbors(lg, n_neighbors=15, use_rep=rep)
            sc.tl.umap(lg)
            try:  # PAGA trajectory over leiden clusters
                try:
                    sc.tl.leiden(lg, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
                except (TypeError, ImportError):
                    sc.tl.leiden(lg, resolution=1.0)
                sc.tl.paga(lg, groups="leiden")
                _save_paga(lg, pop, f"{DATA}/{out_name}_paga.json")
            except Exception as e:
                print(f"[paga] skipped: {type(e).__name__}: {e}")
            xy = lg.obsm["X_umap"]
            idx = np.sort(rng.choice(lg.n_obs, size=min(lg.n_obs, 8000), replace=False))
            gpos = {g: i for i, g in enumerate(genes)}
            ud = pd.DataFrame({"UMAP1": xy[idx, 0], "UMAP2": xy[idx, 1],
                               "cell_type": lg.obs["cell_type"].astype(str).to_numpy()[idx],
                               "pop": pop[idx], "donor": donor[idx]})
            for g in df.head(12)["gene"].tolist():
                if g in gpos:
                    ud[g] = Xc[idx, gpos[g]]
            ud.to_parquet(f"{DATA}/{out_name}_umap.parquet")
            umap_saved = True
        except Exception as e:
            print(f"[umap] skipped: {type(e).__name__}: {e}")

    ucounts = sorted(set(pop_detail))
    counts = {u: int((pop_detail == u).sum()) for u in ucounts}
    ndon = {u: int(len(set(donor[pop_detail == u]))) for u in ucounts}

    df.to_parquet(f"{DATA}/{out_name}_racs.parquet")
    meta = {"disease": disease, "tissue": tissue, "pathogenic_cell_types": pct,
            "pathogenic_labels": p_labels, "related_diseases": related,
            "cell_counts": counts, "n_donors": ndon, "n_candidate_genes": int(len(genes)),
            "n_candidate_prefilter": int(n_before), "reference_pop": ref,
            "has_umap": umap_saved, "has_de": de is not None, "subpop": subpop_info,
            "census_version": CENSUS_VERSION, "top25": df.head(25).to_dict(orient="records")}
    with open(f"{DATA}/{out_name}_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=float)

    # update the dashboard manifest of scored diseases
    man_path = f"{DATA}/manifest.json"
    try:
        manifest = json.load(open(man_path))
    except Exception:
        manifest = {"diseases": []}
    manifest["diseases"] = [d for d in manifest["diseases"] if d.get("out_name") != out_name]
    manifest["diseases"].append({"disease": disease, "out_name": out_name,
                                 "cell_type": pathogenic_cell_types, "n_genes": int(len(genes)),
                                 "n_donors_P": ndon.get("P", 0), "reference_pop": ref,
                                 "has_umap": umap_saved})
    json.dump(manifest, open(man_path, "w"), indent=2)
    vol.commit()

    print(json.dumps({"disease": disease, "cell_counts": counts, "n_donors": ndon,
                      "n_candidate_genes": int(len(genes)), "reference": ref,
                      "umap": umap_saved, "de": de is not None}, indent=2))
    print("TOP 15 RADAR targets:")
    for r in df.head(15).to_dict("records"):
        print(f"  {str(r['gene'])[:16]:16s} RACS={r['RACS']:.3f}  Sep={r['Sep']:.3f}  "
              f"Feas={r['Feas']:.3f}  Off={r['OffMax']:.3f}  log2FC={r.get('log2FC', float('nan')):.2f}")
    return meta


@app.function(image=image, volumes={DATA: vol}, timeout=5400, memory=49152, cpu=8.0)
def ingest_and_score_geo(gse: str = "GSE163973", disease_label: str = "keloid",
                         out_name: str = "keloid_geo", min_detect_frac: float = 0.03) -> dict:
    """Ingest an INDEPENDENT GEO scRNA cohort (10x RAW.tar), gate fibroblasts by
    markers, and run RACS. This is the second keloid cohort (disease vs normal-scar)
    for cross-cohort validation. GSE163973 = Deng 2021 (3 keloid + 3 normal scar).
    """
    import gzip
    import os
    import re
    import shutil
    import tarfile
    import urllib.request

    import anndata as ad
    import numpy as np
    import scanpy as sc
    from scipy.io import mmread

    from radar_scout.scoring import score_matrix
    from radar_scout.genesets import filter_technical

    wd = "/tmp/geo"
    os.makedirs(wd, exist_ok=True)
    stub = gse[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/{gse}_RAW.tar"
    tarp = f"{wd}/{gse}_RAW.tar"
    if not os.path.exists(tarp):
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, tarp)
    with tarfile.open(tarp) as t:
        t.extractall(wd)

    # discover all files recursively, group by GSM id
    allfiles = []
    for root, _, fs in os.walk(wd):
        for f in fs:
            if f.startswith("GSM"):
                allfiles.append(os.path.join(root, f))
    groups: dict = {}
    for f in allfiles:
        m = re.match(r"(GSM\d+)", os.path.basename(f))
        if m:
            groups.setdefault(m.group(1), []).append(f)

    def _read_mtx(path):
        if path.endswith(".gz"):
            tmp = path[:-3]
            with gzip.open(path) as fin, open(tmp, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            path = tmp
        return mmread(path).tocsr()

    def _lines(path):
        op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
        return [ln.rstrip("\n") for ln in op]

    adatas = []
    for gsm, fs in sorted(groups.items()):
        search = list(fs)
        # some series nest each sample as GSM..._matrix.tar.gz (a 10x bundle inside the tar)
        nested = next((f for f in fs if f.endswith(".tar.gz") and ".mtx" not in f.lower()), None)
        if nested and not any(".mtx" in f.lower() for f in fs):
            gdir = os.path.join(wd, gsm + "_x")
            os.makedirs(gdir, exist_ok=True)
            with tarfile.open(nested) as t:
                t.extractall(gdir)
            search = [os.path.join(r, f) for r, _, ffs in os.walk(gdir) for f in ffs]
        mtxf = next((f for f in search if "matrix" in os.path.basename(f).lower() and ".mtx" in f.lower()), None)
        bcf = next((f for f in search if "barcode" in os.path.basename(f).lower()), None)
        ftf = next((f for f in search if "feature" in os.path.basename(f).lower() or "genes" in os.path.basename(f).lower()), None)
        if not (mtxf and bcf and ftf):
            print(f"[skip {gsm}] missing 10x files among: {[os.path.basename(x) for x in search][:6]}")
            continue
        X = _read_mtx(mtxf)  # genes x cells
        barcodes = [ln.split("\t")[0] for ln in _lines(bcf)]
        symbols = [(ln.split("\t")[1] if "\t" in ln else ln) for ln in _lines(ftf)]
        if X.shape[0] == len(symbols) and X.shape[1] == len(barcodes):
            X = X.T.tocsr()  # -> cells x genes
        a = ad.AnnData(X)
        a.var_names = symbols[: a.n_vars]
        a.var_names_make_unique()
        a.obs_names = [f"{gsm}_{b}" for b in barcodes[: a.n_obs]]
        a.obs["donor"] = gsm
        blob = " ".join(os.path.basename(x) for x in fs).upper()
        a.obs["condition"] = "keloid" if ("KL" in blob or "KELOID" in blob) else "normal_scar"
        print(f"[{gsm}] {a.n_obs} cells x {a.n_vars} genes -> {a.obs['condition'][0]}")
        adatas.append(a)
    if not adatas:
        raise RuntimeError("no 10x samples parsed from RAW.tar")

    adata = ad.concat(adatas, join="inner")
    adata.obs_names_make_unique()

    # QC
    sc.pp.filter_cells(adata, min_genes=200)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)  # CP10k
    lg = adata.copy()
    sc.pp.log1p(lg)

    def mscore(genes):
        gg = [g for g in genes if g in lg.var_names]
        return np.asarray(lg[:, gg].X.mean(axis=1)).ravel() if gg else np.zeros(lg.n_obs)

    fib = mscore(["COL1A1", "COL1A2", "LUM", "DCN", "PDGFRA", "COL3A1"])
    imm = mscore(["PTPRC", "CD3D", "CD68", "LYZ"])
    endo = mscore(["PECAM1", "VWF", "CLDN5"])
    kera = mscore(["KRT14", "KRT5", "KRT1"])
    is_fib = (fib > 0.5) & (fib > imm) & (fib > endo) & (fib > kera)
    print(f"fibroblasts gated: {int(is_fib.sum())}/{adata.n_obs}")
    adata = adata[is_fib].copy()

    pop = np.where(adata.obs["condition"].values == disease_label, "P", "H")
    donor = adata.obs["donor"].astype(str).to_numpy()
    genes_all = adata.var_names.to_numpy().astype(str)
    Xp = adata[pop == "P"].X
    detect = np.asarray((Xp > 0).mean(axis=0)).ravel()
    keep = np.where(detect >= min_detect_frac)[0]
    genes = genes_all[keep]
    m = filter_technical(genes)
    keep, genes = keep[m], genes[m]
    Xc = adata[:, keep].X
    Xc = np.asarray(Xc.todense()) if hasattr(Xc, "todense") else np.asarray(Xc)

    df = score_matrix(Xc, genes, donor, pop, pos_label="P")
    df.to_parquet(f"{DATA}/{out_name}_racs.parquet")
    counts = {p: int((pop == p).sum()) for p in ["P", "H"]}
    ndon = {p: int(len(set(donor[pop == p]))) for p in ["P", "H"]}
    meta = {"disease": disease_label, "source": gse, "cohort": "GEO",
            "cell_counts": counts, "n_donors": ndon, "n_candidate_genes": int(len(genes)),
            "reference_pop": "H (normal scar)", "top25": df.head(25).to_dict(orient="records")}
    with open(f"{DATA}/{out_name}_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=float)
    vol.commit()

    print(json.dumps({"gse": gse, "cell_counts": counts, "n_donors": ndon,
                      "n_candidate_genes": int(len(genes))}, indent=2))
    print("TOP 20 (GEO cohort):")
    for r in df.head(20).to_dict("records"):
        print(f"  {str(r['gene'])[:16]:16s} RACS={r['RACS']:.3f}  Sep={r['Sep']:.3f}  "
              f"log2FC={r.get('log2FC', float('nan')):.2f}  detect_P={r.get('detect_P', float('nan')):.0f}%")
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
