"""
analyze_rna.py — Stage 2 of the pipeline (RNA arm).

Steps:
  1. Per-condition Wilcoxon significance testing (each knockdown vs NTC)
  2. Clustering of perturbations by their transcriptional signature
  3. Summary bar charts for known resistance-pathway genes
  4. Global cell-state clustering (PCA -> neighbors -> Leiden -> UMAP),
     with a focused NTC vs JAK/STAT-knockdown comparison per condition

Run directly:  python analyze_rna.py
Or via the orchestrator: python run_all.py
"""

import os
import time

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

from console import stage, step, ok, warn, kv, result, suppress_noisy_warnings

suppress_noisy_warnings()

IN_H5AD = "data/perturb_cite_seq_qc.h5ad"
PCA_CHECKPOINT = "data/perturb_cite_seq_pca.h5ad"
CLUSTERED_H5AD = "data/perturb_cite_seq_clustered.h5ad"
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# ---- significance testing thresholds ----
CONTROL_LABELS = ["NO_SITE", "ONE_NON-GENE_SITE"]
MIN_CELLS_PER_GROUP = 20
PADJ_THRESHOLD = 0.05

# ---- perturbation clustering thresholds ----
MIN_SIG_GENES_TO_INCLUDE = 5

# ---- plotting ----
KNOWN_RESISTANCE_GENES = [
    "JAK1", "JAK2", "STAT1", "STAT3", "IFNGR1", "IFNGR2",
    "B2M", "TAPBP", "TAPBPL", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "CD274",
]
CONDITION_ORDER = ["Control", "Co-culture", "IFNγ"]
CONDITION_COLORS = {"Control": "#999999", "Co-culture": "#4C72B0", "IFNγ": "#C44E52"}

# ---- cell state clustering ----
N_TOP_GENES = 2000
N_PCS = 30
JAK_STAT_GENES = ["JAK1", "JAK2", "STAT1", "IFNGR1", "IFNGR2"]


def compute_significance():
    """Per-condition Wilcoxon significance testing, target gene vs NTC."""
    stage("Step 1 — RNA significance testing (per condition, vs NTC)")
    t0 = time.time()
    adata = sc.read_h5ad(IN_H5AD)
    kv("Loaded", f"{adata.n_obs} cells x {adata.n_vars} genes")

    adata.obs["target_gene_merged"] = adata.obs["target_gene"].astype(str)
    adata.obs.loc[
        adata.obs["target_gene_merged"].isin(CONTROL_LABELS), "target_gene_merged"
    ] = "NTC"

    conditions = adata.obs["condition"].unique().tolist()
    kv("Conditions found", ", ".join(conditions))

    all_summary_rows = []
    all_sig_rows = []

    for cond in conditions:
        sub = adata[adata.obs["condition"] == cond].copy()

        group_counts = sub.obs["target_gene_merged"].value_counts()
        valid_groups = group_counts[group_counts >= MIN_CELLS_PER_GROUP].index.tolist()
        if "NTC" not in valid_groups:
            warn(f"{cond}: skipped — not enough NTC control cells")
            continue
        test_groups = [g for g in valid_groups if g != "NTC"]
        step(f"{cond}: testing {len(test_groups)} target genes against {group_counts['NTC']} NTC cells")

        sub = sub[sub.obs["target_gene_merged"].isin(valid_groups)].copy()

        sc.tl.rank_genes_groups(
            sub, groupby="target_gene_merged", groups=test_groups,
            reference="NTC", method="wilcoxon", tie_correct=False,
        )

        for gene in test_groups:
            result_df = sc.get.rank_genes_groups_df(sub, group=gene)
            n_sig = (result_df["pvals_adj"] < PADJ_THRESHOLD).sum()
            n_cells = int(group_counts[gene])

            all_summary_rows.append({
                "target_gene": gene, "condition": cond, "n_cells": n_cells,
                "n_sig_genes": n_sig, "n_genes_tested": len(result_df),
                "pct_sig_genes": 100 * n_sig / len(result_df) if len(result_df) else 0,
            })

            sig_df = result_df[result_df["pvals_adj"] < PADJ_THRESHOLD].copy()
            if len(sig_df) > 0:
                sig_df["target_gene"] = gene
                sig_df["condition"] = cond
                all_sig_rows.append(sig_df)

    summary = pd.DataFrame(all_summary_rows)
    summary_path = os.path.join(RESULTS_DIR, "significance_summary.csv")
    summary.to_csv(summary_path, index=False)
    ok(f"Saved {summary_path} ({len(summary)} rows)")

    if all_sig_rows:
        sig_long = pd.concat(all_sig_rows, ignore_index=True)
        sig_long_path = os.path.join(RESULTS_DIR, "significance_long.csv")
        sig_long.to_csv(sig_long_path, index=False)
        ok(f"Saved {sig_long_path} ({len(sig_long)} rows)")

    top5 = summary.sort_values("n_sig_genes", ascending=False).head(5)
    step("Top 5 (target_gene, condition) pairs by impact:")
    for r in top5.itertuples():
        print(f"      {r.target_gene:<10} {r.condition:<12} n_sig_genes={r.n_sig_genes}")

    kv("Runtime", f"{time.time() - t0:.0f}s")
    top = top5.iloc[0]
    result(f"Strongest hit: {top['target_gene']} in {top['condition']} ({int(top['n_sig_genes'])} genes shifted)")

    return adata


def cluster_perturbations():
    """Cluster (gene, condition) pairs by similarity of their significant-hit profiles."""
    stage("Step 2 — Cluster perturbations by transcriptional signature")
    long_df = pd.read_csv(os.path.join(RESULTS_DIR, "significance_long.csv"))
    summary_df = pd.read_csv(os.path.join(RESULTS_DIR, "significance_summary.csv"))

    valid_pairs = summary_df[summary_df["n_sig_genes"] >= MIN_SIG_GENES_TO_INCLUDE]
    kv(f"Pairs with >= {MIN_SIG_GENES_TO_INCLUDE} significant genes", len(valid_pairs))
    if len(valid_pairs) < 5:
        raise ValueError(
            "Too few (target_gene, condition) pairs have enough significant genes to cluster. "
            "Lower MIN_SIG_GENES_TO_INCLUDE or check the significance results."
        )

    valid_pairs = valid_pairs.copy()
    valid_pairs["pair_id"] = valid_pairs["target_gene"] + " | " + valid_pairs["condition"]
    long_df["pair_id"] = long_df["target_gene"] + " | " + long_df["condition"]
    long_df = long_df[long_df["pair_id"].isin(valid_pairs["pair_id"])]

    # 'scores' is the signed Wilcoxon test statistic — captures direction and
    # strength of each measured gene's shift, not just significance
    matrix = long_df.pivot_table(index="pair_id", columns="names", values="scores", fill_value=0)
    kv("Signature matrix (pairs x measured genes)", str(matrix.shape))

    pair_adata = ad.AnnData(X=matrix.values)
    pair_adata.obs_names = matrix.index.tolist()
    pair_adata.var_names = matrix.columns.tolist()

    meta = valid_pairs.set_index("pair_id").loc[matrix.index]
    pair_adata.obs["target_gene"] = meta["target_gene"].values
    pair_adata.obs["condition"] = meta["condition"].values
    pair_adata.obs["n_sig_genes"] = meta["n_sig_genes"].values

    n_pcs = min(30, pair_adata.n_obs - 1, pair_adata.n_vars - 1)
    sc.pp.pca(pair_adata, n_comps=n_pcs)
    n_neighbors = min(15, pair_adata.n_obs - 1)
    sc.pp.neighbors(pair_adata, n_neighbors=n_neighbors)
    sc.tl.umap(pair_adata)
    sc.tl.leiden(pair_adata, resolution=1.0)

    kv("Leiden clusters found", pair_adata.obs["leiden"].nunique())

    out_df = pair_adata.obs.reset_index().rename(columns={"index": "pair_id"})
    out_csv = os.path.join(RESULTS_DIR, "perturbation_clusters.csv")
    out_df.to_csv(out_csv, index=False)
    ok(f"Saved {out_csv}")

    step("Members by cluster:")
    for cluster in sorted(out_df["leiden"].unique(), key=int):
        members = out_df[out_df["leiden"] == cluster]
        gene_list = ", ".join(f"{r.target_gene} ({r.condition})" for r in members.itertuples())
        print(f"      Cluster {cluster} ({len(members)}): {gene_list}")

    # NOTE: this UMAP is at the (gene, condition) PAIR level, not per-cell —
    # one point per perturbation signature, so a low point count here is
    # expected. Per-cell UMAP lives in umap_cell_state_by_condition.png.
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, color_by, title in [
        (axes[0], "leiden", "Colored by Leiden cluster"),
        (axes[1], "condition", "Colored by condition"),
    ]:
        sc.pl.umap(pair_adata, color=color_by, ax=ax, show=False, title=title,
                   size=200, legend_fontsize=9)
    fig.subplots_adjust(wspace=0.5)
    plt.tight_layout()
    out_png = os.path.join(RESULTS_DIR, "perturbation_umap.png")
    plt.savefig(out_png, dpi=150)
    plt.close(fig)
    ok(f"Saved {out_png}")

    cross = pd.crosstab(out_df["leiden"], out_df["condition"])
    single_condition_clusters = (cross.gt(0).sum(axis=1) == 1).sum()
    if single_condition_clusters:
        warn(f"{single_condition_clusters} cluster(s) are 100% one condition — worth a closer look")

    result(f"{pair_adata.obs['leiden'].nunique()} perturbation clusters found from "
           f"{pair_adata.n_obs} gene×condition pairs")


def visualize_results():
    """Summary bar charts: known genes by condition, and top-15 overall."""
    stage("Step 3 — Summary plots")
    df = pd.read_csv(os.path.join(RESULTS_DIR, "significance_summary.csv"))

    known = df[df["target_gene"].isin(KNOWN_RESISTANCE_GENES)].copy()
    pivot = known.pivot(index="target_gene", columns="condition", values="n_sig_genes").fillna(0)
    pivot = pivot.reindex(columns=[c for c in CONDITION_ORDER if c in pivot.columns])
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(pivot.index))
    width = 0.25
    for i, cond in enumerate(pivot.columns):
        ax.bar(x + i * width, pivot[cond], width, label=cond, color=CONDITION_COLORS.get(cond, None))
    ax.set_xticks(x + width)
    ax.set_xticklabels(pivot.index, rotation=45, ha="right")
    ax.set_ylabel("Number of significantly shifted genes\n(adjusted p < 0.05)")
    ax.set_title("Known resistance-pathway genes: knockdown effect strength by condition")
    ax.legend(title="Condition")
    plt.tight_layout()
    out1 = os.path.join(RESULTS_DIR, "known_genes_by_condition.png")
    plt.savefig(out1, dpi=150)
    plt.close(fig)
    ok(f"Saved {out1}")

    top15 = df.sort_values("n_sig_genes", ascending=False).head(15).copy()
    top15["label"] = top15["target_gene"] + " (" + top15["condition"] + ")"

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = [CONDITION_COLORS.get(c, "#333333") for c in top15["condition"]]
    ax.barh(top15["label"][::-1], top15["n_sig_genes"][::-1], color=colors[::-1])
    ax.set_xlabel("Number of significantly shifted genes (adjusted p < 0.05)")
    ax.set_title("Top 15 knockdown x condition pairs by transcriptional impact")
    plt.tight_layout()
    out2 = os.path.join(RESULTS_DIR, "top15_overall.png")
    plt.savefig(out2, dpi=150)
    plt.close(fig)
    ok(f"Saved {out2}")

    top_gene = pivot.max(axis=1).idxmax()
    result(f"Strongest known resistance-gene signal: {top_gene}")


def cell_state_clustering(adata):
    """Global cell-state clustering (PCA -> neighbors -> Leiden -> UMAP),
    then a focused NTC vs JAK/STAT-knockdown comparison per condition."""
    stage("Step 4 — Global cell-state clustering")
    adata.obs["target_gene_merged"] = adata.obs["target_gene"].astype(str).replace(
        {"NO_SITE": "NTC", "ONE_NON-GENE_SITE": "NTC"}
    )
    kv("Starting cells", adata.n_obs)

    step("Running PCA on top HVGs")
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor="seurat")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=N_PCS)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.write(PCA_CHECKPOINT)
    ok(f"Checkpoint saved to {PCA_CHECKPOINT} — resume from here if anything crashes")

    step("Neighbors + Leiden + UMAP")
    # leidenalg flavor (default) is deliberate here, NOT "igraph" — switching
    # flavors changes cluster boundaries/numbering, breaking exact
    # reproducibility of downstream cluster-identity findings.
    sc.pp.neighbors(adata, n_pcs=N_PCS)
    sc.tl.leiden(adata, resolution=1.0, key_added="cell_state")
    sc.tl.umap(adata)
    adata.write(CLUSTERED_H5AD)
    ok(f"{adata.obs['cell_state'].nunique()} cell states found -> {CLUSTERED_H5AD}")

    step("Cell state composition by condition:")
    comp = pd.crosstab(adata.obs["cell_state"], adata.obs["condition"], normalize="columns").round(3)
    print(comp.to_string())

    adata.obs["group"] = "other"
    adata.obs.loc[adata.obs["target_gene_merged"] == "NTC", "group"] = "NTC"
    adata.obs.loc[adata.obs["target_gene_merged"].isin(JAK_STAT_GENES), "group"] = "JAK_STAT_knockdown"

    for cond in adata.obs["condition"].unique():
        sub = adata.obs[(adata.obs["condition"] == cond) & (adata.obs["group"] != "other")]
        step(f"{cond} — cell state by group (NTC vs JAK/STAT knockdown):")
        print(pd.crosstab(sub["cell_state"], sub["group"], normalize="columns").round(3).to_string())

    # Wider spacing (wspace) + smaller legend font stop the 10-category
    # cell_state legend from overlapping the neighboring panel's UMAP2 label.
    sc.settings.set_figure_params(figsize=(6, 5))
    fig = sc.pl.umap(
        adata, color=["condition", "cell_state", "group"], ncols=3,
        show=False, return_fig=True, wspace=0.6, legend_fontsize=9,
    )
    out_png = os.path.join(RESULTS_DIR, "umap_cell_state_by_condition.png")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    ok(f"Saved {out_png}")

    result(f"{adata.obs['cell_state'].nunique()} cell states identified across {adata.n_obs:,} cells")


def main():
    t0 = time.time()
    adata = compute_significance()
    cluster_perturbations()
    visualize_results()
    cell_state_clustering(adata)
    kv("analyze_rna.py runtime", f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
