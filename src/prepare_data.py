"""
prepare_data.py — Stage 1 of the pipeline.

Filters cells (MOI=1, guide-count threshold), extracts the matching columns
from the raw expression matrix, loads them into a sparse AnnData object, and
applies standard QC + normalization (CP10K + log1p).

Run directly:  python prepare_data.py
Or via the orchestrator: python run_all.py
"""

import os
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt

from console import stage, step, ok, warn, kv, result, suppress_noisy_warnings

suppress_noisy_warnings()

# ---- shared paths ----
RNA_FULL_CSV = "data/RNA_expression_full.csv"
RNA_META_CSV = "data/RNA_metadata.csv"
FILTERED_CELLS_CSV = "data/filtered_cells.csv"
WANTED_COLS_TXT = "data/wanted_columns.txt"
EXPR_FILTERED_CSV = "data/RNA_expression_filtered.csv"
RAW_H5AD = "data/perturb_cite_seq.h5ad"
QC_H5AD = "data/perturb_cite_seq_qc.h5ad"
PLOT_DIR = "results"
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs("data", exist_ok=True)

# ---- step 1 thresholds ----
MIN_CELLS_PER_GENE = 50

# ---- step 4 thresholds ----
CHUNKSIZE = 500  # genes (rows) per chunk; tune down if memory pressure appears

# ---- step 5 thresholds ----
MIN_GENES_PER_CELL = 200
MIN_COUNTS_PER_CELL = 500
MAX_PCT_MT = 20.0


def step1_filter_cells():
    """MOI=1 filtering + require >=50 cells per target gene."""
    stage("Step 1 — Filter cells (MOI=1, min cells per guide)")
    df = pd.read_csv(RNA_META_CSV, skiprows=[1])
    moi1 = df[df["MOI"] == 1].copy()
    moi1["target_gene"] = moi1["sgRNA"].str.rsplit("_", n=1).str[0]

    gene_counts = moi1.groupby("target_gene").size()
    keep_genes = gene_counts[gene_counts >= MIN_CELLS_PER_GENE].index

    filtered = moi1[moi1["target_gene"].isin(keep_genes)]
    filtered[["NAME", "condition", "sgRNA", "target_gene"]].to_csv(FILTERED_CELLS_CSV, index=False)

    kv("Cells at MOI=1", len(moi1))
    kv("Target genes kept (>= {} cells)".format(MIN_CELLS_PER_GENE), len(keep_genes))
    ok(f"{len(filtered)} cells kept -> {FILTERED_CELLS_CSV}")


def step2_get_column_positions():
    """Map filtered cell names to their column positions in the raw matrix, for awk."""
    stage("Step 2 — Locate matching columns in raw matrix")
    with open(RNA_FULL_CSV) as f:
        header = f.readline().strip().split(",")
    col_position = {name: i + 1 for i, name in enumerate(header)}  # header[0] = "GENE"

    filtered = pd.read_csv(FILTERED_CELLS_CSV)
    wanted_cells = filtered["NAME"].tolist()

    positions = [1] + [col_position[c] for c in wanted_cells if c in col_position]
    kv("Columns to extract (incl. GENE)", len(positions))

    with open(WANTED_COLS_TXT, "w") as f:
        f.write(",".join(str(p) for p in positions))
    ok(f"Saved {WANTED_COLS_TXT}")


def step3_extract_columns():
    """Stream the raw matrix through extract.awk, keeping only wanted columns."""
    stage("Step 3 — Extract columns from raw matrix (awk streaming)")
    step("This can take a while....")
    t0 = time.time()
    with open(EXPR_FILTERED_CSV, "w") as out_f:
        __import__("subprocess").run(["awk", "-f", "extract.awk", RNA_FULL_CSV], stdout=out_f, check=True)
    ok(f"Extraction complete in {time.time() - t0:.0f}s -> {EXPR_FILTERED_CSV}")


def step4_load_to_anndata():
    """Chunked CSV -> sparse AnnData."""
    stage("Step 4 — Load filtered matrix into AnnData")
    t0 = time.time()

    with open(EXPR_FILTERED_CSV, "r") as f:
        header = f.readline().strip().split(",")
    gene_col_name, cell_names = header[0], header[1:]
    kv("Cells", len(cell_names))

    gene_names = []
    sparse_blocks = []

    # Dtype map applies float32 only to cell-value columns, leaving the GENE
    # column untyped — a scalar dtype passed to read_csv applies BEFORE
    # index_col splits off the index, which fails on gene name strings.
    dtype_map = {c: np.float32 for c in cell_names}
    reader = pd.read_csv(EXPR_FILTERED_CSV, chunksize=CHUNKSIZE, index_col=0, dtype=dtype_map)

    n_chunks = 0
    for chunk in reader:
        gene_names.extend(chunk.index.tolist())
        sparse_blocks.append(sp.csr_matrix(chunk.values))
        n_chunks += 1
        if n_chunks % 20 == 0:
            step(f"{len(gene_names)} genes read so far ({time.time() - t0:.0f}s elapsed)")

    expr_genes_by_cells = sp.vstack(sparse_blocks, format="csr")
    expr_cells_by_genes = expr_genes_by_cells.transpose().tocsr()
    kv("Assembled matrix (cells x genes)", str(expr_cells_by_genes.shape))

    meta = pd.read_csv(FILTERED_CELLS_CSV).set_index("NAME").loc[cell_names]
    adata = ad.AnnData(X=expr_cells_by_genes, obs=meta, var=pd.DataFrame(index=gene_names))

    adata.write_h5ad(RAW_H5AD)
    ok(f"Saved {RAW_H5AD} ({time.time() - t0:.0f}s total)")


def step5_qc_normalize():
    """QC filtering (genes/counts/%MT) + CP10K + log1p normalization."""
    stage("Step 5 — QC filtering + normalization")
    adata = sc.read_h5ad(RAW_H5AD)

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt_genes = int(adata.var["mt"].sum())
    kv("Mitochondrial genes found", n_mt_genes)

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].hist(adata.obs["total_counts"], bins=100)
    axes[0].axvline(MIN_COUNTS_PER_CELL, color="red", linestyle="--", label="QC cutoff")
    axes[0].set_title("Total counts per cell")
    axes[0].set_xlabel("total_counts")
    axes[0].legend(fontsize=8)

    axes[1].hist(adata.obs["n_genes_by_counts"], bins=100)
    axes[1].axvline(MIN_GENES_PER_CELL, color="red", linestyle="--", label="QC cutoff")
    axes[1].set_title("Genes detected per cell")
    axes[1].set_xlabel("n_genes_by_counts")
    axes[1].legend(fontsize=8)

    if n_mt_genes > 0:
        axes[2].hist(adata.obs["pct_counts_mt"], bins=100)
        axes[2].axvline(MAX_PCT_MT, color="red", linestyle="--", label="QC cutoff")
        axes[2].set_title("% mitochondrial counts")
        axes[2].set_xlabel("pct_counts_mt")
        axes[2].legend(fontsize=8)
    else:
        axes[2].text(0.5, 0.5, "No MT- genes found\nin this gene panel", ha="center", va="center")
        axes[2].set_title("% mitochondrial counts")

    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "qc_plots.png")
    plt.savefig(plot_path, dpi=150)
    plt.close(fig)
    ok(f"Saved {plot_path}")

    # Single combined mask, one copy — avoids the extra intermediate copies
    # that previously caused an OOM kill on this machine.
    n_before = adata.n_obs
    mask = (
        (adata.obs["n_genes_by_counts"] >= MIN_GENES_PER_CELL)
        & (adata.obs["total_counts"] >= MIN_COUNTS_PER_CELL)
    )
    if n_mt_genes > 0:
        mask &= adata.obs["pct_counts_mt"] < MAX_PCT_MT

    adata = adata[mask].copy()
    n_after = adata.n_obs
    kv("Cells before QC", n_before)
    kv("Cells after QC", f"{n_after}  (-{n_before - n_after}, {100 * (n_before - n_after) / n_before:.1f}%)")

    genes_below_threshold = int((adata.obs["target_gene"].value_counts() < MIN_CELLS_PER_GENE).sum())
    if genes_below_threshold:
        warn(f"{genes_below_threshold} target genes fell below the original "
             f"{MIN_CELLS_PER_GENE}-cell threshold after QC")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata.write_h5ad(QC_H5AD)
    ok(f"Saved QC'd + normalized data -> {QC_H5AD}")
    kv("Final shape", str(adata.shape))
    result(f"{n_after:,}/{n_before:,} cells retained after QC ({100 * n_after / n_before:.1f}%)")


def main():
    t0 = time.time()
    step1_filter_cells()
    step2_get_column_positions()
    step3_extract_columns()
    step4_load_to_anndata()
    step5_qc_normalize()
    kv("prepare_data.py runtime", f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
