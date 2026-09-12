"""
guide_replication_check.py — Validation script.

For each top resistance-pathway gene, checks that independent sgRNAs
targeting it produce consistent transcriptional effects, ruling out
single-guide artifacts as the source of a significance pattern.

Run directly:  python validation/guide_replication_check.py
Or via the orchestrator: python run_all.py --with-validation
"""

import os
import sys

import scanpy as sc
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from console import stage, step, ok, warn, kv, result, suppress_noisy_warnings

suppress_noisy_warnings()

IN_H5AD = "data/perturb_cite_seq_qc.h5ad"
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

MIN_CELLS_PER_GUIDE = 25  # lower bar than the 50 used for gene-level testing
TOP_GENES = ["JAK1", "JAK2", "STAT1", "IFNGR1", "IFNGR2"]


def main():
    stage("Validation — Guide-level replication check")
    adata = sc.read_h5ad(IN_H5AD)
    adata.obs["target_gene_merged"] = adata.obs["target_gene"].astype(str).replace(
        {"NO_SITE": "NTC", "ONE_NON-GENE_SITE": "NTC"}
    )
    kv("Genes checked", ", ".join(TOP_GENES))

    results = []
    for condition in adata.obs["condition"].unique():
        sub_cond = adata[adata.obs["condition"] == condition].copy()
        if (sub_cond.obs["target_gene_merged"] == "NTC").sum() < 25:
            continue

        for gene in TOP_GENES:
            guides = sub_cond.obs.loc[sub_cond.obs["target_gene"] == gene, "sgRNA"].astype(str).unique()
            for guide in guides:
                mask = sub_cond.obs["sgRNA"].astype(str).eq(guide) | (sub_cond.obs["target_gene_merged"] == "NTC")
                test_sub = sub_cond[mask].copy()
                test_sub.obs["guide_group"] = test_sub.obs["sgRNA"].astype(str).where(
                    test_sub.obs["sgRNA"].astype(str) == guide, "NTC"
                )
                n_cells = (test_sub.obs["guide_group"] == guide).sum()

                sc.tl.rank_genes_groups(
                    test_sub, groupby="guide_group", groups=[guide], reference="NTC", method="wilcoxon"
                )
                df = sc.get.rank_genes_groups_df(test_sub, group=guide)
                n_sig = (df["pvals_adj"] < 0.05).sum()
                self_row = df[df["names"] == gene]

                results.append({
                    "target_gene": gene, "guide": guide, "condition": condition,
                    "n_cells": n_cells, "underpowered": n_cells < MIN_CELLS_PER_GUIDE,
                    "n_sig_genes": n_sig,
                    "self_padj": self_row["pvals_adj"].values[0] if len(self_row) else None,
                    "self_lfc": self_row["logfoldchanges"].values[0] if len(self_row) else None,
                })
        step(f"{condition}: guides checked for all target genes")

    res_df = pd.DataFrame(results).sort_values(["target_gene", "condition", "guide"])
    out_csv = os.path.join(RESULTS_DIR, "guide_replication_check.csv")
    res_df.to_csv(out_csv, index=False)
    ok(f"Saved {out_csv} ({len(res_df)} rows)")

    step("Within-gene, within-condition consistency check:")
    flagged_pairs = []
    for (gene, cond), grp in res_df.groupby(["target_gene", "condition"]):
        if len(grp) < 2:
            continue
        n_sig_vals = grp["n_sig_genes"].values
        ratio = n_sig_vals.max() / max(n_sig_vals.min(), 1)
        flagged = ratio > 5 and n_sig_vals.max() > 20
        if flagged:
            flagged_pairs.append(f"{gene} ({cond})")
        marker = "FLAG — one guide much weaker" if flagged else "consistent"
        print(f"      {gene:8s} {cond:12s} n_sig per guide: {list(n_sig_vals)}  [{marker}]")

    if not flagged_pairs:
        result("All target genes replicate consistently across independent guides")
    else:
        result(f"{len(flagged_pairs)} pair(s) flagged for guide inconsistency: {', '.join(flagged_pairs)}")


if __name__ == "__main__":
    main()
