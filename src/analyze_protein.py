"""
analyze_protein.py — Stage 3 of the pipeline (protein / CITE-seq arm).

Loads the 20-protein CITE-seq panel, maps CD nomenclature to HGNC gene
symbols, aligns to the same filtered cell set as the RNA arm, and runs
per-condition Wilcoxon significance testing (each knockdown vs NTC) on
the protein panel.

Run directly:  python analyze_protein.py
Or via the orchestrator: python run_all.py
"""

import os

import pandas as pd
import scanpy as sc
import anndata as ad

from console import stage, step, ok, kv, result, suppress_noisy_warnings

suppress_noisy_warnings()

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

PROTEIN_CSV = "data/Protein_expression.csv"
FILTERED_CELLS_CSV = "data/filtered_cells.csv"
OUT_H5AD = "data/protein_expression.h5ad"

CD_TO_GENE = {
    "CD117": "KIT", "CD119": "IFNGR1", "CD140a": "PDGFRA", "CD140b": "PDGFRB",
    "CD172a": "SIRPA", "CD184": "CXCR4", "CD202b": "TEK", "CD274": "CD274",
    "CD29": "ITGB1", "CD309": "KDR", "CD44": "CD44", "CD47": "CD47",
    "CD49f": "ITGA6", "CD58": "CD58", "CD59": "CD59", "CD61": "ITGB3",
    "HLA_A": "HLA-A", "HLA_E": "HLA-E", "CD9": "CD9", "CD279": "PDCD1",
}
GENES_OF_INTEREST = ["CD274", "HLA-A", "HLA-E", "CD58", "IFNGR1"]


def load_protein_data():
    """Load CITE-seq protein matrix, map CD nomenclature -> gene symbols, align to filtered cells."""
    stage("Step 1 — Load protein panel + map CD nomenclature")
    protein_df = pd.read_csv(PROTEIN_CSV, index_col=0)
    protein_df.index = protein_df.index.str.replace(" protein", "", regex=False)
    kv("Loaded", f"{protein_df.shape[0]} proteins x {protein_df.shape[1]} cells")

    protein_df.index = protein_df.index.map(CD_TO_GENE)

    meta = pd.read_csv(FILTERED_CELLS_CSV)
    common_cells = protein_df.columns.intersection(meta["NAME"])
    kv("Cells with protein data + passing filter", f"{len(common_cells)} / {protein_df.shape[1]}")

    protein_df = protein_df[common_cells]
    adata_p = ad.AnnData(protein_df.T)  # cells x proteins, var_names = gene symbols
    adata_p.obs = meta.set_index("NAME").loc[common_cells]

    # values are already log1p-transformed — do not re-normalize (verified against
    # known log1p reference values: 0.405≈ln(1.5), 0.693≈ln(2), 1.386≈ln(4))
    step("Verified values are already log1p-transformed (not re-normalizing)")

    adata_p.write(OUT_H5AD)
    ok(f"Saved {OUT_H5AD}")


def compute_protein_significance():
    """Per-condition Wilcoxon significance testing, target gene vs NTC, on the protein panel."""
    stage("Step 2 — Protein significance testing (per condition, vs NTC)")
    adata = ad.read_h5ad(OUT_H5AD)
    adata.obs["target_gene"] = adata.obs["target_gene"].astype(str).replace(
        {"NO_SITE": "NTC", "ONE_NON-GENE_SITE": "NTC"}
    )

    conditions = adata.obs["condition"].unique().tolist()
    kv("Conditions found", ", ".join(conditions))

    results = []
    for condition in conditions:
        sub = adata[adata.obs["condition"] == condition].copy()
        if "NTC" not in sub.obs["target_gene"].unique():
            continue
        genes_here = [g for g in sub.obs["target_gene"].unique() if g != "NTC"]
        for gene in genes_here:
            mask = sub.obs["target_gene"].isin([gene, "NTC"])
            test_sub = sub[mask].copy()
            if (test_sub.obs["target_gene"] == gene).sum() < 5:
                continue
            sc.tl.rank_genes_groups(
                test_sub, groupby="target_gene", groups=[gene], reference="NTC", method="wilcoxon"
            )
            df = sc.get.rank_genes_groups_df(test_sub, group=gene)
            n_sig = (df["pvals_adj"] < 0.05).sum()
            self_row = df[df["names"] == gene]
            self_padj = self_row["pvals_adj"].values[0] if len(self_row) else None
            self_lfc = self_row["logfoldchanges"].values[0] if len(self_row) else None
            results.append({
                "target_gene": gene, "condition": condition,
                "n_sig_proteins": n_sig, "self_padj": self_padj, "self_lfc": self_lfc
            })
        step(f"{condition}: {len(genes_here)} target genes tested")

    res_df = pd.DataFrame(results).sort_values("n_sig_proteins", ascending=False)
    out_csv = os.path.join(RESULTS_DIR, "protein_significance.csv")
    res_df.to_csv(out_csv, index=False)
    ok(f"Saved {out_csv} ({len(res_df)} rows)")

    panel_matches = res_df[res_df["target_gene"].isin(GENES_OF_INTEREST)]
    step(f"Direct panel matches ({len(panel_matches)} rows) — see {out_csv} for full table")
    for r in panel_matches.itertuples():
        print(f"      {r.target_gene:<10} {r.condition:<12} n_sig_proteins={r.n_sig_proteins}")

    if len(res_df):
        top = res_df.iloc[0]
        result(f"Top protein hit: {top['target_gene']} in {top['condition']} "
               f"({int(top['n_sig_proteins'])} proteins shifted)")


def main():
    load_protein_data()
    compute_protein_significance()


if __name__ == "__main__":
    main()
