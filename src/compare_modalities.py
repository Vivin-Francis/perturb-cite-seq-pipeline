"""
compare_modalities.py — Stage 4 of the pipeline.

Cross-references RNA-level and protein-level significance for the genes
present on the CITE-seq panel, flagging cases where the two modalities
disagree (e.g. RNA-silent but protein-significant) — the project's
clearest original finding.

Run directly:  python compare_modalities.py
Or via the orchestrator: python run_all.py
"""

import os

import pandas as pd

from console import stage, step, ok, kv, result, suppress_noisy_warnings

suppress_noisy_warnings()

RESULTS_DIR = "results"
GENES_IN_PANEL = ["CD274", "HLA-A", "HLA-E", "CD58", "IFNGR1"]


def main():
    stage("Compare RNA vs protein significance across the CITE-seq panel")

    summary = pd.read_csv(os.path.join(RESULTS_DIR, "significance_summary.csv"))
    protein = pd.read_csv(os.path.join(RESULTS_DIR, "protein_significance.csv"))
    sig_long = pd.read_csv(os.path.join(RESULTS_DIR, "significance_long.csv"))

    rna_self = summary[summary["target_gene"].isin(GENES_IN_PANEL)][
        ["target_gene", "condition", "n_cells", "n_sig_genes", "n_genes_tested"]
    ].copy()

    self_hits = sig_long[sig_long["names"] == sig_long["target_gene"]][
        ["target_gene", "condition", "pvals_adj", "logfoldchanges"]
    ].rename(columns={"pvals_adj": "rna_self_padj", "logfoldchanges": "rna_self_lfc"})

    rna_self = rna_self.merge(self_hits, on=["target_gene", "condition"], how="left")
    rna_self["rna_self_sig"] = rna_self["rna_self_padj"].notna()

    prot_self = protein[protein["target_gene"].isin(GENES_IN_PANEL)][
        ["target_gene", "condition", "self_padj", "self_lfc"]
    ]

    merged = rna_self.merge(prot_self, on=["target_gene", "condition"], how="left")
    merged["protein_sig"] = merged["self_padj"] < 0.05
    merged["discordant"] = merged["rna_self_sig"] != merged["protein_sig"]
    merged = merged.sort_values(["target_gene", "condition"])

    out_csv = os.path.join(RESULTS_DIR, "rna_protein_discordance.csv")
    merged.to_csv(out_csv, index=False)
    ok(f"Saved {out_csv} ({len(merged)} rows)")

    kv("Panel genes compared", ", ".join(GENES_IN_PANEL))
    step("RNA vs protein significance by gene x condition:")
    print(merged[["target_gene", "condition", "n_cells", "rna_self_sig", "protein_sig", "discordant"]]
          .to_string(index=False))

    n_discordant = int(merged["discordant"].sum())
    kv("Discordant (RNA and protein disagree)", f"{n_discordant} / {len(merged)}")

    discordant_genes = sorted(merged.loc[merged["discordant"], "target_gene"].unique())
    if discordant_genes:
        step(f"Genes with RNA/protein discordance: {', '.join(discordant_genes)}")

    result(f"{n_discordant}/{len(merged)} pairs discordant"
           + (f": {', '.join(discordant_genes)}" if discordant_genes else ""))


if __name__ == "__main__":
    main()
