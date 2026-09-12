"""
cluster0_enrichment.py — Validation script.

Functional enrichment (GO Biological Process, Reactome) of the genes
driving perturbation cluster 0, to check whether it reflects a coherent
biological pathway or a generic essentiality/housekeeping signature.

Run directly:  python validation/cluster0_enrichment.py
Or via the orchestrator: python run_all.py --with-validation
"""

import os
import sys

import gseapy as gp
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from console import stage, ok, kv, result, suppress_noisy_warnings

suppress_noisy_warnings()

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

CLUSTER0_GENES = [
    "CCND1", "CDK6", "CTPS1", "EIF3K", "MRPL47", "PABPC1",
    "TXNDC17", "TMED10", "NPC1", "NONO", "LY96",
]


def main():
    stage("Validation — Cluster 0 functional enrichment")
    kv("Genes tested", len(CLUSTER0_GENES))

    enr = gp.enrichr(
        gene_list=CLUSTER0_GENES,
        gene_sets=["GO_Biological_Process_2023", "Reactome_2022"],
        outdir=os.path.join(RESULTS_DIR, "cluster0_enrichment"),
    )
    results = enr.results.sort_values("Adjusted P-value")
    out_csv = os.path.join(RESULTS_DIR, "cluster0_enrichment.csv")
    results.to_csv(out_csv, index=False)
    ok(f"Saved {out_csv} ({len(results)} terms)")

    top = results.head(10)
    print(top[["Gene_set", "Term", "Overlap", "Adjusted P-value", "Genes"]].to_string(index=False))

    if len(results):
        top1 = results.iloc[0]
        result(f"Top enriched term: {top1['Term']} (adj p={top1['Adjusted P-value']:.2e})")
    else:
        result("No enriched terms passed threshold — cluster 0 does not resolve to a known pathway")


if __name__ == "__main__":
    main()
