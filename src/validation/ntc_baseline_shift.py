"""
ntc_baseline_shift.py — Validation script.

Checks how much non-targeting-control (NTC) cells shift transcriptionally
just from condition exposure alone (Control vs Co-culture vs IFN-gamma),
confirming that downstream comparisons must be made against
condition-matched controls rather than a pooled baseline.

Run directly:  python validation/ntc_baseline_shift.py
Or via the orchestrator: python run_all.py --with-validation
"""

import os
import sys

import scanpy as sc
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from console import stage, step, ok, kv, result, suppress_noisy_warnings

suppress_noisy_warnings()

IN_H5AD = "data/perturb_cite_seq_qc.h5ad"
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)
JAK_STAT_GENES = ["JAK1", "JAK2", "STAT1", "IFNGR1", "IFNGR2"]


def main():
    stage("Validation — NTC baseline shift across conditions")
    adata = sc.read_h5ad(IN_H5AD)
    adata.obs["target_gene_merged"] = adata.obs["target_gene"].astype(str).replace(
        {"NO_SITE": "NTC", "ONE_NON-GENE_SITE": "NTC"}
    )

    ntc = adata[adata.obs["target_gene_merged"] == "NTC"].copy()
    kv("NTC cells by condition", dict(ntc.obs["condition"].value_counts()))

    results = []
    shift_counts = {}
    pairs = [("Co-culture", "Control"), ("IFNγ", "Control"), ("Co-culture", "IFNγ")]
    for cond_a, cond_b in pairs:
        sub = ntc[ntc.obs["condition"].isin([cond_a, cond_b])].copy()
        sc.tl.rank_genes_groups(sub, groupby="condition", groups=[cond_a], reference=cond_b, method="wilcoxon")
        df = sc.get.rank_genes_groups_df(sub, group=cond_a)
        n_sig = (df["pvals_adj"] < 0.05).sum()
        df["comparison"] = f"{cond_a}_vs_{cond_b}"
        results.append(df)
        shift_counts[f"{cond_a} vs {cond_b}"] = int(n_sig)
        step(f"{cond_a} vs {cond_b} (NTC cells only): {n_sig} significant genes")

    full = pd.concat(results)
    out_csv = os.path.join(RESULTS_DIR, "ntc_baseline_shift.csv")
    full.to_csv(out_csv, index=False)
    ok(f"Saved {out_csv}")

    step("JAK/STAT genes in NTC cells, across conditions:")
    jak_stat_rows = full[full["names"].isin(JAK_STAT_GENES)][
        ["names", "comparison", "logfoldchanges", "pvals_adj"]
    ]
    print(jak_stat_rows.to_string(index=False))

    biggest = max(shift_counts, key=shift_counts.get)
    result(f"Largest condition-driven NTC shift: {biggest} ({shift_counts[biggest]} genes)")


if __name__ == "__main__":
    main()
