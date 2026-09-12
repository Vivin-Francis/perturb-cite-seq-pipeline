"""
run_all.py — runs the entire pipeline end to end, in the correct order,
stopping immediately if any stage fails.

Usage:
    python run_all.py

That's it — one command runs data prep, RNA analysis, protein analysis,
modality comparison, and validation checks.

Expected layout (run this script from the project root):

    project_root/
    ├── run_all.py          <- you are here
    ├── requirements.txt
    ├── extract.awk
    ├── src/
    │   ├── console.py
    │   ├── prepare_data.py
    │   ├── analyze_rna.py
    │   ├── analyze_protein.py
    │   ├── compare_modalities.py
    │   └── validation/
    │       ├── guide_replication_check.py
    │       ├── ntc_baseline_shift.py
    │       └── cluster0_enrichment.py
    ├── data/
    └── results/
"""

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
VALIDATION = SRC / "validation"

sys.path.insert(0, str(SRC))
from console import stage, ok, warn, kv, _WIDTH  # noqa: E402

STAGES = [
    ("prepare_data", SRC / "prepare_data.py"),
    ("analyze_rna", SRC / "analyze_rna.py"),
    ("analyze_protein", SRC / "analyze_protein.py"),
    ("compare_modalities", SRC / "compare_modalities.py"),
    ("guide_replication_check", VALIDATION / "guide_replication_check.py"),
    ("ntc_baseline_shift", VALIDATION / "ntc_baseline_shift.py"),
    ("cluster0_enrichment", VALIDATION / "cluster0_enrichment.py"),
]


def run_stage(name, script_path, index, total):
    stage(f"STAGE {index}/{total}: {name}")
    if not script_path.exists():
        warn(f"{script_path} not found. Aborting.")
        sys.exit(1)

    t0 = time.time()
    # run with cwd = project root, since every script uses paths like
    # "data/..." and "results/..." relative to the project root
    result = subprocess.run([sys.executable, str(script_path)], cwd=str(ROOT))
    elapsed = time.time() - t0

    if result.returncode != 0:
        warn(f"{name} exited with code {result.returncode} after "
             f"{elapsed:.0f}s. Stopping pipeline here.")
        sys.exit(result.returncode)

    ok(f"{name} finished in {elapsed:.0f}s")
    return elapsed


def main():
    total = len(STAGES)
    stage(f"PIPELINE START — {total} stages")
    t_total = time.time()

    timings = []
    for i, (name, path) in enumerate(STAGES, start=1):
        elapsed = run_stage(name, path, i, total)
        timings.append((name, elapsed))

    stage("PIPELINE COMPLETE")
    for name, elapsed in timings:
        kv(name, f"{elapsed:.0f}s")
    print(f"{'─' * _WIDTH}")
    kv("TOTAL", f"{time.time() - t_total:.0f}s")


if __name__ == "__main__":
    main()
