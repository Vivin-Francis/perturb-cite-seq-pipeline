"""
console.py — shared terminal output helpers, used by every stage script
so the pipeline reads the same way end to end instead of each script
printing in its own style.

Import pattern at the top of each stage script:

    from console import stage, step, ok, info, warn, result, suppress_noisy_warnings
    suppress_noisy_warnings()
"""

import warnings

_WIDTH = 60


def suppress_noisy_warnings():
    """Silence the routine FutureWarning/DeprecationWarning spam that
    scanpy/anndata/leidenalg emit on every run. Real errors still surface —
    this only filters warning-level noise."""
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)


def stage(title: str):
    """Top-level banner for a pipeline stage (one per script / major step)."""
    print(f"\n{'─' * _WIDTH}")
    print(f"  {title}")
    print(f"{'─' * _WIDTH}")


def step(msg: str):
    """A sub-step within a stage — quieter than stage()."""
    print(f"  · {msg}")


def ok(msg: str):
    print(f"  ✓ {msg}")


def warn(msg: str):
    print(f"  ! {msg}")


def info(msg: str):
    print(f"    {msg}")


def kv(label: str, value):
    """Aligned label: value line, for quick stats."""
    print(f"    {label:<32} {value}")


def result(msg: str):
    """One crisp finding/number at the end of a stage — the actual
    output of the analysis, not just a completion tick."""
    print(f"  → {msg}")
