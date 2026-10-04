#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 00_paths.py
#
# Repo-relative path resolution for the scripts that run on a clone.
#
# The analysis scripts were written against the original working tree, whose
# layout differs from this repository's. Rather than rewrite every path inside
# them (and risk silently changing which file a figure reads), this module maps
# the original relative paths onto their location here. The scripts keep their
# original expressions; only the root they resolve against changes.
#
# Only the CPU scripts that a cloner can actually run use this (97, 98, 99,
# 99b). The GPU scripts (80-96) address /workspace on a rented instance and are
# left exactly as they were executed.
#
# CONSUMES: nothing
# PRODUCES: nothing (module)
# ---------------------------------------------------------------------------
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# original working-tree location  ->  location in this repository
_MAP = [
    ("replicative_ceiling_test", "results/replicative_ceiling"),
    ("enhancer_test", "results/enhancer_cross_mechanism"),
    ("ois_enhancer_trainset", "results/ois_trainset"),
    ("ois_enhancer_run/results", "results/ois_model"),
    ("ois_enhancer_run/ois_cache", "results/ois_model"),
    ("ois_gates/results", "results/design_gates"),
    ("audit_fixes", "audit_fixes"),
    ("number_audit", "audit"),
    ("figures_manuscript", "figures"),
    ("infra_borzoi_grelu/execution_results_permutations", "results/expression_phase"),
]


def repo_path(rel) -> Path:
    """Map an original-tree relative path onto this repository."""
    rel = str(rel).lstrip("/")
    if rel.startswith("output/"):
        rel = rel[len("output/"):]
    for old, new in _MAP:
        if rel == old:
            return REPO / new
        if rel.startswith(old + "/"):
            return REPO / new / rel[len(old) + 1:]
    return REPO / rel


class _Dir:
    """Stands in for a directory so `BASE / "legacy/rel/path"` still works."""

    def __init__(self, base: str = ""):
        self.base = base

    def __truediv__(self, rel) -> Path:
        return repo_path(f"{self.base}/{rel}".strip("/"))

    def __fspath__(self):
        return str(REPO / self.base) if self.base else str(REPO)

    def __str__(self):
        return self.__fspath__()


ROOT = _Dir()
OUTPUT = _Dir("output")
FIGURES = REPO / "figures"
