#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 63_shared_denominator_null.py
#
# Constructed nulls showing a shared control arm manufactures correlation.
#
# CONSUMES: output/replicative_ceiling_test/signal_matrix_replicative.csv
# PRODUCES: output/replicative_ceiling_test/shared_denominator_null.txt
# ---------------------------------------------------------------------------
"""Null test for the shared-denominator artifact in the IR-vs-replicative comparison.

Both GSE106146 responses are computed against the SAME single proliferating control, so
correlating (IR - Pro) with (Sen - Pro) correlates the shared -Pro term. This quantifies
how much correlation that alone produces, using numerators with no shared biology.
"""
import pandas as pd, numpy as np
from scipy import stats
P = 0.10
sm = pd.read_csv('output/replicative_ceiling_test/signal_matrix_replicative.csv')
n = sm.copy()
for c in n.columns:
    n[c] = n[c] / np.nanmean(n[c].to_numpy())
def lfc(on, off):
    return np.log2((n[on].mean(axis=1).to_numpy() + P) / (n[off].mean(axis=1).to_numpy() + P))
rep   = lfc(['Sen2019_Sen'], ['Sen2019_Pro'])
ir    = lfc(['Sen2019_IR_R1', 'Sen2019_IR_R2'], ['Sen2019_Pro'])
fake  = lfc(['IMR90rep_Young_R1', 'IMR90rep_Young_R2'], ['Sen2019_Pro'])
fake2 = lfc(['BJrep_Young'], ['Sen2019_Pro'])
def sp(a, b, t=0.5):
    m = pd.DataFrame({'a': a, 'b': b}).dropna()
    m = m[(m.a.abs() >= t) | (m.b.abs() >= t)]
    return stats.spearmanr(m.a, m.b).statistic, len(m)
for lab, pair in [("I1 IR vs replicative (shared denom)", (ir, rep)),
                  ("NULL unrelated numerator vs replicative", (fake, rep)),
                  ("NULL 2 BJ Young vs replicative", (fake2, rep)),
                  ("NULL 3 two unrelated numerators", (fake, fake2))]:
    r, nn = sp(*pair)
    print(f"{lab:44s} rho {r:+.4f}  n={nn:,}")
