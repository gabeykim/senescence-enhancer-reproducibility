#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 41_parse_expansion_metadata.py
#
# Assign samples to ON/OFF arms across all expansion studies.
#
# CONSUMES: data/*/ expression tables
# PRODUCES: output/tables/expansion_sample_assignments.csv, expansion_arm_rationale.json
# ---------------------------------------------------------------------------
"""
Arm assignment for the multi-study senescence expansion.

Parses SOFT metadata for each candidate expansion series, maps GEO sample records onto
the columns actually present in each deposited expression matrix, and emits one unified
sample-assignment table.

Every judgment call is recorded in ARM_RATIONALE (written verbatim into the AnnData .uns
by script 42) rather than being buried in the code.

Inputs : data/<ACC>/<ACC>_family.soft, deposited expression matrices
Outputs: output/tables/expansion_sample_assignments.csv
         output/tables/expansion_series_status.csv
"""
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
TAB = ROOT / "output" / "tables"
TAB.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------------------
# SOFT parsing (hand-rolled; no GEOparse, no HTML scraping -- same approach as script 02).
# --------------------------------------------------------------------------------------
def parse_soft_samples(path):
    """Return {gsm: {title, characteristics:[...], library_strategy, platform}}.

    Series-level fields are deliberately NOT propagated onto samples: GEO copies
    Sample_treatment_protocol_ch1 boilerplate identically onto every sample, so per-sample
    status must never be inferred from it (established in scripts/02_parse_metadata.py).
    """
    out, cur = {}, None
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = line.split("=", 1)[1].strip()
                out[cur] = {"gsm": cur, "title": "", "characteristics": [],
                            "library_strategy": "", "platform": ""}
            elif cur is None or not line.startswith("!Sample_"):
                continue
            elif line.startswith("!Sample_title"):
                out[cur]["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1"):
                out[cur]["characteristics"].append(line.split("=", 1)[1].strip())
            elif line.startswith("!Sample_library_strategy"):
                out[cur]["library_strategy"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_platform_id"):
                out[cur]["platform"] = line.split("=", 1)[1].strip()
    return out


def char(rec, key):
    for c in rec["characteristics"]:
        if c.lower().startswith(key.lower() + ":"):
            return c.split(":", 1)[1].strip()
    return ""


# --------------------------------------------------------------------------------------
# Recorded judgment calls -- these are the assumptions, stated, with their reasoning.
# --------------------------------------------------------------------------------------
ARM_RATIONALE = {
    "GSE99028": {
        "verdict": "INCLUDE -- as held-out transfer test, not training",
        "cell_line": "IMR90 (foetal lung fibroblast)",
        "mechanism": "therapy-induced senescence (TIS), etoposide",
        "ON": "SC-etoposide-rep1, SC-etoposide-rep2  (n=2)",
        "OFF": "PD29, PD34  (n=2)",
        "EXCLUDED": "shcGas-etoposide-rep1/2 -- cGAS knockdown suppresses the SASP arm of "
                    "the senescence programme, so these are a deliberately partial "
                    "phenotype, not a clean ON arm.",
        "assumption_column_mapping":
            "The deposited matrix columns are PD29, PD34, SC-etoposide-rep{1,2}, "
            "shcGas-etoposide-rep{1,2}; the SOFT titles are 'Proliferating control (PD29)', "
            "'Proliferating control (PD34)', 'sh-NTC etoposide-replicate {1,2}', "
            "'sh-cGAS etoposide-replicate {1,2}'. PD29/PD34 map by literal string. "
            "'SC' is read as 'scrambled control' == the SOFT's 'sh-NTC' (non-targeting "
            "control). Evidence: it is the only remaining pairing (the shcGas columns are "
            "named explicitly), GSM order in the SOFT matches column order left-to-right, "
            "and SOFT 'mutation:' is 'Wild-type' for the SC pair and 'cGAS shRNA knockdown' "
            "for the shcGas pair. ASSUMPTION -- GEO does not state the column-to-GSM map.",
        "caveat_vector_control":
            "The two proliferating controls carry SOFT 'treatment: None' and no shRNA "
            "designation, while the ON arm carries a non-targeting shRNA. The contrast "
            "therefore confounds etoposide with shRNA vector + selection. Small relative to "
            "the etoposide effect, but not zero, and it is not correctable from these data.",
        "why_holdout":
            "Etoposide and doxorubicin are both topoisomerase-II poisons, so this series is "
            "the closest available proxy for the project's eventual doxorubicin validation. "
            "It is a 2-vs-2 contrast: as training data it would add almost nothing (2 ON "
            "tasks) but as a test set it is the only way to measure whether a model trained "
            "on replicative + oncogene-induced senescence transfers to therapy-induced "
            "senescence. Holding it out keeps one induction mechanism entirely unseen, which "
            "is the only honest test of the shared-regulatory-logic assumption.",
    },
    "GSE74324": {
        "verdict": "INCLUDE in training",
        "cell_line": "IMR90 (foetal lung fibroblast)",
        "mechanism": "oncogene-induced senescence (OIS), HRAS-V12",
        "ON": "H_Ren_1..3  (n=3)",
        "OFF": "P_Ren_1..3  (n=3)",
        "QUIESCENT": "PQ_Ren_1..3  (n=3) -- third class, matching the existing build's "
                     "contact-inhibition/quiescent handling",
        "EXCLUDED": "63 of 72 samples. The series is an 8-arm perturbation panel "
                    "(shBrd4, shp65, shp53, shp53Rb, shp16p21, DMSO vehicle, JQ1 drug) "
                    "x 3 conditions x 3 replicates. Only the shRenilla no-treatment arm "
                    "(_Ren) is an unperturbed P/Q/S baseline. The _DMSO arm is the vehicle "
                    "control of the JQ1 sub-experiment -- a different protocol branch, not "
                    "a second replicate of _Ren -- so it is not pooled with it.",
        "assumption_prefix_decoding":
            "Column prefixes decoded as P=proliferating, PQ=quiescent (contact-inhibited), "
            "H=HRAS-V12 senescent, consistent with the series' P/Q/S triad design. "
            "'Ren' = shRenilla, the standard non-targeting control shRNA in this lab's "
            "system. ASSUMPTION -- decoded from the design description, not from an "
            "explicit GEO column legend.",
    },
    "GSE205692": {
        "verdict": "INCLUDE in training, with the cross-platform caveat below",
        "cell_line": "GM21 (adult skin fibroblast)",
        "mechanism": "oncogene-induced senescence (OIS), H-RAS-G12V retroviral",
        "ON": "RAS arm, days 13-25 only (see why_day_window)",
        "OFF": "empty vector, all days (8/32/56), 3 biological replicates",
        "EXCLUDED": "RAS days 8, 32, 45, 46, 50, 56 -- see why_day_window.",
        "why_day_window":
            "This SuperSeries' entire subject is ESCAPE from OIS: a subpopulation resumes "
            "proliferation at late timepoints. Late RAS samples are therefore NOT reliably "
            "senescent and labelling them ON would inject mislabelled positives -- exactly "
            "the poisoning failure mode the task brief warns about. Confirmed by "
            "measurement (script 42 step 5b): days 32-56 show CXCL8 +0.60 vs +5.43 in the "
            "chosen window and MKI67 -0.53 vs -0.83, i.e. a markedly attenuated senescence "
            "signature, consistent with a partially escaped population.",
        "why_day8_excluded_CORRECTED":
            "Day 8 was excluded a priori on the reasoning that arrest was 'not yet "
            "established'. THAT REASONING WAS WRONG and is corrected here rather than left "
            "standing. Measurement shows day 8 is the MOST arrested timepoint in the series "
            "(LMNB1 -0.71, MKI67 -0.77 vs empty vector; arrest z-score +1.54, higher than "
            "any ON day). What day 8 lacks is the SASP: IL6 is -1.05, i.e. BELOW the "
            "empty-vector control. The exclusion therefore stands, but on the opposite "
            "ground: day 8 is arrest-without-SASP -- early, incomplete senescence -- and it "
            "fails this project's pre-registered 'CDKN1A/IL6/CXCL8 up' criterion outright. "
            "Including it would put a phenotypically different arm under the same ON label "
            "as the three fully-SASP-positive ON arms in the other studies.",
        "caveat_weak_CDKN1A":
            "GSE205692 clears the CDKN1A-up criterion only marginally (+0.146 log2FC, vs "
            "+0.52 / +1.50 / +2.53 in the three RNA-seq studies). Its marker evidence rests "
            "mainly on CXCL8 (+5.43), IL6 (+1.02) and MKI67 (-0.83). Part of this is the "
            "microarray's compressed dynamic range and the submitters' sva correction, but "
            "it means this study contributes a weaker-magnitude senescence contrast than "
            "its sample count suggests.",
        "caveat_platform":
            "This is an Affymetrix HTA 2.0 microarray (GPL17586), not RNA-seq. Values are "
            "log2 RMA intensities, already sva batch-corrected by the submitters. Pooling "
            "microarray with RNA-seq is the most severe cross-study harmonisation problem "
            "in this build and is handled explicitly by the per-study correction in "
            "script 42, with the PCA reported before and after.",
        "why_this_subseries":
            "GSE206402 (ATAC) is the accession named in the brief, but ATAC peaks cannot "
            "supply per-gene expression targets. Its named matched RNA-seq, GSE206493, "
            "turns out to contain NO proliferating control (all 12 samples are "
            "RAS-overexpressing, differing only by shPOU2F2 vs pLKO). The empty-vector vs "
            "RAS contrast in this study exists ONLY in GSE205692, the undisclosed 4th "
            "subseries found in the prior audit.",
    },
    "GSE206493": {
        "verdict": "EXCLUDE -- no control arm exists",
        "reason":
            "All 12 samples are H-RAS-G12V-overexpressing GM21 at day 14 or day 27, "
            "differing only by shRNA (pLKO non-targeting vs shPOU2F2 #324 vs #325). There "
            "is no empty-vector, no proliferating, and no untreated arm anywhere in the "
            "series, so no ON-vs-OFF senescence contrast can be formed from it. The single "
            "deposited processed file is GSE206493_Counts_POU2F2_KD.txt.gz, whose name "
            "states the same thing. This contradicts the task brief's description of "
            "GSE206493 as the matched RNA-seq carrying empty-vector proliferating controls.",
    },
    "GSE210285": {
        "verdict": "EXCLUDE -- no expression data exists",
        "reason":
            "The series contains 4 samples, all ATAC-seq (2 Growing, 2 Senescence, 2BS "
            "cells). Despite the series summary describing RNA-seq and ChIP-seq, none was "
            "deposited under this accession, and a GEO search for a sibling 2BS senescence "
            "RNA-seq series returned no such record. The only supplementary data are bigWig "
            "signal tracks (no peak files either, as the prior audit found). With no "
            "expression matrix there is nothing to contribute to an expression-target "
            "training set.",
    },
}


def main():
    rows, status = [], []

    # ---------------- GSE99028 : IMR90, etoposide TIS ----------------
    s = parse_soft_samples(ROOT / "data/GSE99028/GSE99028_family.soft")
    mat = pd.read_csv(ROOT / "data/GSE99028/GSE99028_read_count.txt", sep="\t", nrows=1)
    cols = [c for c in mat.columns if c != "refGene"]
    order = sorted(s)                                   # GSM order == column order (assumption)
    col_map = dict(zip(order, cols))
    for gsm in order:
        rec, c = s[gsm], col_map[gsm]
        t = rec["title"]
        if "Proliferating control" in t:
            cls, role = "proliferating", "test_transfer"
        elif "sh-NTC" in t:
            cls, role = "senescent", "test_transfer"
        else:
            cls, role = "senescent_perturbed", "excluded"
        rows.append(dict(study="GSE99028", gsm=gsm, matrix_column=c, title=t,
                         cell_line="IMR90", mechanism="etoposide_TIS",
                         platform_kind="RNA-seq", quantification="raw counts",
                         cls=cls, role=role,
                         detail=char(rec, "treatment") + " | " + char(rec, "mutation")))
    status.append(dict(study="GSE99028", decision="include_as_transfer_test",
                       cell_line="IMR90", mechanism="etoposide_TIS",
                       platform_kind="RNA-seq", quantification="raw counts",
                       n_on=2, n_off=2, n_excluded=2))

    # ---------------- GSE74324 : IMR90, OIS (clean shRen baseline only) ----------------
    hdr = pd.read_csv(ROOT / "data/GSE74324/GSE74324_ALL_samples_rpkm.txt",
                      sep="\t", nrows=0)
    sample_cols = [c for c in hdr.columns if c not in
                   ("", "entrezgene", "ensembl_gene_id", "external_gene_name",
                    "chromosome_name", "band", "start_position", "end_position",
                    "transcript_length")]
    PREFIX = {"P": ("proliferating", "train_eligible"),
              "PQ": ("quiescent", "train_eligible"),
              "H": ("senescent", "train_eligible")}
    for c in sample_cols:
        m = re.match(r"^(PQ|P|H)_(.+?)_(\d)$", c)
        if not m:
            continue
        pref, background, rep = m.groups()
        clean = background == "Ren"
        cls, role = PREFIX[pref]
        rows.append(dict(study="GSE74324", gsm="", matrix_column=c, title=c,
                         cell_line="IMR90", mechanism="HRASV12_OIS",
                         platform_kind="RNA-seq", quantification="RPKM",
                         cls=cls if clean else cls + "_perturbed",
                         role=role if clean else "excluded",
                         detail=f"background={background} rep={rep}"))
    status.append(dict(study="GSE74324", decision="include_in_training",
                       cell_line="IMR90", mechanism="HRASV12_OIS",
                       platform_kind="RNA-seq", quantification="RPKM",
                       n_on=3, n_off=3, n_excluded=63))

    # ---------------- GSE205692 : GM21, OIS microarray ----------------
    s = parse_soft_samples(ROOT / "data/GSE205692/GSE205692_family.soft")
    ON_DAYS = {13, 14, 18, 19, 21, 22, 23, 25}
    for gsm, rec in sorted(s.items()):
        t = rec["title"]
        day = int(re.search(r"day (\d+)", t).group(1))
        rep = int(re.search(r"biol rep (\d+)", t).group(1))
        is_ras = "RAS" in t
        if not is_ras:
            cls, role = "proliferating", "train_eligible"
        elif day in ON_DAYS:
            cls, role = "senescent", "train_eligible"
        else:
            cls = "ras_pre_senescent" if day < min(ON_DAYS) else "ras_post_escape"
            role = "excluded"
        rows.append(dict(study="GSE205692", gsm=gsm, matrix_column=gsm, title=t,
                         cell_line="GM21", mechanism="HRASG12V_OIS",
                         platform_kind="microarray_HTA2.0",
                         quantification="log2 RMA (sva-corrected)",
                         cls=cls, role=role, detail=f"day={day} biol_rep={rep}"))
    n_on = sum(1 for r in rows if r["study"] == "GSE205692" and r["role"] == "train_eligible"
               and r["cls"] == "senescent")
    status.append(dict(study="GSE205692", decision="include_in_training",
                       cell_line="GM21", mechanism="HRASG12V_OIS",
                       platform_kind="microarray_HTA2.0",
                       quantification="log2 RMA (sva-corrected)",
                       n_on=n_on, n_off=9, n_excluded=34 - n_on - 9))

    # ---------------- excluded series, recorded so the report can state why ----------------
    status.append(dict(study="GSE206493", decision="EXCLUDE_no_control_arm",
                       cell_line="GM21", mechanism="HRASG12V_OIS",
                       platform_kind="RNA-seq", quantification="raw counts",
                       n_on=0, n_off=0, n_excluded=12))
    status.append(dict(study="GSE210285", decision="EXCLUDE_no_expression_data",
                       cell_line="2BS", mechanism="replicative",
                       platform_kind="ATAC-seq only", quantification="none",
                       n_on=0, n_off=0, n_excluded=4))

    df = pd.DataFrame(rows)
    df.to_csv(TAB / "expansion_sample_assignments.csv", index=False)
    st = pd.DataFrame(status)
    st.to_csv(TAB / "expansion_series_status.csv", index=False)
    (TAB / "expansion_arm_rationale.json").write_text(json.dumps(ARM_RATIONALE, indent=2))

    print(df.groupby(["study", "cls", "role"]).size().to_string())
    print()
    print(st.to_string(index=False))
    print(f"\nwrote {TAB/'expansion_sample_assignments.csv'} ({len(df)} rows)")


if __name__ == "__main__":
    main()
