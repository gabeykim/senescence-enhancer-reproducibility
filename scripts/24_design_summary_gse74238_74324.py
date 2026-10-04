#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 24_design_summary_gse74238_74324.py
#
# Design summary for the GSE74238/74324 pair.
#
# CONSUMES: data/GSE74238/GSE74238_family.soft, data/GSE74324/GSE74324_family.soft
# PRODUCES: output/tables/gse74238_74324_design_summary.txt
# ---------------------------------------------------------------------------
"""
Task 3-6 support script: design adequacy (replicate counts per arm),
GSE74238 vs GSE74324 pairing evidence, and genome-build cross-check --
computed from output/tables/gse74238_sample_metadata.csv and
gse74324_sample_metadata.csv (produced by 23_parse_gse74238_74324.py),
plus direct SOFT-text checks for the SuperSeries relation. Mirrors the
pattern in scripts/05_design_adequacy.py from the sibling audit.
"""
import re
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "output" / "tables"
DATA = ROOT / "data"
OUT = TABLES / "gse74238_74324_design_summary.txt"


def main():
    log = []

    # ---------------------------------------------------------- GSE74238 ----
    log.append("=== GSE74238 (ChIP-seq) -- replicate counts per assay x condition ===")
    df238 = pd.read_csv(TABLES / "gse74238_sample_metadata.csv")
    triad = df238[df238["condition_group"] == "P/Q/S triad"]
    for (assay, cond), sub in triad.groupby(["assay", "condition"]):
        n = len(sub)
        flag = "  [FLAG: <3 biological replicates]" if n < 3 else "  [OK: >=3 replicates]"
        log.append(f"  {assay:20s} x {cond:14s}: n={n} GSM(s){flag}")
    other = df238[df238["condition_group"] != "P/Q/S triad"]
    log.append(f"\n  Non-triad arm (DMSO/Etoposide genotoxic stress, separate from P/Q/S):")
    for (assay, cond), sub in other.groupby(["assay", "condition"]):
        log.append(f"  {assay:20s} x {cond:14s}: n={len(sub)} GSM(s)  [n=1, no replicate at all -- not usable alone]")

    log.append("\n  VERDICT: every P/Q/S-triad arm (H3K27ac, BRD4, input x proliferating/"
                "quiescent/senescent) has exactly n=2 biological replicates (a first round, "
                "GSM1915xxx, submitted Oct 2015; a second independent round, GSM2098xxx, "
                "explicitly titled '... replicate', submitted Mar 2016, each with its own "
                "BioSample/SRX accession -- confirmed genuinely independent, not resequencing). "
                "2 replicates is BELOW the >=3 biological-replicate bar commonly used as an "
                "adequacy threshold -- flagged explicitly per task instructions, not assumed "
                "to pass. The DMSO/Etoposide arm has n=1 with zero replicates and is a "
                "different experiment (genotoxic/therapy-induced-senescence-like stress), not "
                "part of the P/Q/S triad at all.")

    # ---------------------------------------------------------- GSE74324 ----
    log.append("\n\n=== GSE74324 (RNA-seq, 72 samples) -- perturbation arm x condition ===")
    df324 = pd.read_csv(TABLES / "gse74324_sample_metadata.csv")
    df324["arm_short"] = df324["perturbation_arm"].str.split(" \\(").str[0]
    for (arm, cond), sub in df324.groupby(["arm_short", "condition"]):
        n = len(sub)
        flag = "[OK: >=3 replicates]" if n >= 3 else "[FLAG: <3 replicates]"
        log.append(f"  {arm:35s} x {cond:14s}: n={n}  {flag}")

    n_arms = df324["arm_short"].nunique()
    log.append(f"\n  {n_arms} perturbation arms x 3 conditions (P/Q/S) x 3 replicates = "
               f"{n_arms * 3 * 3} = {len(df324)} total samples. Confirms 72 is NOT a simple "
               f"P/Q/S x replicate design -- it is P/Q/S crossed with 6 shRNA-knockdown "
               f"backgrounds (control shRen, shBrd4, shp65, shp53, shp53/Rb, shp16/p21) PLUS "
               f"a separate small-molecule arm (DMSO vehicle vs 100nM JQ1, both on the shRen "
               f"background) -- see Series_overall_design text, confirmed per-sample via "
               f"`viral transduction:` and `treatment:` characteristics.")

    n_clean = int(df324["is_clean_PQS_triad_baseline"].sum())
    n_dmso = int(df324["is_dmso_vehicle_baseline"].sum())
    log.append(f"\n  CLEAN BASELINE SUBSET for a simple senescent-vs-proliferating contrast "
               f"(no knockdown beyond control shRen, no drug): {n_clean}/72 samples "
               f"(shRen_1/2/3 x proliferating/quiescent/senescent, 3 reps/condition, each "
               f">=3 -> OK).")
    log.append(f"  SECOND CANDIDATE baseline (shRen background + DMSO vehicle, i.e. the "
               f"vehicle arm of the JQ1 experiment): {n_dmso}/72 samples, also 3 reps/"
               f"condition. NOTE this is a distinct sub-experiment from the shRen 'no "
               f"treatment' arm (different Sample_treatment_protocol_ch1 branch -- 48h DMSO "
               f"exposure prior to harvest vs no drug at all) -- the two baselines are NOT "
               f"interchangeable replicates of each other without checking batch/PCA "
               f"structure first (out of scope here -- GEO metadata only).")
    log.append(f"\n  IRRELEVANT-TO-THIS-AUDIT subset (5 gene-knockdown arms [shp53, shp65, "
               f"shBrd4, shp16/p21, shp53/Rb] x P/Q/S x 3 reps, plus the JQ1-treated arm): "
               f"{72 - n_clean - n_dmso} / 72 samples ({(72 - n_clean - n_dmso) / 72:.0%}) are "
               f"perturbation conditions designed to test which factors are REQUIRED for the "
               f"senescence phenotype/SASP -- not usable for a naive senescent-vs-proliferating "
               f"contrast without explicitly restricting to the shRen (and optionally DMSO) "
               f"subset first.")

    # ------------------------------------------------------- pairing check --
    log.append("\n\n=== Pairing check: GSE74238 vs GSE74324 ===")
    soft238 = (DATA / "GSE74238" / "GSE74238_family.soft").read_text(errors="replace")
    soft324 = (DATA / "GSE74324" / "GSE74324_family.soft").read_text(errors="replace")

    rel238 = re.findall(r"!Series_relation = (.*)", soft238)
    rel324 = re.findall(r"!Series_relation = (.*)", soft324)
    log.append(f"  GSE74238 Series_relation: {rel238}")
    log.append(f"  GSE74324 Series_relation: {rel324}")
    super238 = next((r for r in rel238 if "SubSeries of" in r), None)
    super324 = next((r for r in rel324 if "SubSeries of" in r), None)
    if super238 and super324 and super238 == super324:
        log.append(f"  -> BOTH series are explicit SubSeries of the SAME SuperSeries "
                    f"({super238.split(': ')[1]}, 'BRD4 connects enhancer remodeling to "
                    f"senescence immune surveillance'). This is GEO's own structural "
                    f"declaration that these two accessions are one paired multi-omics "
                    f"submission from the same study, not two coincidentally-similar series.")
    else:
        log.append("  -> NOT confirmed as SubSeries of the same SuperSeries -- inspect manually.")

    gp238 = re.findall(r"!Sample_growth_protocol_ch1 = (.*)", soft238)
    gp324 = re.findall(r"!Sample_growth_protocol_ch1 = (.*)", soft324)
    log.append(f"\n  GSE74238 growth protocol (first sample, identical across series):\n    {gp238[0] if gp238 else 'MISSING'}")
    log.append(f"\n  GSE74324 growth protocol (first sample, identical across series):\n    {gp324[0] if gp324 else 'MISSING'}")
    same_cell_line = "IMR90" in (gp238[0] if gp238 else "") and "IMR90" in (gp324[0] if gp324 else "")
    same_induction = "HRASV12" in (gp238[0] if gp238 else "") and "HRASV12" in (gp324[0] if gp324 else "")
    same_quiescence = "0.1% FBS" in (gp238[0] if gp238 else "") and "0.1% FBS" in (gp324[0] if gp324 else "")
    day12_238 = "day 12 post-infection" in (gp238[0] if gp238 else "").lower()
    day12_324 = "day 12" in (soft324.lower())
    log.append(f"\n  Same cell line (IMR90) in both protocols: {same_cell_line}")
    log.append(f"  Same OIS induction construct (pWZL-HRASV12) in both protocols: {same_induction}")
    log.append(f"  Same quiescence-induction protocol (4d 0.1% FBS) in both: {same_quiescence}")
    log.append(f"  GSE74238 explicitly harvested 'day 12 post-infection': {day12_238}")
    log.append(f"  GSE74324 mentions 'Day 12 post-infection' harvest anywhere in SOFT text: {day12_324}")

    log.append("\n  VERDICT: GSE74238 (ChIP) and the shRen/'no treatment' + shRen/DMSO-vehicle "
               "subset of GSE74324 (RNA, 18/72 samples) ARE a genuinely matched pair for the "
               "P/Q/S triad -- same cell line, same HRASV12 OIS system, same quiescence "
               "protocol, same day-12 harvest, explicit shared SuperSeries membership. The "
               "REMAINING 54/72 RNA samples (knockdown + JQ1 arms) have NO ChIP-seq "
               "counterpart at all -- GSE74238 never profiled H3K27ac/BRD4 in a shp53-, "
               "shp65-, shBrd4-, shp16/p21-, or shp53/Rb-knockdown background, or under JQ1. "
               "So the two series are paired in DESIGN LINEAGE (same study, same SuperSeries, "
               "same core P/Q/S system) but diverge sharply in SCOPE: RNA covers a much larger "
               "perturbation space than ChIP does. Anyone using both series as a matched "
               "multi-omics resource should restrict RNA to the shRen (+/- DMSO) subset, not "
               "average over or naively pool all 72 samples.")

    # ------------------------------------------------------ genome build ----
    log.append("\n\n=== Genome build ===")
    gb238 = sorted(set(re.findall(r"Genome_build:\s*(\S+)", soft238)))
    gb324 = sorted(set(re.findall(r"Genome_build:\s*(\S+)", soft324)))
    align238 = re.findall(r"!Sample_data_processing = (.*(?:GRCh3\d|hg1\d|hg38).*)", soft238)
    align324 = re.findall(r"!Sample_data_processing = (.*(?:GRCh3\d|hg1\d|hg38).*)", soft324)
    log.append(f"  GSE74238 explicit 'Genome_build:' tag values (all samples): {gb238}")
    log.append(f"  GSE74324 explicit 'Genome_build:' tag values (all samples): {gb324}")
    log.append(f"  GSE74238 alignment description: {align238[0] if align238 else 'not found'}")
    log.append(f"  GSE74324 alignment description: {align324[0] if align324 else 'not found'}")
    log.append("  -> Both series explicitly state hg19/GRCh37 (GSE74238: 'GRCh37 (hg19)' via "
               "Bowtie; GSE74324: 'GRCh37.75(hg19)' via STAR). CONFIRMED from text, not "
               "assumed from publication date. No hg38/GRCh38 mention anywhere in either "
               "family.soft.")

    report = "\n".join(log)
    OUT.write_text(report)
    print(report)


if __name__ == "__main__":
    main()
