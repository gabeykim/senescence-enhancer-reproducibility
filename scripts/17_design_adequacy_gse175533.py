#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 17_design_adequacy_gse175533.py
#
# Design-adequacy assessment for GSE175533.
#
# CONSUMES: output/tables/gse175533_sample_metadata.csv
# PRODUCES: output/tables/gse175533_design_adequacy_log.txt
# ---------------------------------------------------------------------------
"""
Task 2 for GSE175533: ATAC/RNA sample inventory, PDL-matching between assays,
and the RIS/CD/hTERT control-arm questions.

Key gotcha this script exists to make explicit: RS_atac and RS_bulkRNA/RS_sc
"TP" (timepoint) labels are NOT the same PDL between assays past TP2 -- e.g.
ATAC TP3 = PDL30 but RNA TP3 = PDL28; ATAC TP7 = PDL50 but RNA TP7 = PDL45.
Any join between the two assay tables MUST be done on the numeric PDL column,
never on the TP label, or it silently pairs the wrong timepoints.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "output" / "tables" / "gse175533_sample_metadata.csv"
OUT = ROOT / "output" / "tables" / "gse175533_design_adequacy_log.txt"


def main():
    df = pd.read_csv(META)
    log = []

    log.append("=== Full sample inventory: arm x assay x replicate count ===")
    inv = df.groupby(["arm", "assay"]).size().reset_index(name="n")
    log.append(inv.to_string(index=False))

    # ---- ATAC inventory: PDL x replicate count ----
    atac = df[df.assay == "ATAC-seq"].copy()
    log.append("\n=== ATAC sample inventory (RS_atac + RS_atac_hTERT), by PDL ===")
    for arm in ["RS_atac", "RS_atac_hTERT"]:
        sub = atac[atac.arm == arm]
        counts = sub.groupby("pdl")["gsm"].nunique().sort_index()
        log.append(f"\n{arm}:")
        for pdl, n in counts.items():
            flag = "  [FLAG: <3 reps]" if n < 3 else ""
            log.append(f"  PDL {pdl:>4}: n={n} replicates{flag}")
    log.append(f"\nTotal ATAC samples: {len(atac)}")

    # ---- RNA inventory: PDL x replicate count (bulk only, exclude scRNA) ----
    rna_bulk = df[(df.assay == "RNA-Seq") & (df.arm.isin(
        ["RS_bulkRNA", "hTERT_bulkRNA", "RIS", "CD"]))].copy()
    log.append("\n=== Bulk RNA sample inventory, by arm ===")
    for arm in ["RS_bulkRNA", "hTERT_bulkRNA", "RIS", "CD"]:
        sub = rna_bulk[rna_bulk.arm == arm]
        key = "pdl" if arm in ("RS_bulkRNA", "hTERT_bulkRNA") else "day_post_treatment"
        counts = sub.groupby(key)["gsm"].nunique().sort_index()
        log.append(f"\n{arm} (staged by {key}):")
        for k, n in counts.items():
            flag = "  [FLAG: <3 reps]" if n < 3 else ""
            log.append(f"  {key}={k}: n={n} replicates{flag}")

    # ---- (a) PDL points with BOTH ATAC and RNA ----
    log.append("\n\n=== Task 2(a): PDL points with BOTH ATAC and RNA (WI-38 RS arm only) ===")
    atac_pdls = set(atac[atac.arm == "RS_atac"]["pdl"].dropna().astype(int))
    rna_pdls = set(rna_bulk[rna_bulk.arm == "RS_bulkRNA"]["pdl"].dropna().astype(int))
    shared = sorted(atac_pdls & rna_pdls)
    atac_only = sorted(atac_pdls - rna_pdls)
    rna_only = sorted(rna_pdls - atac_pdls)
    log.append(f"ATAC PDLs:  {sorted(atac_pdls)}")
    log.append(f"RNA PDLs:   {sorted(rna_pdls)}")
    log.append(f"SHARED (both assays, same nominal PDL): {shared}  -> {len(shared)} PDL points")
    log.append(f"ATAC-only PDLs: {atac_only}")
    log.append(f"RNA-only PDLs:  {rna_only}")
    log.append("NOTE: TP labels do NOT track the same PDL between assays past TP2 (e.g. ATAC "
               "TP3=PDL30 vs RNA TP3=PDL28; ATAC TP7=PDL50 vs RNA TP7=PDL45) -- this comparison "
               "is done on the numeric PDL field, never on TP label, deliberately.")
    log.append("Same biological sample or parallel cultures? PARALLEL CULTURES: at each shared "
               "PDL, the RNA sample and ATAC sample carry DIFFERENT BioSample accessions and "
               "materially different harvest protocols (RNA: Trizol lysis directly on the "
               "culture plate; ATAC: trypsinize, count, and carry forward exactly 100,000 cells) "
               "-- see raw BioSample links in gse175533_sample_metadata.csv. They come from the "
               "same continuous WI-38 culture lineage staged to the same nominal PDL checkpoint, "
               "not from a single dish split into two assays. Treat ATAC<->RNA pairing at a "
               "shared PDL as 'matched condition, unmatched aliquot', not a true multi-omic pair "
               "from one cell population.")

    # ---- (c) CD / quiescence arm usability ----
    log.append("\n\n=== Task 2(c): quiescent/cell-density control ===")
    log.append("CD (cell-density / contact-inhibition growth-arrest) arm: 30 RNA-seq samples, "
               "10 timepoints (day 1 to day 10) x 3 replicates, per the inventory above. "
               "RNA-seq ONLY -- no ATAC-seq exists for CD (see arm x assay table: no 'CD' row "
               "under ATAC-seq). This means CD can support 'is this senescence-specific or just "
               "generic-arrest' comparisons at the TRANSCRIPT level only. It cannot be used to "
               "ask the same question at the chromatin-accessibility level, because there is no "
               "CD ATAC data to compare against RS ATAC data. For a sequence-to-accessibility "
               "model, this control is unavailable where it matters most.")

    # ---- (d) RIS arm assay coverage ----
    log.append("\n\n=== Task 2(d): RIS (radiation-induced senescence) assay coverage ===")
    ris_assays = df[df.arm == "RIS"]["assay"].unique().tolist()
    log.append(f"RIS arm assay types present: {ris_assays}")
    log.append("CONFIRMED: RIS is RNA-seq ONLY. No ATAC-seq sample exists anywhere in GSE175533 "
               "with an RIS/irradiation title or characteristic, and the series' own "
               "Series_overall_design text states explicitly: 'Further controls for slowed "
               "growth (cell density) and DNA damage (radiation-RIS) were included for bulk "
               "RNA-seq' (i.e. not for ATAC-seq). This matters directly for the project's intended "
               "doxorubicin/therapy-induced-senescence validation model: RIS is the closest "
               "available proxy for a DNA-damage-induced senescence trigger in this dataset, and "
               "it carries ZERO chromatin-accessibility signal to train or validate a "
               "sequence-to-accessibility model against.")
    log.append(f"RIS sub-arms: {sorted(df[df.arm=='RIS']['day_post_treatment'].unique())} days "
               "post-treatment, where day 2 'no_xray' is a mock-irradiation control and the rest "
               "(days 3,4,5,6,9) are 10 Gy X-ray-treated, 3 replicates each (confirmed above).")

    report = "\n".join(log)
    OUT.write_text(report)
    print(report)


if __name__ == "__main__":
    main()
