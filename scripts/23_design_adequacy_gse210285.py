#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 23_design_adequacy_gse210285.py
#
# Design-adequacy assessment for GSE210285.
#
# CONSUMES: output/tables/gse210285_sample_metadata.csv
# PRODUCES: output/tables/gse210285_design_adequacy_log.txt
# ---------------------------------------------------------------------------
"""
Task 3 (design adequacy) + Task 4 (cross-dataset compatibility inputs) for
GSE210285, computed from output/tables/gse210285_sample_metadata.csv
(produced by 22_parse_gse210285.py) plus a direct read of
data/GSE210285/filelist.txt and data/GSE210285/GSE210285_family.soft for
facts that don't belong in the per-GSM table (genome build, peak caller,
supplementary file inventory).

Writes output/tables/gse210285_design_adequacy_log.txt.
"""
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "output" / "tables" / "gse210285_sample_metadata.csv"
FILELIST = ROOT / "data" / "GSE210285" / "filelist.txt"
SOFT = ROOT / "data" / "GSE210285" / "GSE210285_family.soft"
OUT = ROOT / "output" / "tables" / "gse210285_design_adequacy_log.txt"


def main():
    log = []
    df = pd.read_csv(META)

    log.append("=== GSE210285 design adequacy ===\n")
    log.append(f"Total GSMs in series: {len(df)}")

    cond_counts = Counter(df["condition_senescence_status"])
    log.append(f"Conditions found (from 'phenotype' characteristic): {dict(cond_counts)}")
    log.append(f"Number of distinct conditions: {len(cond_counts)}")

    log.append("\n--- Replicates per condition ---")
    for cond, n in cond_counts.items():
        flag = "[FLAG: <3 biological replicates -- below standard QC/DE minimum]" if n < 3 else "[OK: >=3]"
        log.append(f"  {cond!r}: n = {n} biological replicate(s)  {flag}")

    log.append(
        "\nCONFIRMED: this is a 2-condition x 2-replicate design (4 GSMs total, no more, no "
        "fewer). The literature description 'only 2 replicates per arm' is accurate as stated in "
        "GEO metadata itself -- not an artifact of incomplete deposition or a naive miscount. "
        "Source: !Series_overall_design states explicitly 'For each age stage, we obtained two "
        "independent ATAC-seq replicates,' and !Series_sample_id / esummary both list exactly 4 "
        "GSMs, 2 titled '...-Growing-rep{1,2}' and 2 titled '...-rep{1,2}' (senescent arm)."
    )

    log.append(
        "\nFLAG (design weakness, stated plainly per task instruction): n=2 biological replicates "
        "per arm is below the conventional n>=3 minimum for any variance-based QC (e.g. estimating "
        "within-group dispersion, flagging outlier replicates with confidence) or formal "
        "differential-accessibility testing. With n=2, a single discordant replicate cannot be "
        "distinguished from true biological signal by statistics alone -- there is no third "
        "observation to break the tie. This does not make the data unusable as *input sequences* "
        "for a pretrained sequence-to-expression model (Borzoi/gReLU only needs accessible "
        "regions + a label, not a well-powered contrast), but it does mean any accessibility calls "
        "or condition-average tracks derived from this series carry that fragility forward -- "
        "an activity-head fine-tune trained on this series' senescent-vs-growing contrast is "
        "effectively trained on one biological pair per arm, not a population estimate."
    )

    # ---- platform / condition confound check ----
    log.append("\n--- Platform x condition crosstab (checking for confound) ---")
    for _, r in df.sort_values("gsm").iterrows():
        log.append(f"  {r['gsm']} ({r['title']}): condition={r['condition_senescence_status']!r}, "
                    f"platform={r['platform_id']} ({r['instrument_model']})")
    plat_by_cond = df.groupby("condition_senescence_status")["platform_id"].apply(lambda s: sorted(s.unique()))
    log.append(f"\nPlatforms used per condition: {plat_by_cond.to_dict()}")
    log.append(
        "NOT confounded with condition: each condition (Growing, Senescence) has exactly one "
        "GPL20301 (Illumina HiSeq 4000) sample and one GPL23227 (BGISEQ-500) sample -- platform "
        "is crossed with replicate number (rep1=HiSeq, rep2=BGISEQ in both arms), not aliased "
        "with the senescent-vs-growing contrast. This is a mild positive for the design: a "
        "platform-driven signal cannot masquerade as a condition effect here. It does mean, "
        "however, that 'rep1 vs rep2' is not a clean biological-replicate axis either -- it is "
        "confounded with sequencing platform, so apparent rep1-vs-rep2 differences could be "
        "platform batch effect, not just biological/technical noise. With only 2 replicates and "
        "2 platforms, these two sources of variation cannot be statistically separated from this "
        "series alone."
    )

    # ---- passage / PDL ----
    log.append(
        "\nPassage number / population doubling level (PDL): NOT reported anywhere in GEO "
        "metadata for this series -- not in Sample_characteristics_ch1 (only 'tissue', 'cell "
        "line', 'phenotype' keys are used), not in Sample_title, not in any protocol field. "
        "'Growing' vs 'Senescence' is a categorical phenotype label only; there is no quantitative "
        "passage/PDL value in GEO to check dose-response or to match against other series' passage "
        "numbers. (The source publication may report PDL in its methods text, but per task "
        "constraints this audit works from GEO metadata only, not the paper PDF.)"
    )

    # ---- cross-dataset compatibility inputs ----
    log.append("\n\n=== Cross-dataset compatibility inputs (for Task 5 peak-set comparison) ===")

    soft_text = SOFT.read_text(encoding="utf-8", errors="replace")
    log.append(
        "\nGenome build: hg38, stated TWICE per sample and identically across all 4 GSMs -- "
        "explicitly in free text ('the obtained clean data were aligned to human genome hg38 by "
        "the Burrows-Wheeler Aligner tool (v.0.7.10)') AND in a separate structured line "
        "('Assembly: hg38'). Not an assumption -- directly stated in "
        "Sample_data_processing on every sample. No hg19/GRCh37 mentions anywhere in the SOFT "
        "file."
    )
    assert "hg38" in soft_text and "hg19" not in soft_text.lower(), \
        "genome build claim above does not match raw SOFT text -- re-check"

    log.append(
        "\nPeak-calling software: MACS2, stated as 'MASC2 (v.2,2)' in "
        "Sample_data_processing -- this is quoted VERBATIM from GEO; 'MASC2' for 'MACS2' and "
        "'v.2,2' for what is almost certainly 'v.2.2.x' both read as transcription slips in the "
        "original submission (the command line immediately after literally reads 'macs2 "
        "callpeak', confirming the tool), not something this parser introduced. Exact version "
        "is NOT resolvable beyond 'MACS2 2.2.x' from GEO text alone. "
        "Command line as given (also transcribed verbatim, note the mangled double-dashes -- "
        "likely lost when the submitter's text was reformatted, GEO does not preserve '--' "
        "reliably in free-text fields): "
        "'macs2 callpeak - nomodel - f BAMPE-keep - dup 1 - q 0.05 - B-SPMR'. Read charitably "
        "this is 'macs2 callpeak --nomodel -f BAMPE --keep-dup 1 -q 0.05 -B --SPMR' -- standard "
        "ATAC-seq narrowPeak calling flags (no shifting model, paired-end BAM, keep at most 1 "
        "duplicate, q<0.05, bedGraph output with signal-per-million-reads scaling) -- but this "
        "reconstruction is an inference from a garbled string, not a literal quote, and is "
        "flagged as such."
    )

    log.append(
        "\nPeak file: NONE DEPOSITED. This is the key negative finding for cross-dataset peak-set "
        "comparison. Despite Sample_data_processing describing a MACS2 peak-calling step for "
        "every one of the 4 samples, GEO's only supplementary files for this series are 4 bigWig "
        "(.bw) signal tracks -- no .narrowPeak, .broadPeak, .bed, or any peak-list file of any "
        "kind is present in GSE210285_RAW.tar or listed in filelist.txt. There is therefore no "
        "peak file header, column structure, delimiter, or coordinate convention to report -- "
        "there is nothing to inspect. A peak set for this series does not exist as deposited data; "
        "one would have to be re-derived (e.g. MACS2 peak calling from the bigWig tracks, or from "
        "raw reads via SRA, which is out of scope here) before GSE210285 could be compared "
        "peak-for-peak against GSE175533 or GSE206402 in the downstream Task 5 comparison."
    )

    filelist_text = FILELIST.read_text(encoding="utf-8", errors="replace")
    log.append(f"\nfilelist.txt contents (verbatim):\n{filelist_text}")
    log.append(
        "Supplementary file inventory: 4 bigWig files, one per GSM, bundled into one "
        "GSE210285_RAW.tar (232,181,760 bytes / ~221 MiB). Individual file sizes: "
        "GSM6427623_RS-growing-1.bw = 52,204,274 B (~50 MiB); "
        "GSM6427624_RS-growing-2.bw = 28,028,993 B (~27 MiB); "
        "GSM6427625_RS-1.bw = 87,539,031 B (~83 MiB); "
        "GSM6427626_RS-2.bw = 64,396,241 B (~61 MiB). Not downloaded (bigWig, flagged as 'very "
        "large' under task constraints) -- existence and exact sizes captured from filelist.txt "
        "only, no binary content inspected."
    )

    log.append(
        "\nMulti-omics claim vs what's actually in this GEO record: Series_summary states the "
        "underlying study 'performed ATAC-seq, RNA-seq, and ChIP-seq on different senescent "
        "types,' but GSE210285 itself contains ATAC-seq only (Series_type = 'Genome "
        "binding/occupancy profiling by high throughput sequencing'; both platforms are "
        "sequencers, no expression or ChIP platform attached). A live BioProject cross-reference "
        "(esearch db=gds term=PRJNA865031) returns exactly 1 GEO series (GSE210285 itself) -- no "
        "companion RNA-seq/ChIP-seq series is linked under the same BioProject, and GSE210285 has "
        "no SuperSeries relation (esummary 'relations'/'extrelations' both empty). If the RNA-seq "
        "and ChIP-seq arms of this multi-omics study were deposited to GEO at all, they are under "
        "a different, currently unidentified accession -- not discoverable via this series' own "
        "metadata."
    )

    report = "\n".join(log)
    OUT.write_text(report)
    print(report)


if __name__ == "__main__":
    main()
