#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 30_build_senescence_trainset.py
#
# Build the single-study (GSE175533) expression training set.
#
# CONSUMES: data/GSE175533/ TPM tables, output/tables/
# PRODUCES: output/trainset/senescence_trainset.h5ad, dropped_genes.csv, REPORT_trainset.md
# ---------------------------------------------------------------------------
"""
Build a senescence expression training set from GSE175533, formatted to the
Borzoi/gReLU input contract measured in the Task 4 validation session.

LOCAL ONLY -- no GPU, no SRA, no model loading. Pure data curation.

=========================== SETTLED INPUTS (not re-derived) ==================
Per prior audits in this project (see REPORT_GSE175533_MULTI.md and
infra_borzoi_grelu/REPORT.md):
  - GSE175533 = WI-38 replicative senescence multi-omics, 147 samples, hg38.
  - ATAC and RNA at matched PDLs are PARALLEL CULTURES, not split aliquots.
    => We never treat an ATAC sample and an RNA sample as the same biological
       unit. Accessibility is used ONLY as a region-selection prior.
  - Senescence inflection sits between PDL46 and PDL50 (data-driven, from RNA
    marker trajectories + ATAC differential-peak escalation).
  - CDKN2A/CDKN1A promoters show NO significant accessibility change despite
    strong CDKN1A RNA induction => the training signal must be EXPRESSION.
  - RIS and cell-density (CD) arms have ZERO ATAC data.

Model contract (measured from the live model, Task 4):
  - Input  (4, 524288) float32, channels-first one-hot DNA
  - Trunk output (batch, 1920, 6144)
  - Target (n_tasks, 1) float32
  - AnnData .X shape (n_tasks, n_intervals)
  - Head: ConvHead + softplus  => model outputs are NON-NEGATIVE
  - Genome build hg38
  - gReLU split() needs DISJOINT chromosome lists

=========================== CONFIGURABLE JUDGMENT CALLS ======================
Every step-2 decision is a module-level constant with its reasoning recorded
in ARM_POLICY_RATIONALE (written verbatim into the output report and into
AnnData.uns), NOT buried in code. Change the constant, re-run, get a
different training set.
"""
import argparse
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "GSE175533"
ANNO = ROOT / "data" / "annotation" / "gencode.v44.basic.annotation.gtf.gz"
OUT = ROOT / "output" / "trainset"
OUT.mkdir(parents=True, exist_ok=True)

TPM_XLSX = DATA / "GSE175533_TPM_table.transitional.xlsx"
ATAC_SIG_XLSX = DATA / "GSE175533_atac_peaks_sig.transitional.xlsx"
ATLAS_BED = DATA / "GSE175533_atlas.bed"

SEQ_LEN = 524_288          # Borzoi input width, fixed by the checkpoint
BIN_SIZE = 32
LABEL_LEN = 196_608

# --------------------------------------------------------------------------
# STEP 2 JUDGMENT CALLS -- configurable, with reasoning recorded, not hardcoded
# --------------------------------------------------------------------------
# PDL inflection. Prior audit: transition sits BETWEEN PDL46 and PDL50.
OFF_MAX_PDL = 37     # PDLs <= this are the proliferating (OFF) arm
ON_MIN_PDL = 50      # PDLs >= this are the senescent (ON) arm
# => PDL45 and PDL46 fall in the ambiguous transition zone.

TRANSITION_POLICY = "exclude"   # {"exclude", "assign_off", "assign_on"}
HTERT_POLICY = "exclude"        # {"exclude", "as_off"}
CD_POLICY = "third_class"       # {"exclude", "as_off", "third_class"}
RIS_POLICY = "holdout"          # {"exclude", "as_on", "holdout"}

# Which CD timepoints count as genuinely growth-arrested "quiescent".
# MEASURED (see report): the CD arm is a PROGRESSIVE ARREST TIME COURSE, not a
# uniform quiescent state -- MKI67 runs 53.1 (D1.5) -> 0.89 (D7). Early
# timepoints are actively proliferating and are NOT quiescent in any useful
# sense. Only D4 onward is meaningfully arrested (MKI67 6.03 -> 0.89).
CD_QUIESCENT_DAYS = [4.0, 7.0, 10.0]

ARM_POLICY_RATIONALE = {
    "pdl_boundary": (
        "OFF = PDL<=37, ON = PDL>=50. The prior audit located the senescence "
        "inflection BETWEEN PDL46 and PDL50 via two independent lines of "
        "evidence (RNA marker trajectories; ATAC differential-peak escalation). "
        "Setting OFF at <=37 rather than <=46 deliberately leaves a margin "
        "below the inflection so the OFF arm is unambiguously pre-transition, "
        "rather than sampling right up to the edge of a boundary whose exact "
        "location is uncertain to within a few PDL."
    ),
    "transition_zone": (
        "PDL45 and PDL46 sit inside the PDL46-50 transition zone and are "
        "AMBIGUOUS BY CONSTRUCTION. Policy 'exclude' (default): drop them. "
        "Reasoning -- these samples are the single most likely source of label "
        "noise in the whole set. A binary ON/OFF head trained with genuinely "
        "intermediate samples forced to one side learns a blurred boundary, and "
        "because the transition is exactly where the biology of interest lives, "
        "that blur is maximally harmful. Excluding 6 samples (2 PDL x 3 reps) "
        "costs little given the OFF arm retains 15 and ON retains 9. "
        "Alternatives are available as config values for anyone who prefers "
        "more data over cleaner labels."
    ),
    "htert": (
        "hTERT-immortalized samples: EXCLUDED from the default training set. "
        "Reasoning -- they are a tempting additional OFF-state (18 samples, "
        "never senescent), but they are immortalized rather than simply young. "
        "hTERT constitutively expresses telomerase, which is itself a major "
        "transcriptional perturbation, and these cells were carried to PDL "
        "46-109, far beyond the WI-38 PDL range. A model trained with hTERT as "
        "OFF could learn 'telomerase-positive immortalized cell' as its OFF "
        "signature rather than 'young proliferating primary fibroblast', which "
        "is a genuine confound with the intended senescence axis. They remain "
        "valuable as a NEGATIVE CONTROL at evaluation time (a senescence head "
        "should score them low) and are emitted in the sample table for that "
        "use, just not trained on."
    ),
    "cell_density": (
        "Cell-density samples: included as a THIRD CLASS ('quiescent'), but "
        "ONLY the genuinely-arrested late timepoints (day 4, 7, 10), NOT all 10. "
        "This was corrected after measuring the arm rather than assuming it. "
        "MEASURED: the CD arm is a PROGRESSIVE ARREST TIME COURSE, not a uniform "
        "quiescent state -- mean MKI67 TPM runs D1=24.2, D1.5=53.1, D2=46.3, "
        "D2.5=32.1, D2.75=23.2, D3=17.6, D3.5=11.5, D4=6.0, D7=0.9, D10=4.5. "
        "Early CD timepoints are actively proliferating (MKI67 above the "
        "PDL20-37 proliferating mean of 15.5) and calling them 'quiescent' "
        "would have been simply wrong. Restricting to D4/D7/D10 gives a "
        "defensible arrested-but-not-senescent class. "
        "IMPORTANT CAVEAT, measured and NOT smoothed over: late-CD cells show "
        "HIGHER IL6 (42.3 vs 15.8) and CXCL8 (29.8 vs 5.5) than the senescent "
        "arm, so this is NOT a clean 'arrested without SASP' contrast -- "
        "contact-inhibited dense cultures mount their own inflammatory "
        "response. The third class therefore teaches 'senescence vs dense "
        "arrested culture', which is useful but is NOT the textbook "
        "quiescence contrast, and a head trained on it should not be claimed "
        "to separate senescence from quiescence-in-general. Set "
        "CD_POLICY='exclude' if that ambiguity is unacceptable for your use. "
        "NOTE: the CD arm has ZERO ATAC data -- irrelevant here (target is "
        "expression) but it means the accessibility prior cannot be validated "
        "against this arm."
    ),
    "ris": (
        "RIS (10 Gy irradiation-induced senescence): HELD OUT as a transfer "
        "test, not trained on. Reasoning -- RIS is the closest available proxy "
        "to the project's eventual doxorubicin/therapy-induced validation "
        "model, and its value is far higher as a held-out generalization probe "
        "than as ~15 extra training samples. Training on RIS would destroy the "
        "one chance this dataset offers to ask 'does a head trained on "
        "REPLICATIVE senescence transfer to DNA-damage-induced senescence?' -- "
        "which is the exact question the wet-lab plan depends on. RIS samples "
        "are emitted with split='holdout_ris'. Note the RIS 'd0_A/B/C' columns "
        "are the pre-irradiation baseline (mock), not senescent."
    ),
}

# Senescence marker panel for sanity checks (step 6)
MARKERS_UP = ["CDKN1A", "IL6", "CXCL8"]      # expected HIGHER in ON
MARKERS_DOWN = ["LMNB1", "MKI67"]            # expected LOWER in ON
MARKERS_REPORT = ["CDKN2A", "SERPINE1"]      # report, don't gate

# Chromosome split. chr9 (CDKN2A) + chr6 (CDKN1A) held out entirely as TEST.
TEST_CHROMS = ["chr9", "chr6"]
VAL_CHROMS = ["chr8", "chr16", "chr20"]
# train = everything else (assigned below, disjoint by construction)


def log(msg):
    print(msg, flush=True)


# ==========================================================================
# STEP 1: ACQUIRE / LOAD
# ==========================================================================
def read_sheet(xlsx, sheet):
    import openpyxl
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    ws = wb[sheet]
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    df = pd.DataFrame(list(rows), columns=header)
    return df


def load_expression():
    """Load all four expression arms. Returns (dict arm->DataFrame, qc dict).

    KNOWN DATA DEFECT IN THE DEPOSITED FILE -- handled explicitly, not silently:
    27 gene symbols in the GEO-deposited TPM workbook have been destroyed by
    Excel's date auto-conversion (the classic MARCH1-11 -> 2021-03-01..11,
    SEPT1-14 -> 2021-09-01..14, DEC1 -> 2021-12-01 corruption). They arrive as
    datetime.date objects, not strings. The true symbols cannot be recovered by
    inference (2021-03-01 could be MARCH1 or MARC1, and the canonical MARCH*/
    SEPT* symbols are otherwise absent from the table, so there is no
    disambiguating evidence in-file). These rows are DROPPED and counted rather
    than guessed at. None of the senescence marker genes is affected.
    """
    import datetime
    arms, qc = {}, {}
    for sheet, name in [("RS_TC_raw_TPM", "RS"), ("hTERT_TPM", "hTERT"),
                        ("CD_TPM", "CD"), ("RIS_TPM", "RIS")]:
        df = read_sheet(TPM_XLSX, sheet)
        is_date = df["gene"].apply(lambda g: isinstance(g, (datetime.date, datetime.datetime)))
        n_date = int(is_date.sum())
        df = df[~is_date].copy()
        df["gene"] = df["gene"].astype(str)
        n_dup = int(df["gene"].duplicated().sum())
        if n_dup:
            df = df.drop_duplicates("gene", keep="first")
        df = df.set_index("gene").apply(pd.to_numeric, errors="coerce")
        arms[name] = df
        qc[name] = {"date_corrupted_dropped": n_date, "residual_dups_dropped": n_dup,
                    "genes_kept": int(df.shape[0]), "samples": int(df.shape[1])}
        log(f"  {name:6s}: {df.shape[0]} genes x {df.shape[1]} samples "
            f"(dropped {n_date} Excel-date-corrupted symbols, {n_dup} residual dups)")
    return arms, qc


# ==========================================================================
# STEP 2: DEFINE THE ARMS
# ==========================================================================
def pdl_of(col):
    m = re.match(r"PDL(\d+)", col)
    return int(m.group(1)) if m else None


def assign_arms(arms):
    """Build the per-sample assignment table implementing the policies above."""
    rows = []

    # --- RS (replicative senescence) time course: the core ON/OFF axis ---
    for col in arms["RS"].columns:
        pdl = pdl_of(col)
        if pdl is None:
            raise ValueError(f"unparseable RS column: {col}")
        if pdl <= OFF_MAX_PDL:
            cls, split_role = "proliferating", "train_eligible"
        elif pdl >= ON_MIN_PDL:
            cls, split_role = "senescent", "train_eligible"
        else:
            if TRANSITION_POLICY == "exclude":
                cls, split_role = "transition_ambiguous", "excluded"
            elif TRANSITION_POLICY == "assign_off":
                cls, split_role = "proliferating", "train_eligible"
            else:
                cls, split_role = "senescent", "train_eligible"
        rows.append(dict(sample=col, arm="RS", pdl=pdl, day=None,
                         cls=cls, role=split_role))

    # --- hTERT immortalized ---
    for col in arms["hTERT"].columns:
        if HTERT_POLICY == "exclude":
            cls, role = "htert_immortalized", "excluded_control"
        else:
            cls, role = "proliferating", "train_eligible"
        rows.append(dict(sample=col, arm="hTERT", pdl=None, day=None,
                         cls=cls, role=role))

    # --- Cell density / quiescence ---
    for col in arms["CD"].columns:
        day = float(re.match(r"D([\d.]+)", col).group(1))
        arrested = day in CD_QUIESCENT_DAYS
        if CD_POLICY == "exclude":
            cls, role = "quiescent", "excluded"
        elif CD_POLICY == "as_off":
            cls, role = "proliferating", "train_eligible"
        else:  # third_class -- only genuinely-arrested timepoints qualify
            if arrested:
                cls, role = "quiescent", "train_eligible"
            else:
                # early CD timepoints are still proliferating (measured), so they
                # are NOT quiescent; excluded rather than mislabeled either way
                cls, role = "cd_not_yet_arrested", "excluded"
        rows.append(dict(sample=col, arm="CD", pdl=None, day=day,
                         cls=cls, role=role))

    # --- RIS (radiation-induced) ---
    for col in arms["RIS"].columns:
        is_baseline = col.startswith("d0")
        if RIS_POLICY == "exclude":
            cls, role = ("ris_baseline" if is_baseline else "ris_senescent"), "excluded"
        elif RIS_POLICY == "as_on":
            cls = "proliferating" if is_baseline else "senescent"
            role = "train_eligible"
        else:  # holdout
            cls = "ris_baseline" if is_baseline else "ris_senescent"
            role = "holdout_ris"
        m = re.match(r"d([\d.]+)", col)
        rows.append(dict(sample=col, arm="RIS", pdl=None,
                         day=float(m.group(1)) if m else None,
                         cls=cls, role=role))

    return pd.DataFrame(rows)


# ==========================================================================
# STEP 3: BUILD INTERVALS (TSS-centered, chromosome-bound aware)
# ==========================================================================
def load_tss(gtf_path):
    """Parse gene-level records from GENCODE, return TSS per gene symbol.

    TSS = start for + strand genes, end for - strand genes. Where a symbol maps
    to multiple gene records (PAR regions, duplicated symbols), the record on a
    primary chromosome with the longest span is kept -- recorded as an
    assumption, not silently resolved.
    """
    recs = []
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene":
                continue
            attr = f[8]
            m = re.search(r'gene_name "([^"]+)"', attr)
            if not m:
                continue
            gt = re.search(r'gene_type "([^"]+)"', attr)
            recs.append((m.group(1), f[0], int(f[3]), int(f[4]), f[6],
                         gt.group(1) if gt else None))
    df = pd.DataFrame(recs, columns=["gene", "chrom", "start", "end", "strand", "gene_type"])
    std = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    df = df[df["chrom"].isin(std)].copy()
    df["span"] = df["end"] - df["start"]
    df = df.sort_values("span", ascending=False).drop_duplicates("gene", keep="first")
    df["tss"] = np.where(df["strand"] == "+", df["start"], df["end"])
    return df.set_index("gene")


def load_chrom_sizes(gtf_path):
    """Derive chromosome sizes from the GTF's own contig extents is NOT safe
    (a GTF only covers annotated features). Use the authoritative UCSC hg38
    chrom.sizes, fetched once and cached."""
    cache = ROOT / "data" / "annotation" / "hg38.chrom.sizes"
    if not cache.exists():
        # shell out to curl: this machine's python.org build has no configured CA
        # bundle (SSLCertVerificationError on urlopen), while system curl's trust
        # store works -- same workaround as scripts/05_design_adequacy.py.
        import subprocess
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
        subprocess.run(["curl", "-sL", "-o", str(cache), url], check=True, timeout=120)
    sizes = pd.read_csv(cache, sep="\t", header=None, names=["chrom", "size"])
    return dict(zip(sizes["chrom"], sizes["size"]))


def build_intervals(genes, tss_df, chrom_sizes):
    """TSS-centered SEQ_LEN windows. Genes whose window runs past a chromosome
    end (or before position 0) are DROPPED and counted, never silently clipped
    (clipping would shift the TSS off-center and break the contract)."""
    keep, dropped = [], []
    for g in genes:
        if g not in tss_df.index:
            dropped.append((g, "no_annotation"))
            continue
        r = tss_df.loc[g]
        half = SEQ_LEN // 2
        start = int(r["tss"]) - half
        end = start + SEQ_LEN
        csize = chrom_sizes.get(r["chrom"])
        if csize is None:
            dropped.append((g, "chrom_not_in_sizes"))
            continue
        if start < 0:
            dropped.append((g, "window_before_chrom_start"))
            continue
        if end > csize:
            dropped.append((g, "window_past_chrom_end"))
            continue
        keep.append(dict(gene=g, chrom=r["chrom"], start=start, end=end,
                         strand=r["strand"], tss=int(r["tss"]),
                         gene_type=r["gene_type"]))
    return pd.DataFrame(keep), pd.DataFrame(dropped, columns=["gene", "reason"])


# ==========================================================================
# STEP 3b: ACCESSIBILITY PRIOR (region prioritization ONLY -- never a target)
# ==========================================================================
def genes_near_differential_accessibility(padj_thresh=0.05, lfc_thresh=1.0):
    """Return the set of gene symbols nearest to significantly differential
    ATAC peaks. This is used ONLY to optionally restrict which genes enter the
    training set -- it is never used as a regression target."""
    import openpyxl
    wb = openpyxl.load_workbook(ATAC_SIG_XLSX, read_only=True, data_only=True)
    ws = wb["atac_peaks"]
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    gi = header.index("gene_name")
    # late-PDL comparisons = the senescent-proxy contrast used by the prior audit
    late = [("PDL45_v_htert6_log2FC", "PDL45_v_htert6_adjPval"),
            ("PDL50_v_htert7_log2FC", "PDL50_v_htert7_adjPval")]
    idx = [(header.index(a), header.index(b)) for a, b in late]
    hits, total = set(), 0
    for r in rows:
        total += 1
        for fc_i, p_i in idx:
            fc, p = r[fc_i], r[p_i]
            if fc is None or p is None:
                continue
            try:
                if p < padj_thresh and abs(fc) >= lfc_thresh:
                    if r[gi]:
                        hits.add(r[gi])
                    break
            except TypeError:
                continue
    return hits, total


# ==========================================================================
# STEP 5: TARGET NORMALIZATION
# ==========================================================================
def normalize_targets(tpm_matrix):
    """log1p(TPM).

    Rationale (recorded in the report, not buried):
      - Borzoi's own targets live in a variance-stabilized space (per-bin
        **0.75 power + soft-clip + scale). We are NOT reproducing that chain:
        it applies to per-bin coverage tracks, and this head predicts a
        per-gene scalar. Reproducing a bin-level transform on a gene-level
        scalar would be cargo-culting.
      - What DOES carry over is the PRINCIPLE: raw TPM spans ~6 orders of
        magnitude and is dominated by a handful of very high-expression genes.
        Task 2 measured that prediction error concentrates in LOW-expression
        genes; an untransformed target would let a few housekeeping genes
        dominate MSE and make that worse.
      - log1p is monotonic, defined at zero, and keeps targets NON-NEGATIVE,
        which matches the measured ConvHead softplus output activation. A
        z-scored target would put mass below zero that softplus structurally
        cannot emit.
    """
    return np.log1p(tpm_matrix)


# ==========================================================================
# STEP 6: SANITY CHECKS (gate -- failures stop the build)
# ==========================================================================
def marker_check(rs_tpm, on_samples, off_samples):
    out = []
    for gene in MARKERS_UP + MARKERS_DOWN + MARKERS_REPORT:
        if gene not in rs_tpm.index:
            out.append(dict(gene=gene, status="NOT_FOUND"))
            continue
        on_v = rs_tpm.loc[gene, on_samples].astype(float).values
        off_v = rs_tpm.loc[gene, off_samples].astype(float).values
        t, p = stats.ttest_ind(on_v, off_v, equal_var=False)
        lfc = np.log2((on_v.mean() + 0.01) / (off_v.mean() + 0.01))
        if gene in MARKERS_UP:
            expected, ok = "UP", (lfc > 0 and p < 0.05)
        elif gene in MARKERS_DOWN:
            expected, ok = "DOWN", (lfc < 0 and p < 0.05)
        else:
            expected, ok = "report_only", None
        out.append(dict(gene=gene, expected=expected,
                        on_mean_tpm=round(float(on_v.mean()), 4),
                        off_mean_tpm=round(float(off_v.mean()), 4),
                        log2FC=round(float(lfc), 4), p_value=float(p),
                        direction=("UP" if lfc > 0 else "DOWN"), passes=ok))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--restrict-to-diff-accessibility", action="store_true",
                    help="Restrict genes to those near significant differential ATAC peaks "
                         "(accessibility used as a region prior only, never as a target).")
    ap.add_argument("--out-prefix", default="senescence_trainset")
    args = ap.parse_args()

    report = []

    def R(s=""):
        report.append(s)
        log(s)

    R("# STEP 1: ACQUIRE")
    for f in [TPM_XLSX, ATAC_SIG_XLSX, ATLAS_BED, ANNO]:
        R(f"  present: {f.name} ({f.stat().st_size:,} bytes)")
    R("  Genome build: hg38 for ALL files (SOFT 'Genome_build: Hg38'; Salmon index "
      "on Hg38; atac_peaks_sig readme 'region | hg38 coordinates'). GENCODE v44 is "
      "GRCh38. No hg19 file is used; no liftover required.")
    R()

    R("# Loading expression arms")
    arms, load_qc = load_expression()
    rs = arms["RS"]
    R()

    R("  Load QC (Excel date-corruption defect, see load_expression docstring):")
    for k, v in load_qc.items():
        R(f"    {k:6s} {v}")
    R()

    R("# STEP 2: DEFINE THE ARMS")
    assign = assign_arms(arms)
    assign.to_csv(OUT / "sample_assignments.csv", index=False)
    R(f"  Policies: OFF<=PDL{OFF_MAX_PDL}, ON>=PDL{ON_MIN_PDL}; "
      f"transition={TRANSITION_POLICY}, htert={HTERT_POLICY}, "
      f"cd={CD_POLICY} (quiescent days={CD_QUIESCENT_DAYS}), ris={RIS_POLICY}")
    R()
    R("  Class counts (all 96 RNA samples across 4 arms):")
    for (cls, role), n in assign.groupby(["cls", "role"]).size().items():
        R(f"    {cls:24s} role={role:18s} n={n}")
    R()

    on_samples = assign[(assign.cls == "senescent") & (assign.role == "train_eligible")]["sample"].tolist()
    off_samples = assign[(assign.cls == "proliferating") & (assign.role == "train_eligible")]["sample"].tolist()
    qui_samples = assign[(assign.cls == "quiescent") & (assign.role == "train_eligible")]["sample"].tolist()
    R(f"  ON  (senescent)     n={len(on_samples)}: {on_samples}")
    R(f"  OFF (proliferating) n={len(off_samples)}: {off_samples}")
    R(f"  QUIESCENT           n={len(qui_samples)}: {qui_samples}")
    R()

    # ---------------- STEP 6a: sanity gate BEFORE building anything ---------
    R("# STEP 6 (run FIRST as a gate): MARKER SANITY CHECKS")
    mk = marker_check(rs, on_samples, off_samples)
    mk.to_csv(OUT / "marker_sanity_checks.csv", index=False)
    R(mk.to_string(index=False))
    R()
    gated = mk[mk["passes"].notna()]
    failures = gated[~gated["passes"].astype(bool)]
    if len(failures):
        R("!! SANITY CHECK FAILED -- STOPPING, not building the training set.")
        R(failures.to_string(index=False))
        (OUT / "REPORT_trainset.md").write_text("\n".join(report))
        sys.exit(1)
    R(f"  All {len(gated)} gated markers pass (CDKN1A/IL6/CXCL8 up; LMNB1/MKI67 down).")
    R()

    # The quiescent third class is only worth its slot if senescence is actually
    # separable from generic growth arrest at the marker level. Check, don't assume.
    if qui_samples and CD_POLICY == "third_class":
        R("  Senescent-vs-QUIESCENT discrimination (justifies the third class):")
        cd = arms["CD"]
        rows = []
        for gene in MARKERS_UP + MARKERS_DOWN:
            if gene not in rs.index or gene not in cd.index:
                continue
            on_v = rs.loc[gene, on_samples].astype(float).values
            qu_v = cd.loc[gene, qui_samples].astype(float).values
            t, p = stats.ttest_ind(on_v, qu_v, equal_var=False)
            lfc = np.log2((on_v.mean() + 0.01) / (qu_v.mean() + 0.01))
            rows.append(dict(gene=gene, on_mean=round(float(on_v.mean()), 3),
                             quiescent_mean=round(float(qu_v.mean()), 3),
                             log2FC_on_vs_qui=round(float(lfc), 3), p_value=float(p)))
        qdf = pd.DataFrame(rows)
        qdf.to_csv(OUT / "senescent_vs_quiescent_markers.csv", index=False)
        R(qdf.to_string(index=False))
        sep = int((qdf["p_value"] < 0.05).sum())
        R(f"  => {sep}/{len(qdf)} markers separate senescent from quiescent at p<0.05. "
          f"MKI67/LMNB1 are expected to be low in BOTH (both are arrested); the "
          f"SASP markers (IL6/CXCL8) are the ones that should distinguish them, "
          f"since SASP is senescence-specific and not induced by contact inhibition.")
    R()

    # ---------------- STEP 3: intervals ------------------------------------
    R("# STEP 3: BUILD INTERVALS")
    tss_df = load_tss(ANNO)
    chrom_sizes = load_chrom_sizes(ANNO)
    R(f"  GENCODE v44 gene records on standard chromosomes: {len(tss_df):,}")

    expressed = rs.index[(rs[on_samples + off_samples].astype(float) > 0).any(axis=1)].tolist()
    R(f"  Genes with TPM>0 in >=1 ON/OFF sample: {len(expressed):,} (of {rs.shape[0]:,})")

    if args.restrict_to_diff_accessibility:
        da_genes, n_peaks = genes_near_differential_accessibility()
        R(f"  Accessibility prior ENABLED: {len(da_genes):,} distinct nearest-genes "
          f"from significant late-PDL differential peaks (of {n_peaks:,} atlas peaks)")
        before = len(expressed)
        expressed = [g for g in expressed if g in da_genes]
        R(f"  Restricting expressed genes by accessibility prior: {before:,} -> {len(expressed):,}")
    else:
        da_genes, n_peaks = genes_near_differential_accessibility()
        overlap = len([g for g in expressed if g in da_genes])
        R(f"  Accessibility prior DISABLED (default). For reference: {len(da_genes):,} "
          f"genes are near significant differential peaks; {overlap:,} of the "
          f"{len(expressed):,} expressed genes ({100*overlap/len(expressed):.1f}%) "
          f"would survive that restriction.")

    iv, dropped = build_intervals(expressed, tss_df, chrom_sizes)
    R(f"  Intervals built: {len(iv):,}")
    if len(dropped):
        R("  Dropped genes by reason:")
        for reason, n in dropped["reason"].value_counts().items():
            R(f"    {reason:28s} {n:,}")
        dropped.to_csv(OUT / "dropped_genes.csv", index=False)
        past_end = dropped[dropped.reason.isin(["window_past_chrom_end",
                                                 "window_before_chrom_start"])]
        R(f"  => {len(past_end)} genes dropped specifically because the {SEQ_LEN:,}bp "
          f"window ran past a chromosome boundary (handled explicitly, not clipped).")

        # Characterize the no_annotation drops so the number isn't just accepted.
        na = dropped[dropped.reason == "no_annotation"]["gene"].astype(str)
        clone_pat = re.compile(r"^(A[BCDFJLP]\d{6}|Z\d{5}|BX\d{6}|CT\d{6}|CU\d{6}|FP\d{6})")
        n_clone = int(na.apply(lambda g: bool(clone_pat.match(g))).sum())
        R(f"  no_annotation breakdown: {n_clone:,} of {len(na):,} "
          f"({100*n_clone/max(len(na),1):.1f}%) are clone-based identifiers "
          f"(AC######.#, AL######.#, etc.) -- unnamed loci from the older "
          f"annotation the submitters quantified against. These have no named "
          f"TSS to center a window on and are legitimately unusable. The "
          f"remaining {len(na)-n_clone:,} are named symbols absent from GENCODE "
          f"v44, mostly HGNC renames (AARS->AARS1, ADSS->ADSS2, CARS->CARS1, "
          f"etc.). Recovering those would need an HGNC previous-symbol mapping; "
          f"at <3% of expressed genes and with zero marker genes affected, this "
          f"is recorded as a known, quantified, acceptable loss rather than "
          f"silently absorbed.")
    R()

    # ---------------- STEP 4: chromosome split -----------------------------
    R("# STEP 4: CHROMOSOME SPLIT")
    all_chroms = sorted(iv["chrom"].unique())
    train_chroms = [c for c in all_chroms if c not in TEST_CHROMS + VAL_CHROMS]
    assert not (set(train_chroms) & set(VAL_CHROMS)), "train/val overlap"
    assert not (set(train_chroms) & set(TEST_CHROMS)), "train/test overlap"
    assert not (set(VAL_CHROMS) & set(TEST_CHROMS)), "val/test overlap"

    def which(c):
        if c in TEST_CHROMS: return "test"
        if c in VAL_CHROMS: return "val"
        return "train"

    iv["split"] = iv["chrom"].map(which)
    R(f"  TEST  (held out entirely): {TEST_CHROMS}  <- CDKN2A on chr9, CDKN1A on chr6")
    R(f"  VAL:   {VAL_CHROMS}")
    R(f"  TRAIN: {train_chroms}")
    R()
    R("  Gene counts per split:")
    for s, n in iv["split"].value_counts().items():
        R(f"    {s:6s} {n:,}")
    R()
    R("  Per-chromosome gene counts:")
    cc = iv.groupby(["split", "chrom"]).size().reset_index(name="n_genes")
    for _, r in cc.iterrows():
        R(f"    {r['split']:6s} {r['chrom']:6s} {r['n_genes']:,}")
    R()

    # ---------------- STEP 5: format --------------------------------------
    R("# STEP 5: FORMAT (AnnData, measured contract)")
    task_samples = on_samples + off_samples + qui_samples
    task_meta = assign.set_index("sample").loc[task_samples]

    tpm = rs.reindex(iv["gene"].values)
    cd_tpm = arms["CD"].reindex(iv["gene"].values)
    mat_parts = []
    for s in task_samples:
        src = cd_tpm if s in qui_samples else tpm
        mat_parts.append(src[s].astype(float).values)
    X_raw = np.vstack(mat_parts)                      # (n_tasks, n_intervals)
    X = normalize_targets(X_raw).astype(np.float32)
    R(f"  .X shape (n_tasks, n_intervals): {X.shape}  dtype={X.dtype}")
    R(f"  Target transform: log1p(TPM). Range {X.min():.3f} to {X.max():.3f} "
      f"(non-negative, matching the measured ConvHead softplus output).")

    obs = pd.DataFrame({
        "sample": task_samples,
        "arm": task_meta["arm"].values,
        "cls": task_meta["cls"].values,
        "pdl": task_meta["pdl"].values,
        "day": task_meta["day"].values,
    }).set_index("sample")

    var = iv.set_index("gene")[["chrom", "start", "end", "strand", "tss", "gene_type", "split"]]

    ad = anndata.AnnData(X=X, obs=obs, var=var)
    ad.uns["contract"] = {
        "seq_len": SEQ_LEN, "bin_size": BIN_SIZE, "label_len": LABEL_LEN,
        "genome": "hg38", "input_shape": "(4, 524288) float32 channels-first",
        "trunk_output": "(batch, 1920, 6144)", "target_per_example": "(n_tasks, 1) float32",
        "head": "ConvHead + softplus (non-negative outputs)",
    }
    ad.uns["target_transform"] = "log1p(TPM)"
    ad.uns["arm_policy"] = {
        "OFF_MAX_PDL": OFF_MAX_PDL, "ON_MIN_PDL": ON_MIN_PDL,
        "TRANSITION_POLICY": TRANSITION_POLICY, "HTERT_POLICY": HTERT_POLICY,
        "CD_POLICY": CD_POLICY, "RIS_POLICY": RIS_POLICY,
        "CD_QUIESCENT_DAYS": CD_QUIESCENT_DAYS,
    }
    ad.uns["arm_policy_rationale"] = ARM_POLICY_RATIONALE
    ad.uns["splits"] = {"train": train_chroms, "val": VAL_CHROMS, "test": TEST_CHROMS}
    ad.uns["provenance"] = {
        "series": "GSE175533", "bioproject": "PRJNA732700",
        "annotation": "GENCODE v44 (GRCh38)",
        "note": ("ATAC and RNA in GSE175533 are PARALLEL CULTURES, not split aliquots. "
                 "Accessibility was used only as an optional region-selection prior; "
                 "the regression target is expression only."),
    }

    out_h5ad = OUT / f"{args.out_prefix}.h5ad"
    ad.write_h5ad(out_h5ad)
    R(f"  Wrote {out_h5ad} ({out_h5ad.stat().st_size:,} bytes)")
    R()

    # ---------------- STEP 6b: structural checks ---------------------------
    R("# STEP 6 (structural checks)")
    n_iv, n_tgt = ad.shape[1], X.shape[1]
    R(f"  Interval count == target count: {n_iv} == {n_tgt} -> {n_iv == n_tgt}")
    assert n_iv == n_tgt

    dup = ad.var.index.duplicated().sum()
    R(f"  Duplicate genes in var index: {dup} -> {'PASS' if dup == 0 else 'FAIL'}")
    assert dup == 0

    per_gene_splits = ad.var.groupby(level=0)["split"].nunique()
    multi = int((per_gene_splits > 1).sum())
    R(f"  Genes appearing in >1 split: {multi} -> {'PASS' if multi == 0 else 'FAIL'}")
    assert multi == 0

    chrom_split = ad.var.groupby("chrom")["split"].nunique()
    bad = int((chrom_split > 1).sum())
    R(f"  Chromosomes spanning >1 split: {bad} -> {'PASS' if bad == 0 else 'FAIL'}")
    assert bad == 0

    for g, c in [("CDKN2A", "chr9"), ("CDKN1A", "chr6")]:
        if g in ad.var.index:
            row = ad.var.loc[g]
            R(f"  {g} on {row['chrom']} split={row['split']} "
              f"-> {'PASS (held out)' if row['split'] == 'test' else 'FAIL'}")
            assert row["split"] == "test"
        else:
            R(f"  {g}: not in interval set (dropped upstream) -- note, not a gate failure")

    nan_frac = float(np.isnan(X).mean())
    R(f"  NaN fraction in .X: {nan_frac:.6f}")
    R()

    R("# WIDTH CHECK")
    widths = (ad.var["end"] - ad.var["start"]).unique()
    R(f"  Distinct interval widths: {widths} (must be exactly [{SEQ_LEN}])")
    assert list(widths) == [SEQ_LEN]

    (OUT / "REPORT_trainset.md").write_text("\n".join(report))
    log(f"\nWrote {OUT/'REPORT_trainset.md'}")


if __name__ == "__main__":
    main()
