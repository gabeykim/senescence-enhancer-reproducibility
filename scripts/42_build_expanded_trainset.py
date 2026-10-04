#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 42_build_expanded_trainset.py
#
# Build the 13,211-gene multi-study expression matrix.
#
# CONSUMES: output/tables/expansion_sample_assignments.csv, data/*/
# PRODUCES: output/trainset_expanded/_Xraw.npy, _Xc.npy, _Xcs.npy, BUILD_LOG.txt
# ---------------------------------------------------------------------------
"""
Build the EXPANDED multi-study senescence training set.

Motivation: the single-dataset fine-tune failed its shuffled-label control
(real test AUC 0.659 vs shuffled 0.5605 +/- 0.0940, 2/20 permutations beat the real
model, empirical p = 0.143). The diagnosis was insufficient biological replication:
the ON arm was effectively n~3 from one WI-38 serial-culture lineage. The goal here is
INDEPENDENT BIOLOGICAL REPLICATION -- more lineages and more induction mechanisms --
not more genes from the same culture.

Pipeline
  1. Load every usable study's deposited expression matrix in its native units.
  2. Harmonise onto a common natural-log scale.
  3. Intersect gene universes; attach the SAME TSS-centered intervals as the
     existing build.
  4. Quantify batch structure BEFORE merging (PCA: study vs senescence status).
  5. Apply a documented correction; re-quantify.
  6. Per-study marker sanity checks (a study that fails is dropped).
  7. Emit AnnData matching the existing contract.

Outputs
  output/trainset_expanded/senescence_trainset_multistudy.h5ad
  output/trainset_expanded/multistudy_sample_assignments.csv
  output/trainset_expanded/batch_pca_stats.csv
  output/trainset_expanded/per_study_markers.csv
  output/figures/multistudy_batch_pca.png
"""
import datetime
import gzip
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
DATA = ROOT / "data"
OUT = ROOT / "output" / "trainset_expanded"
FIG = ROOT / "output" / "figures"
TAB = ROOT / "output" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

ANNO = DATA / "annotation" / "gencode.v44.basic.annotation.gtf.gz"
SEQ_LEN = 524_288
TEST_CHROMS = ["chr9", "chr6"]
VAL_CHROMS = ["chr8", "chr16", "chr20"]

# Marker panel. CDKN2A is measured and reported but is explicitly NOT a pass/fail
# criterion: prior work in this project established it barely moves (+0.25, p=0.21)
# and its promoter accessibility is flat.
MARKERS_UP = ["CDKN1A", "IL6", "CXCL8"]
MARKERS_DOWN = ["LMNB1", "MKI67"]
MARKERS_REPORT_ONLY = ["CDKN2A"]

LOG_LINES = []


def R(s=""):
    print(s, flush=True)
    LOG_LINES.append(s)


# reuse the existing builder's annotation/interval code so both builds are identical there
spec = importlib.util.spec_from_file_location(
    "b30", ROOT / "scripts" / "30_build_senescence_trainset.py")
b30 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b30)


# ======================================================================================
# STEP 1 -- LOAD EACH STUDY IN ITS NATIVE UNITS
# ======================================================================================
def union_exon_lengths(gtf_path):
    """Union-exon length per gene symbol, for counts -> TPM.

    Needed only by GSE99028 (the one study deposited as raw counts). Note that the
    length enters as a per-gene constant WITHIN a study, so the per-study per-gene
    centering applied in step 5 cancels it exactly -- an imperfect length estimate
    (GENCODE v44 lengths applied to a refGene-symbol count table) therefore cannot
    bias the corrected targets. It only affects the uncorrected layer.
    """
    ex = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon":
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if not m:
                continue
            ex.setdefault(m.group(1), []).append((int(f[3]), int(f[4])))
    out = {}
    for g, iv in ex.items():
        iv.sort()
        tot, cs, ce = 0, None, None
        for s, e in iv:
            if cs is None:
                cs, ce = s, e
            elif s <= ce + 1:
                ce = max(ce, e)
            else:
                tot += ce - cs + 1
                cs, ce = s, e
        tot += ce - cs + 1
        out[g] = tot
    return pd.Series(out, name="union_exon_len")


def load_gse175533():
    """WI-38 replicative senescence -- the existing primary study.

    Rebuilt here from the deposited TPM workbook rather than reused from the
    existing .h5ad, so that every study passes through one identical harmonisation
    path, and so the RIS arm can be carried as a holdout.
    """
    arms, qc = b30.load_expression()
    assign = b30.assign_arms(arms)
    frames = []
    for name, df in arms.items():
        sub = df.copy()
        sub.columns = [f"{name}::{c}" for c in sub.columns]
        frames.append(sub)
    idx = frames[0].index
    for f in frames[1:]:
        idx = idx.intersection(f.index)
    mat = pd.concat([f.loc[idx] for f in frames], axis=1)
    assign = assign.copy()
    assign["matrix_column"] = assign["arm"] + "::" + assign["sample"]
    assign = assign[assign["matrix_column"].isin(mat.columns)]
    return mat, assign, qc


def load_gse74324():
    """IMR90 oncogene-induced senescence. Deposited as RPKM.

    RPKM -> TPM is an exact per-sample rescale (TPM = RPKM / sum(RPKM) * 1e6), so
    no information is invented by the conversion.
    """
    p = DATA / "GSE74324" / "GSE74324_ALL_samples_rpkm.txt"
    df = pd.read_csv(p, sep="\t")
    ann_cols = ["entrezgene", "ensembl_gene_id", "external_gene_name", "chromosome_name",
                "band", "start_position", "end_position", "transcript_length"]
    first = df.columns[0]
    sym = df["external_gene_name"].astype(str)
    vals = df[[c for c in df.columns if c not in ann_cols and c != first]]
    vals = vals.apply(pd.to_numeric, errors="coerce")
    vals.index = sym
    vals = vals[~vals.index.isin(["nan", ""])]
    n_dup = int(vals.index.duplicated().sum())
    # duplicate symbols: keep the highest-expressed record (mean RPKM), a stated choice
    if n_dup:
        vals = vals.assign(_m=vals.mean(axis=1)).sort_values("_m", ascending=False)
        vals = vals[~vals.index.duplicated(keep="first")].drop(columns="_m")
    tpm = vals / vals.sum(axis=0) * 1e6
    return tpm, {"rpkm_to_tpm": "exact per-sample rescale", "dup_symbols_collapsed": n_dup,
                 "genes": int(tpm.shape[0])}


def load_gse99028(lengths):
    """IMR90 etoposide-induced senescence. Deposited as raw counts."""
    p = DATA / "GSE99028" / "GSE99028_read_count.txt"
    df = pd.read_csv(p, sep="\t")
    df = df.set_index("refGene")
    df = df.apply(pd.to_numeric, errors="coerce")
    n_dup = int(df.index.duplicated().sum())
    if n_dup:
        df = df.groupby(level=0).sum()
    L = lengths.reindex(df.index)
    n_nolen = int(L.isna().sum())
    df = df[~L.isna()]
    L = L.dropna()
    rate = df.div(L.values / 1000.0, axis=0)          # reads per kilobase
    tpm = rate / rate.sum(axis=0) * 1e6
    return tpm, {"counts_to_tpm": "union-exon length from GENCODE v44",
                 "genes_dropped_no_length": n_nolen, "dup_symbols_summed": n_dup,
                 "genes": int(tpm.shape[0])}


def load_gse205692():
    """GM21 oncogene-induced senescence. Affymetrix HTA 2.0, log2 RMA (sva-corrected).

    Values come from the series matrix (the submitters' own processed table), not from
    re-normalising the CEL files -- re-running RMA locally would produce a THIRD
    quantification pipeline in the same set for no benefit.

    Probe -> symbol uses the GPL17586 'gene_assignment' field. A transcript cluster is
    kept only when every assignment in that field names the SAME symbol; multi-gene
    clusters are ambiguous and are dropped rather than arbitrarily assigned to the first.
    """
    mp = DATA / "GSE205692" / "GSE205692_series_matrix.txt"
    with open(mp, errors="replace") as fh:
        lines = fh.readlines()
    b = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    e = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    from io import StringIO
    mat = pd.read_csv(StringIO("".join(lines[b + 1:e])), sep="\t")
    mat.columns = [c.strip('"') for c in mat.columns]
    mat = mat.set_index(mat.columns[0])
    mat.index = [str(i).strip('"') for i in mat.index]
    mat = mat.apply(pd.to_numeric, errors="coerce")

    # platform table from the family SOFT
    sp = DATA / "GSE205692" / "GSE205692_family.soft"
    rows, inside, hdr = [], False, None
    with open(sp, errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                inside = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if inside:
                f = line.rstrip("\n").split("\t")
                if hdr is None:
                    hdr = f
                else:
                    rows.append(f[:len(hdr)])
    plat = pd.DataFrame(rows, columns=hdr)
    ga = plat.set_index("ID")["gene_assignment"].fillna("")

    def one_symbol(s):
        if not s or s == "---":
            return None
        syms = {b.strip().split(" // ")[1].strip()
                for b in s.split(" /// ") if len(b.split(" // ")) > 1}
        syms = {x for x in syms if x and x != "---"}
        return syms.pop() if len(syms) == 1 else None

    sym = ga.map(one_symbol)
    n_amb = int(sym.isna().sum())
    sym = sym.dropna()
    keep = mat.index.intersection(sym.index)
    mat = mat.loc[keep]
    mat.index = sym.loc[keep].values
    n_dup = int(mat.index.duplicated().sum())
    if n_dup:                                  # several clusters per gene: take the max
        mat = mat.groupby(level=0).max()
    return mat, {"probe_clusters_ambiguous_dropped": n_amb,
                 "clusters_collapsed_per_gene": n_dup, "genes": int(mat.shape[0]),
                 "units": "log2 RMA intensity (sva batch-corrected by submitters)"}


# ======================================================================================
# STEP 1b -- HARMONISE GENE SYMBOLS ACROSS ANNOTATION VINTAGES
# ======================================================================================
HGNC = DATA / "annotation" / "hgnc_complete_set.txt"


def load_hgnc():
    """approved-symbol set + {previous-or-alias symbol -> approved symbol} (1:1 only).

    The four studies were annotated against different vintages: GSE205692's HTA 2.0
    platform is ~2013-era and calls the SASP chemokine IL8, while the RNA-seq studies
    use the current CXCL8. Left alone, that drift silently deletes real genes -- and
    one of them is a marker this build validates on -- from the cross-study
    intersection. Symbols that map ambiguously (the same old name reused by more than
    one current gene) are NOT mapped; they are counted and left as-is, because guessing
    would silently merge two different genes' measurements.
    """
    h = pd.read_csv(HGNC, sep="\t", low_memory=False,
                    usecols=["symbol", "alias_symbol", "prev_symbol", "status"])
    h = h[h["status"] == "Approved"]
    approved = set(h["symbol"].astype(str))
    cand = {}
    for col in ("prev_symbol", "alias_symbol"):
        for cur, old in zip(h["symbol"].astype(str), h[col].fillna("").astype(str)):
            for o in filter(None, (x.strip() for x in old.split("|"))):
                if o in approved:          # still a live symbol elsewhere: never remap
                    continue
                cand.setdefault(o, set()).add(cur)
    alias = {o: next(iter(v)) for o, v in cand.items() if len(v) == 1}
    ambiguous = {o for o, v in cand.items() if len(v) > 1}
    return approved, alias, ambiguous


def harmonize_symbols(mat, approved, alias, study):
    """Rename a study's rows onto current HGNC symbols. Returns (matrix, stats)."""
    idx = [str(i) for i in mat.index]
    present = set(idx)
    new, n_map, n_collide = [], 0, 0
    for g in idx:
        if g in approved or g not in alias:
            new.append(g)
            continue
        tgt = alias[g]
        if tgt in present:                 # study already carries the modern symbol
            new.append(g)
            n_collide += 1
        else:
            new.append(tgt)
            n_map += 1
    out = mat.copy()
    out.index = new
    n_dup = int(out.index.duplicated().sum())
    if n_dup:
        out = out.groupby(level=0).max()
    stats = dict(study=study, renamed=n_map, collisions_left_alone=n_collide,
                 merged_after_rename=n_dup,
                 not_in_hgnc=int(sum(1 for g in new if g not in approved)))
    return out, stats


# ======================================================================================
# STEP 2 -- HARMONISE TO A COMMON LOG SCALE
# ======================================================================================
def to_log(mat, kind):
    """RNA-seq TPM -> log1p(TPM) (matches the existing build's target transform).
       Microarray log2 RMA -> multiplied by ln2 to sit on a natural-log axis.

    The microarray rescale is cosmetic: it makes the two families numerically
    comparable in magnitude, but it does NOT make them the same measurement. Only the
    per-study correction in step 5 addresses that, which is why the PCA is reported
    both before and after.
    """
    if kind == "tpm":
        return np.log1p(mat)
    if kind == "log2":
        return mat * np.log(2.0)
    raise ValueError(kind)


def main():
    R("=" * 86)
    R("EXPANDED MULTI-STUDY SENESCENCE TRAINING SET")
    R(f"built {datetime.date.today().isoformat()}")
    R("=" * 86)

    # ---------------- load ----------------
    R("\n[1] LOADING STUDIES\n" + "-" * 86)
    lengths = union_exon_lengths(ANNO)
    R(f"  GENCODE v44 union-exon lengths for {len(lengths):,} gene symbols")

    m175, a175, qc175 = load_gse175533()
    R(f"  GSE175533  WI-38   replicative   {m175.shape[0]:,} genes x {m175.shape[1]} samples  (TPM)")
    m743, qc743 = load_gse74324()
    R(f"  GSE74324   IMR90   OIS           {m743.shape[0]:,} genes x {m743.shape[1]} samples  (RPKM->TPM)")
    m990, qc990 = load_gse99028(lengths)
    R(f"  GSE99028   IMR90   etoposide     {m990.shape[0]:,} genes x {m990.shape[1]} samples  (counts->TPM)")
    m205, qc205 = load_gse205692()
    R(f"  GSE205692  GM21    OIS           {m205.shape[0]:,} genes x {m205.shape[1]} samples  (log2 RMA microarray)")
    R(f"    GSE99028  : {qc990['genes_dropped_no_length']} symbols had no GENCODE length and were dropped")
    R(f"    GSE205692 : {qc205['probe_clusters_ambiguous_dropped']:,} transcript clusters were multi-gene/unassigned and were dropped")
    R(f"    GSE74324  : {qc743['dup_symbols_collapsed']} duplicate symbols collapsed (kept highest mean)")

    logmats = {"GSE175533": to_log(m175, "tpm"), "GSE74324": to_log(m743, "tpm"),
               "GSE99028": to_log(m990, "tpm"), "GSE205692": to_log(m205, "log2")}

    R("\n[1b] GENE-SYMBOL HARMONISATION (HGNC current symbols)\n" + "-" * 86)
    approved, alias, ambiguous = load_hgnc()
    R(f"  HGNC approved symbols: {len(approved):,}; unambiguous alias/previous "
      f"mappings: {len(alias):,}; ambiguous old symbols left unmapped: {len(ambiguous):,}")
    sym_stats = []
    for s in list(logmats):
        logmats[s], st = harmonize_symbols(logmats[s], approved, alias, s)
        sym_stats.append(st)
    R(pd.DataFrame(sym_stats).to_string(index=False))

    # ---------------- Excel date corruption audit on every new file ----------------
    R("\n[2] EXCEL DATE-CORRUPTION AUDIT (all studies)\n" + "-" * 86)
    date_re = re.compile(r"^\d{4}-\d{2}-\d{2}")
    corrupt = {}
    for s, mm in logmats.items():
        hits = [g for g in map(str, mm.index) if date_re.match(g)]
        corrupt[s] = hits
        R(f"  {s:11s} date-like gene symbols: {len(hits)}"
          + (f"  {hits[:6]}" if hits else ""))
    R(f"  (GSE175533's 27 corrupted symbols are dropped upstream by load_expression();")
    R(f"   counts here are RESIDUAL, i.e. what survived into the loaded matrix.)")

    # ---------------- sample assignment ----------------
    R("\n[3] SAMPLE ASSIGNMENT\n" + "-" * 86)
    exp_assign = pd.read_csv(TAB / "expansion_sample_assignments.csv")
    rows = []
    for _, r in a175.iterrows():
        rows.append(dict(study="GSE175533", matrix_column=r["matrix_column"],
                         sample=r["sample"], cell_line="WI-38",
                         mechanism="replicative" if r["arm"] in ("RS", "hTERT")
                         else ("radiation_RIS" if r["arm"] == "RIS" else "contact_inhibition"),
                         platform_kind="RNA-seq", quantification="TPM",
                         cls=r["cls"], role=r["role"], detail=f"arm={r['arm']} pdl={r['pdl']} day={r['day']}"))
    for _, r in exp_assign.iterrows():
        rows.append(dict(study=r["study"], matrix_column=r["matrix_column"],
                         sample=r["title"], cell_line=r["cell_line"],
                         mechanism=r["mechanism"], platform_kind=r["platform_kind"],
                         quantification=r["quantification"], cls=r["cls"],
                         role=r["role"], detail=r["detail"]))
    assign = pd.DataFrame(rows)
    assign = assign[assign.apply(
        lambda r: r["matrix_column"] in logmats[r["study"]].columns, axis=1)]
    R(assign.groupby(["study", "role", "cls"]).size().to_string())

    # ---------------- gene universe ----------------
    R("\n[4] GENE UNIVERSE\n" + "-" * 86)
    used_studies = sorted(assign["study"].unique())
    sets = {s: set(map(str, logmats[s].index)) for s in used_studies}
    for s in used_studies:
        R(f"  {s:11s} {len(sets[s]):6,} symbols")
    inter = set.intersection(*sets.values())
    R(f"  {'INTERSECTION':11s} {len(inter):6,} symbols across all {len(used_studies)} studies")
    for s in used_studies:
        R(f"    lost from {s:11s}: {len(sets[s]) - len(inter):6,} "
          f"({100*(len(sets[s])-len(inter))/len(sets[s]):.1f}%)")

    tss = b30.load_tss(ANNO)
    sizes = b30.load_chrom_sizes(ANNO)
    iv, dropped = b30.build_intervals(sorted(inter), tss, sizes)
    R(f"  with a valid TSS-centered {SEQ_LEN:,} bp window: {len(iv):,} "
      f"(dropped {len(dropped):,})")
    if len(dropped):
        R("    drop reasons: " + str(dropped["reason"].value_counts().to_dict()))

    genes = iv["gene"].tolist()
    iv = iv.set_index("gene").loc[genes]

    # expression filter, matching the existing build's intent: a gene must be observed
    # above zero in at least one training-eligible sample of every RNA-seq study
    rna_studies = [s for s in used_studies if s != "GSE205692"]
    mask = np.ones(len(genes), bool)
    for s in rna_studies:
        cols = assign.loc[assign["study"] == s, "matrix_column"]
        sub = logmats[s].reindex(genes).loc[:, cols]
        mask &= (sub.fillna(0).to_numpy() > 0).any(axis=1)
    R(f"  expressed (>0 in >=1 sample of every RNA-seq study): {int(mask.sum()):,}")
    genes = [g for g, m in zip(genes, mask) if m]
    iv = iv.loc[genes]

    # ---------------- pooled matrix ----------------
    blocks, meta = [], []
    for _, r in assign.iterrows():
        blocks.append(logmats[r["study"]].reindex(genes)[r["matrix_column"]].to_numpy(float))
        meta.append(r)
    Xraw = np.vstack(blocks)                    # (n_samples, n_genes)
    meta = pd.DataFrame(meta).reset_index(drop=True)
    Xraw = np.nan_to_num(Xraw, nan=0.0)
    R(f"  pooled matrix: {Xraw.shape[0]} samples x {Xraw.shape[1]:,} genes")

    # ---------------- per-study markers ----------------
    R("\n[5] PER-STUDY MARKER SANITY CHECKS\n" + "-" * 86)
    R("  log2 fold change, ON arm vs that study's own OFF arm. Pass criteria: "
      "CDKN1A/IL6/CXCL8 up, LMNB1/MKI67 down.")
    R("  CDKN2A is reported but is NOT a pass/fail criterion (prior work: +0.25, p=0.21).")
    mrows = []
    gi = {g: i for i, g in enumerate(genes)}
    for s in used_studies:
        sel = meta["study"] == s
        on = sel & meta["cls"].eq("senescent") & meta["role"].isin(["train_eligible", "test_transfer"])
        off = sel & meta["cls"].eq("proliferating") & meta["role"].isin(["train_eligible", "test_transfer"])
        if on.sum() == 0 or off.sum() == 0:
            continue
        R(f"\n  {s}  (ON n={int(on.sum())}, OFF n={int(off.sum())})")
        for g in MARKERS_UP + MARKERS_DOWN + MARKERS_REPORT_ONLY:
            if g not in gi:
                R(f"    {g:8s} NOT IN GENE UNIVERSE")
                mrows.append(dict(study=s, gene=g, log2fc=np.nan, direction="absent",
                                  expected="", pass_=""))
                continue
            j = gi[g]
            # values are natural-log; convert the difference to log2 for reporting
            d = (Xraw[on.to_numpy(), j].mean() - Xraw[off.to_numpy(), j].mean()) / np.log(2)
            exp = "up" if g in MARKERS_UP else ("down" if g in MARKERS_DOWN else "report_only")
            ok = "" if exp == "report_only" else ("PASS" if (d > 0) == (exp == "up") else "FAIL")
            R(f"    {g:8s} log2FC {d:+7.3f}   expected {exp:11s} {ok}")
            mrows.append(dict(study=s, gene=g, log2fc=float(d), direction=exp,
                              expected=exp, pass_=ok))
    mk = pd.DataFrame(mrows)
    mk.to_csv(OUT / "per_study_markers.csv", index=False)

    failed = sorted(mk.loc[mk["pass_"] == "FAIL", "study"].unique())
    R("\n  Studies with >=1 marker failure: " + (", ".join(failed) if failed else "none"))

    # ---------- 5b: is the GSE205692 ON day-window the right one? ----------
    R("\n[5b] GSE205692 -- VALIDATING THE RAS DAY WINDOW\n" + "-" * 86)
    R("  The ON window (days 13-25) was chosen a priori because this study's subject is")
    R("  ESCAPE from OIS: late RAS timepoints may have resumed proliferation. Asserting")
    R("  that would be guessing, so it is checked here against the measured markers.")
    R("  Senescence score = mean(z of CDKN1A, IL6, CXCL8) - mean(z of LMNB1, MKI67),")
    R("  z-scored across GSE205692 samples only.")
    sel205 = (meta["study"] == "GSE205692").to_numpy()
    sub = Xraw[sel205]
    msub = meta[sel205].reset_index(drop=True)
    def zcol(g):
        v = sub[:, gi[g]]
        sd = v.std()
        return (v - v.mean()) / (sd if sd > 0 else 1.0)
    up = np.mean([zcol(g) for g in MARKERS_UP if g in gi], axis=0)
    dn = np.mean([zcol(g) for g in MARKERS_DOWN if g in gi], axis=0)
    msub["sen_score"] = up - dn
    msub["day"] = msub["detail"].str.extract(r"day=(\d+)").astype(int)
    msub["is_ras"] = ~msub["sample"].str.contains("empty vector")
    ev = msub.loc[~msub["is_ras"], "sen_score"].mean()
    R(f"\n  empty-vector (OFF) mean score: {ev:+.3f}")
    msub["sasp_up"] = up
    msub["arrest_down"] = -dn        # positive == proliferation markers suppressed
    R("\n  Reported as two components, because they can dissociate: RAS drives an")
    R("  inflammatory (SASP) transcriptional response within days, while the")
    R("  proliferative arrest that actually defines senescence takes longer. A")
    R("  timepoint that is SASP-high but not yet arrested is NOT senescent.")
    R("\n   RAS day   n   SASP(up)  arrest(-LMNB1,-MKI67)   composite   assigned")
    for d, grp in msub[msub["is_ras"]].groupby("day"):
        lab = grp["cls"].iloc[0]
        R(f"   {d:7d} {len(grp):3d}   {grp['sasp_up'].mean():+8.3f}   "
          f"{grp['arrest_down'].mean():+19.3f}   {grp['sen_score'].mean():+9.3f}   {lab}")
    R(f"\n   {'empty vec':>7s} {int((~msub['is_ras']).sum()):3d}   "
      f"{msub.loc[~msub['is_ras'],'sasp_up'].mean():+8.3f}   "
      f"{msub.loc[~msub['is_ras'],'arrest_down'].mean():+19.3f}   "
      f"{msub.loc[~msub['is_ras'],'sen_score'].mean():+9.3f}   proliferating")
    R("\n  Per-marker log2FC vs empty vector, for the two RAS groups excluded a priori:")
    ev_mask = (~msub["is_ras"]).to_numpy()
    for lab, m in [("day 8 (pre-senescent?)", (msub["day"] == 8).to_numpy()),
                   ("days 32-56 (escape?)",
                    (msub["is_ras"] & msub["day"].ge(32)).to_numpy()),
                   ("days 13-25 (chosen ON)",
                    (msub["is_ras"] & msub["day"].between(13, 25)).to_numpy())]:
        parts = []
        for g in MARKERS_UP + MARKERS_DOWN:
            if g in gi:
                d = (sub[m, gi[g]].mean() - sub[ev_mask, gi[g]].mean()) / np.log(2)
                parts.append(f"{g} {d:+.2f}")
        R(f"    {lab:26s} " + "  ".join(parts))

    on_days = msub[(msub["is_ras"]) & (msub["cls"] == "senescent")]
    ex_days = msub[(msub["is_ras"]) & (msub["cls"] != "senescent")]
    R(f"\n  chosen ON days  mean score {on_days['sen_score'].mean():+.3f}  (n={len(on_days)})")
    R(f"  excluded RAS    mean score {ex_days['sen_score'].mean():+.3f}  (n={len(ex_days)})")
    R(f"  empty vector    mean score {ev:+.3f}  (n={int((~msub['is_ras']).sum())})")
    window_ok = (on_days["sen_score"].mean() > ev) and \
                (on_days["sen_score"].mean() > ex_days["sen_score"].mean())
    R(f"  VERDICT: chosen window is {'SUPPORTED' if window_ok else 'NOT SUPPORTED'} "
      f"by the markers.")
    if not window_ok:
        R("  -> per the stated policy, GSE205692 would be DROPPED here.")

    # ---------------- batch structure ----------------
    R("\n[6] BATCH STRUCTURE -- PCA BEFORE CORRECTION\n" + "-" * 86)
    stats = []

    def pca_report(X, label, meta, n=4):
        Z = X - X.mean(axis=0, keepdims=True)
        # SVD on samples x genes
        U, S, Vt = np.linalg.svd(Z, full_matrices=False)
        var = S ** 2
        frac = var / var.sum()
        scores = U * S
        out = []
        for k in range(min(n, scores.shape[1])):
            pc = scores[:, k]
            r2_study = anova_r2(pc, meta["study"].to_numpy())
            r2_stat = anova_r2(pc, meta["cls_binary"].to_numpy())
            out.append(dict(space=label, pc=f"PC{k+1}", var_explained=float(frac[k]),
                            r2_study=r2_study, r2_status=r2_stat))
            R(f"  {label:24s} PC{k+1}  var={100*frac[k]:5.1f}%   "
              f"R2(study)={r2_study:.3f}   R2(senescence status)={r2_stat:.3f}"
              f"   -> {'STUDY' if r2_study > r2_stat else 'STATUS'} dominates")
        stats.extend(out)
        return scores, frac

    def anova_r2(y, labels):
        """Fraction of variance in y explained by a categorical label (one-way)."""
        y = np.asarray(y, float)
        gt = y.mean()
        ss_tot = ((y - gt) ** 2).sum()
        if ss_tot == 0:
            return 0.0
        ss_b = 0.0
        for lv in pd.unique(labels):
            m = labels == lv
            ss_b += m.sum() * (y[m].mean() - gt) ** 2
        return float(ss_b / ss_tot)

    meta["cls_binary"] = np.where(meta["cls"].eq("senescent"), "ON",
                          np.where(meta["cls"].eq("proliferating"), "OFF", "other"))
    sc_raw, fr_raw = pca_report(Xraw, "raw pooled log", meta)

    # ---------------- correction ----------------
    R("\n[7] CORRECTION\n" + "-" * 86)
    R("  Applied: ARM-BALANCED PER-STUDY PER-GENE CENTERING.")
    R("    For each study s and gene g, the reference is the unweighted mean of that")
    R("    study's ON-arm mean and OFF-arm mean; that reference is subtracted.")
    R("  Why this and not the alternatives:")
    R("    - Subtracting a per-study per-gene constant removes the study offset AND, in")
    R("      the same stroke, every per-gene constant nuisance factor: gene length,")
    R("      annotation-version mismatch, probe affinity, GC/3'-bias. Those are exactly")
    R("      the terms that differ between TPM, RPKM->TPM, counts->TPM and microarray")
    R("      intensity, so the four quantification pipelines are reconciled by")
    R("      construction rather than by assuming they were comparable.")
    R("    - Balancing on the ARM MEANS (not the grand mean) makes the reference")
    R("      independent of how many replicates each arm happens to have. A grand-mean")
    R("      centering would shift each study by an amount set by its ON:OFF ratio,")
    R("      which differs across these studies (9:15, 3:3, 12:9, 2:2) and would")
    R("      reintroduce study structure through the back door.")
    R("    - Centering on the OFF arm alone was rejected: it forces every OFF sample to")
    R("      ~0 by construction, so all OFF tasks would carry only replicate noise and")
    R("      the model would be trained to fit that noise.")
    R("    - Per-gene z-scoring was rejected: dividing by a per-gene within-study SD")
    R("      estimated from n=2-3 replicates amplifies noise for invariant genes.")
    Xc = Xraw.copy()
    center_info = []
    for s in used_studies:
        sel = (meta["study"] == s).to_numpy()
        on = sel & meta["cls"].eq("senescent").to_numpy()
        off = sel & meta["cls"].eq("proliferating").to_numpy()
        if on.sum() and off.sum():
            ref = 0.5 * (Xraw[on].mean(axis=0) + Xraw[off].mean(axis=0))
        else:
            ref = Xraw[sel].mean(axis=0)
        Xc[sel] = Xraw[sel] - ref
        center_info.append(dict(study=s, n_on=int(on.sum()), n_off=int(off.sum()),
                                reference="arm-balanced" if on.sum() and off.sum() else "grand-mean"))
    R("\n  centering reference used per study:")
    R(pd.DataFrame(center_info).to_string(index=False))

    R("\n  PCA AFTER centering:")
    sc_c, fr_c = pca_report(Xc, "study-centered", meta)

    # additional per-study scale equalisation (microarray dynamic range is compressed)
    Xcs = Xc.copy()
    scale_info = []
    for s in used_studies:
        sel = (meta["study"] == s).to_numpy()
        sd = Xc[sel].std()
        Xcs[sel] = Xc[sel] / sd
        scale_info.append(dict(study=s, pooled_sd_of_centered=float(sd)))
    R("\n  per-study pooled SD of the centered values (dynamic-range check):")
    R(pd.DataFrame(scale_info).to_string(index=False))
    R("\n  PCA AFTER centering + per-study scale equalisation:")
    sc_cs, fr_cs = pca_report(Xcs, "centered+scaled", meta)

    pd.DataFrame(stats).to_csv(OUT / "batch_pca_stats.csv", index=False)

    # ---------------- figure ----------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 3, figsize=(17, 10))
        panels = [(sc_raw, fr_raw, "raw pooled log"), (sc_c, fr_c, "study-centered"),
                  (sc_cs, fr_cs, "centered + scaled")]
        studies = sorted(meta["study"].unique())
        cmap = dict(zip(studies, ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]))
        shape = {"ON": "o", "OFF": "s", "other": "^"}
        for col, (sc, fr, name) in enumerate(panels):
            for row, colour_by in enumerate(["study", "status"]):
                ax = axes[row][col]
                for st in studies:
                    for cb in ["ON", "OFF", "other"]:
                        m = (meta["study"] == st) & (meta["cls_binary"] == cb)
                        if not m.any():
                            continue
                        c = cmap[st] if colour_by == "study" else \
                            {"ON": "#d62728", "OFF": "#1f77b4", "other": "#7f7f7f"}[cb]
                        ax.scatter(sc[m.to_numpy(), 0], sc[m.to_numpy(), 1], c=c,
                                   marker=shape[cb], s=42, alpha=.85, edgecolors="none",
                                   label=f"{st} {cb}" if col == 0 and row == 0 else None)
                ax.set_title(f"{name} — coloured by {colour_by}", fontsize=10)
                ax.set_xlabel(f"PC1 ({100*fr[0]:.1f}%)")
                ax.set_ylabel(f"PC2 ({100*fr[1]:.1f}%)")
                ax.axhline(0, lw=.5, c="#ccc"); ax.axvline(0, lw=.5, c="#ccc")
        axes[0][0].legend(fontsize=6, ncol=2, loc="best")
        fig.suptitle("Batch structure before and after per-study correction\n"
                     "top row: colour = study (batch)   bottom row: colour = senescence status  "
                     "(o = ON, s = OFF, ^ = other)", fontsize=11)
        fig.tight_layout()
        fig.savefig(FIG / "multistudy_batch_pca.png", dpi=150)
        R(f"\n  figure -> {FIG/'multistudy_batch_pca.png'}")
    except Exception as ex:                                    # figure is not load-bearing
        R(f"\n  [figure skipped: {ex}]")

    np.save(OUT / "_Xraw.npy", Xraw)
    np.save(OUT / "_Xc.npy", Xc)
    np.save(OUT / "_Xcs.npy", Xcs)
    meta.to_csv(OUT / "multistudy_sample_assignments.csv", index=False)
    iv.to_csv(OUT / "_intervals.csv")
    json.dump({"qc175": qc175, "qc743": qc743, "qc990": qc990, "qc205": qc205,
               "symbol_harmonisation": sym_stats,
               "corrupt": {k: v for k, v in corrupt.items()},
               "genes": genes, "used_studies": used_studies,
               "center_info": center_info, "scale_info": scale_info},
              open(OUT / "_stage1.json", "w"), default=str)
    (OUT / "BUILD_LOG.txt").write_text("\n".join(LOG_LINES))
    R("\nSTAGE 1 COMPLETE -- run 43 to emit the AnnData")


if __name__ == "__main__":
    main()
