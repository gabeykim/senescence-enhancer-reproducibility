#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# verify_figure_typography.py
#
# Reads the ACTUAL rendered font sizes out of each figure PDF rather than
# trusting the source, by decompressing the content streams and parsing the
# `/F<n> <size> Tf` text-font operators.
#
# CONSUMES: figures/*.pdf
# PRODUCES: stdout report; non-zero exit if anything is below the floor
# ---------------------------------------------------------------------------
import re, sys, zlib
from pathlib import Path

FLOOR = 9.0
FIG = Path(__file__).resolve().parents[1] / "figures"
STREAM = re.compile(rb"stream\r?\n(.*?)endstream", re.S)
TF = re.compile(rb"/F\d+\s+([\d.]+)\s+Tf")

bad = 0
print(f"{'figure':46s} {'min pt':>7s} {'sizes used':>34s}  {'mm':>4s}")
for pdf in sorted(FIG.glob("*.pdf")):
    raw = pdf.read_bytes()
    sizes = set()
    for m in STREAM.finditer(raw):
        blob = m.group(1)
        try:
            blob = zlib.decompress(blob)
        except zlib.error:
            pass
        sizes.update(round(float(x), 2) for x in TF.findall(blob))
    png = pdf.with_suffix(".png")
    mm = ""
    if png.exists():
        from PIL import Image
        mm = f"{Image.open(png).size[0] / 300 * 25.4:.0f}"
    if not sizes:
        print(f"{pdf.name:46s} {'(no text found)':>7s}"); continue
    primary = {x for x in sizes if x >= FLOOR - 1e-6}
    # mathtext renders superscripts/subscripts at 0.7x the surrounding size;
    # those are correct to be smaller and are not a legibility violation.
    def exempt(x):
        return any(abs(x - 0.7 * p) < 0.06 for p in primary)
    viol = sorted(x for x in sizes if x < FLOOR - 1e-6 and not exempt(x))
    lo = min(primary) if primary else min(sizes)
    flag = "" if not viol else f"  <-- BELOW FLOOR: {viol}"
    if viol:
        bad += 1
    print(f"{pdf.name:46s} {lo:7.2f} {str(sorted(sizes)):>34s}  {mm:>4s}{flag}")
print(f"\nfloor {FLOOR} pt (mathtext super/subscripts at 0.7x exempt) — "
      f"{'ALL PASS' if bad == 0 else f'{bad} FIGURE(S) BELOW FLOOR'}")
sys.exit(1 if bad else 0)
