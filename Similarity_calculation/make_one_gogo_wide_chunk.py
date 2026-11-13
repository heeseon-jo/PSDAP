# =============================================================================
# merge_gogo_wide_chunk.py
#
# Merges per-ontology GOGO results (BPO, CCO, MFO) into one wide table with a
# row per signature pair.
#
# The three input files must list the same signature pairs in the same order.
# Pair IDs are checked across files; mismatches are counted in a summary file
# (or abort the run with --strict). Deterministic; no random seed required.
#
# Dependencies: Python 3 (standard library only).
# =============================================================================


import argparse
import gzip
import re
from pathlib import Path
from itertools import zip_longest



# Marks where the GO term list starts in a pair string.
_GO_RE = re.compile(r'\bGO:')




# Open a text file for reading; gzip-compressed files (.gz) are handled
# transparently.
# @param path  Path to a plain-text or .gz file.
# @return      File object opened in text mode.
def open_in(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")



# Extract the signature ID from one side of a pair.
# The ID is everything before the first "GO:" token
# @param s  One side of a pair string (signature ID followed by GO terms).
# @return   Signature ID as a string.
def extract_sig_id(s):
    s = s.strip()
    m = _GO_RE.search(s)
    if m:
        return s[:m.start()].strip()
    return s.strip()



# Parse one GOGO output line: '<sigA ...>;<sigB ...> <ontology> <score>'.
# The score is kept as a string to preserve the original precision.
# @param line  One line of a GOGO result file.
# @return      Tuple (sig_a, sig_b, ontology, score).
def parse_line(line):
    line = line.rstrip("\n")
    pair_part, ont, score = line.rsplit(None, 2)

    if ";" not in pair_part:
        raise ValueError(f"No semicolon in line: {line[:300]}")

    left, right = pair_part.split(";", 1)

    sig_a = extract_sig_id(left)
    sig_b = extract_sig_id(right)

    return sig_a, sig_b, ont, score



# Merge the BPO, CCO and MFO files of one base name line by line.
# Writes <base>.wide.tsv.gz (sigA, sigB, GOGO_BPO, GOGO_CCO, GOGO_MFO) and
# <base>.wide.summary.tsv (pair count, number of pair-ID mismatches).
# Command-line arguments:
#   --base     Base name of the input files (<base>_{BPO,CCO,MFO}.txt[.gz]).
#   --raw-dir  Directory containing the input files.
#   --out-dir  Output directory (created if missing).
#   --strict   Abort on a pair-ID mismatch instead of only counting it.
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = args.base

    bpo = raw_dir / f"{base}_BPO.txt"
    cco = raw_dir / f"{base}_CCO.txt"
    mfo = raw_dir / f"{base}_MFO.txt"

    if not bpo.exists():
        bpo = raw_dir / f"{base}_BPO.txt.gz"
        cco = raw_dir / f"{base}_CCO.txt.gz"
        mfo = raw_dir / f"{base}_MFO.txt.gz"

    if not bpo.exists():
        raise RuntimeError(f"missing BPO: {bpo}")
    if not cco.exists():
        raise RuntimeError(f"missing CCO: {cco}")
    if not mfo.exists():
        raise RuntimeError(f"missing MFO: {mfo}")

    out_file = out_dir / f"{base}.wide.tsv.gz"
    summary_file = out_dir / f"{base}.wide.summary.tsv"

    n = 0
    mismatch = 0

    with open_in(bpo) as fb, open_in(cco) as fc, open_in(mfo) as fm, gzip.open(out_file, "wt") as out:
        out.write("sigA\tsigB\tGOGO_BPO\tGOGO_CCO\tGOGO_MFO\n")

        for line_no, (lb, lc, lm) in enumerate(zip_longest(fb, fc, fm), start=1):
            if lb is None or lc is None or lm is None:
                raise RuntimeError(f"line count mismatch: {base} line {line_no}")

            a1, b1, o1, s1 = parse_line(lb)
            a2, b2, o2, s2 = parse_line(lc)
            a3, b3, o3, s3 = parse_line(lm)

            if (o1, o2, o3) != ("BPO", "CCO", "MFO"):
                raise RuntimeError(
                    f"ontology mismatch: {base} line {line_no}: {o1},{o2},{o3}"
                )

            if not (a1 == a2 == a3 and b1 == b2 == b3):
                mismatch += 1
                if args.strict:
                    raise RuntimeError(f"pair mismatch: {base} line {line_no}")

            out.write(f"{a1}\t{b1}\t{s1}\t{s2}\t{s3}\n")
            n += 1

    with open(summary_file, "w") as s:
        s.write("base\twide_file\tn_pairs\tn_pair_mismatch\n")
        s.write(f"{base}\t{out_file}\t{n}\t{mismatch}\n")

    print(f"DONE\t{base}\t{n}\t{mismatch}", flush=True)

if __name__ == "__main__":
    main()
