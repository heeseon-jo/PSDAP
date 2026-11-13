# =============================================================================
# make_pairs_100k.py
#
# Generates signature pair files for GOGO, split into files of a fixed number
# of pairs (default 100,000): every CMAP signature x every NEW signature,
# followed by all unordered NEW x NEW pairs (each pair once).
#
# Input signatures are read from GOGO input files (one signature per line,
# signature ID first). Duplicated IDs are removed, and the number of pairs
# written is checked against the expected number. Existing pair files with
# the same prefix and the manifest/summary in the output directory are
# removed first. Deterministic; no random seed required.
#
# Dependencies: Python 3 (standard library only).
# =============================================================================

import argparse
from pathlib import Path
import sys
import time


# Read a GOGO input file, removing signatures with duplicated IDs.
# @param path  Path to a GOGO input file.
# @return      Tuple (rows, n_duplicates): kept lines and number of removed
#              duplicates.
def read_gogo_input(path):
    rows = []
    seen = set()
    dup = 0

    path = Path(path)

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                raise ValueError(
                    f"Bad line with <2 tokens in {path}, line {line_no}: {line}"
                )

            sig = parts[0]

            if sig in seen:
                dup += 1
                continue

            seen.add(sig)
            rows.append(line)

    return rows, dup


# Write pairs to numbered files of at most `chunk_size` pairs and keep a
# manifest of the files.
# @param out_dir     Output directory.
# @param prefix      Output file prefix.
# @param chunk_size  Maximum number of pairs per file.
class ChunkWriter:
    def __init__(self, out_dir, prefix, chunk_size):
        self.out_dir = Path(out_dir)
        self.prefix = prefix
        self.chunk_size = chunk_size

        self.sub_idx = 0
        self.current_n = 0
        self.current_fh = None
        self.current_path = None

        self.manifest_rows = []
        self.total_pairs = 0

    def _open_new(self, pair_type):
        if self.current_fh is not None:
            self.current_fh.close()
            self.manifest_rows[-1]["n_pair"] = self.current_n

        sub_id = f"{self.sub_idx:06d}"
        self.current_path = self.out_dir / f"{self.prefix}_sub_{sub_id}.pair"
        self.current_fh = self.current_path.open("w", encoding="utf-8")
        self.current_n = 0

        self.manifest_rows.append({
            "sub_id": sub_id,
            "pair_type": pair_type,
            "file": str(self.current_path),
            "n_pair": 0,
        })

        self.sub_idx += 1

    def write_pair(self, left, right, pair_type):
        if self.current_fh is None or self.current_n >= self.chunk_size:
            self._open_new(pair_type)

        self.current_fh.write(left + "; " + right + "\n")
        self.current_n += 1
        self.total_pairs += 1

    def close(self):
        if self.current_fh is not None:
            self.current_fh.close()
            self.manifest_rows[-1]["n_pair"] = self.current_n
            self.current_fh = None


# Generate all pair files, check the pair count, and write the manifest and
# run summary.
def main():
    ap = argparse.ArgumentParser(
        description="Create CMAP_NEW + NEW_NEW pair files split into 100k subchunks."
    )
    ap.add_argument("--cmap", required=True, help="CMAP gogo input file")
    ap.add_argument("--new", required=True, help="NEW gogo input file")
    ap.add_argument("--out-dir", required=True, help="Output pair directory")
    ap.add_argument("--prefix", default="pairs", help="Output pair file prefix")
    ap.add_argument("--chunk-size", type=int, default=100000)
    args = ap.parse_args()

    t0 = time.time()

    cmap_path = Path(args.cmap)
    new_path = Path(args.new)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not cmap_path.exists():
        sys.exit(f"ERROR: CMAP input not found: {cmap_path}")
    if not new_path.exists():
        sys.exit(f"ERROR: NEW input not found: {new_path}")

    print("============================================================", flush=True)
    print("make_pairs_100k.py", flush=True)
    print("CMAP       :", cmap_path, flush=True)
    print("NEW        :", new_path, flush=True)
    print("OUT_DIR    :", out_dir, flush=True)
    print("PREFIX     :", args.prefix, flush=True)
    print("CHUNK_SIZE :", args.chunk_size, flush=True)
    print("============================================================", flush=True)

    cmap_rows, cmap_dup = read_gogo_input(cmap_path)
    new_rows, new_dup = read_gogo_input(new_path)

    n_cmap = len(cmap_rows)
    n_new = len(new_rows)

    expected_cmap_new = n_cmap * n_new
    expected_new_new = n_new * (n_new - 1) // 2
    expected_total = expected_cmap_new + expected_new_new

    print(f"N_CMAP={n_cmap}", flush=True)
    print(f"N_NEW={n_new}", flush=True)
    print(f"DUP_CMAP_REMOVED={cmap_dup}", flush=True)
    print(f"DUP_NEW_REMOVED={new_dup}", flush=True)
    print(f"EXPECTED_CMAP_NEW={expected_cmap_new}", flush=True)
    print(f"EXPECTED_NEW_NEW={expected_new_new}", flush=True)
    print(f"EXPECTED_TOTAL={expected_total}", flush=True)

    manifest_path = out_dir / "manifest_100k.tsv"
    summary_path = out_dir / "pair_generation_summary.tsv"

    for old in out_dir.glob(f"{args.prefix}_sub_*.pair"):
        old.unlink()

    for old in [manifest_path, summary_path]:
        if old.exists():
            old.unlink()

    writer = ChunkWriter(
        out_dir=out_dir,
        prefix=args.prefix,
        chunk_size=args.chunk_size
    )

    print("START CMAP_NEW:", time.ctime(), flush=True)
    for c in cmap_rows:
        for n in new_rows:
            writer.write_pair(c, n, "CMAP_NEW")
    print("END CMAP_NEW:", time.ctime(), flush=True)

    print("START NEW_NEW:", time.ctime(), flush=True)
    for i in range(n_new):
        left = new_rows[i]
        for j in range(i + 1, n_new):
            writer.write_pair(left, new_rows[j], "NEW_NEW")
    print("END NEW_NEW:", time.ctime(), flush=True)

    writer.close()

    actual_total = writer.total_pairs

    if actual_total != expected_total:
        sys.exit(
            f"ERROR: pair count mismatch actual={actual_total} expected={expected_total}"
        )

    with manifest_path.open("w", encoding="utf-8") as f:
        f.write("sub_id\tpair_type\tfile\tn_pair\n")
        for row in writer.manifest_rows:
            f.write(
                f"{row['sub_id']}\t{row['pair_type']}\t{row['file']}\t{row['n_pair']}\n"
            )

    with summary_path.open("w", encoding="utf-8") as f:
        f.write("key\tvalue\n")
        f.write(f"cmap_input\t{cmap_path}\n")
        f.write(f"new_input\t{new_path}\n")
        f.write(f"out_dir\t{out_dir}\n")
        f.write(f"prefix\t{args.prefix}\n")
        f.write(f"chunk_size\t{args.chunk_size}\n")
        f.write(f"n_cmap\t{n_cmap}\n")
        f.write(f"n_new\t{n_new}\n")
        f.write(f"dup_cmap_removed\t{cmap_dup}\n")
        f.write(f"dup_new_removed\t{new_dup}\n")
        f.write(f"expected_cmap_new\t{expected_cmap_new}\n")
        f.write(f"expected_new_new\t{expected_new_new}\n")
        f.write(f"expected_total\t{expected_total}\n")
        f.write(f"actual_total\t{actual_total}\n")
        f.write(f"n_subchunks\t{len(writer.manifest_rows)}\n")
        f.write(f"elapsed_sec\t{time.time() - t0:.2f}\n")

    print("============================================================", flush=True)
    print("DONE:", time.ctime(), flush=True)
    print(f"ACTUAL_TOTAL={actual_total}", flush=True)
    print(f"N_SUBCHUNKS={len(writer.manifest_rows)}", flush=True)
    print("MANIFEST:", manifest_path, flush=True)
    print("SUMMARY :", summary_path, flush=True)
    print("============================================================", flush=True)


if __name__ == "__main__":
    main()