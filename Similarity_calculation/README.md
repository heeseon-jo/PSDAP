# Similarity calculation

Pairwise similarity between signatures based on GOGO scores for the three Gene
Ontology branches (BPO, CCO and MFO). The scripts prepare the signature pairs
that are given to GOGO and merge the GOGO results of the three branches into one
table.

| Script | Description |
|---|---|
| `make_pairs_100k.py` | Generates the signature pair files for GOGO: every CMAP signature x every OTHER signature, followed by all unordered OTHER x OTHER pairs. Pairs are split into files of 100,000 pairs (default), and the pair count is checked against the expected number. |
| `make_one_gogo_wide_chunk.py` | Merges the GOGO results of BPO, CCO and MFO into one table with a row per signature pair. Pair IDs are checked across the three files, and mismatches are counted in a summary file. |

## Usage
### make_pairs_100k.py
```
python make_pairs_100k.py \
    --cmap cmap_gogo_input.txt \
    --new new_gogo_input.txt \
    --out-dir pairs \
    --prefix pairs \
    --chunk-size 100000
```
`--cmap(CMAP)` and `--new(OTHER)` are GOGO input files (one signature per line, signature ID
first); signatures with duplicated IDs are removed. The pair files
`<prefix>_sub_<number>.pair` are written to `--out-dir` together with
`manifest_100k.tsv` (file list) and `pair_generation_summary.tsv`. Existing pair
files with the same prefix and the two summary files in `--out-dir` are deleted
before writing.

### make_one_gogo_wide_chunk.py
```
python make_one_gogo_wide_chunk.py \
    --base pairs_sub_000000 \
    --raw-dir gogo_results \
    --out-dir wide \
    --strict
```

`--raw-dir` must contain `<base>_BPO.txt`, `<base>_CCO.txt` and `<base>_MFO.txt`
(or the `.gz` versions) with the same pairs in the same order. The output is
`<base>.wide.tsv.gz` (columns `sigA`, `sigB`, `GOGO_BPO`, `GOGO_CCO`,
`GOGO_MFO`) and `<base>.wide.summary.tsv`. With `--strict`, the run stops at the
first pair-ID mismatch instead of only counting it.

## Requirements

Python 3 (standard library only).
