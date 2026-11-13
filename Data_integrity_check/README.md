# Data integrity check

Automated integrity check of the four released metadata files
(`signature_info.csv`, `perturbagen_info.csv`, `model_info.csv` and
`source_info.csv`).

| Script | Description |
|---|---|
| `check_metadata_integrity.R` | Checks the integrity of the four metadata files and reports the number of problem records of each check. |

## Usage

```
Rscript check_metadata_integrity.R
```

The four metadata files are read from `output/Metadata` (set in the last line of
the script). The number of problem records of each check is printed, followed by
`Integrity check: PASS` when all checks report 0 problems and
`Integrity check: FAIL` otherwise.

### Checks

- Unnamed columns, duplicated rows, and empty or duplicated primary IDs in each file.
- Every `Pertur_id`, `Model_id` and `Source_id` used in `signature_info.csv` is defined in `perturbagen_info.csv`, `model_info.csv` and `source_info.csv`, respectively.
- The `Source` name of each signature in `signature_info.csv` matches the name of its `Source_id` in `source_info.csv`.
- Perturbagens, models and sources that no signature refers to (reported as notes, not as errors).

## Requirements

R (v4.4.1); base R only.
