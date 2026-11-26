# RNA-seq Batch Effect Correction with pyComBat-seq (Docker)

Batch effect correction tool for RNA-seq count data (htseq-count output), based
on [pyComBat-seq](https://inmoose.readthedocs.io/en/stable/pycombatseq.html)
from the inmoose package. The script is fully parameterized through a **design
file**: no sample is hard-coded and the same container can be reused for any
experiment.

## Directory contents

| File | Purpose |
|---|---|
| `Dockerfile` | Python image with inmoose (pyComBat-seq) and pandas |
| `batch_correction.py` | Correction script (bundled into the image) |
| `design_<Experiment>.csv` | Design table: one sample per row |
| `<sample>.tsv` | htseq-count output: 2 columns `gene` / `count`, no header |

## 1. Building the image

```bash
docker build -t xgrand/pycombatseq:1.0.0 -f Dockerfile.txt .
```

`batch_correction.py` must sit in the same folder as the `Dockerfile` at build
time. The script is copied to `/app` (outside `/data`).

## 2. Design file

One row per sample, CSV or TSV (separator auto-detected). Example
(`design_Huh75.csv`):

| prefix | sample | condition | repetition |
|---|---|---|---|
| 1-7-A2-NP | 1-7-A2-NP | A2_NP | 1 |
| 2-7-A2-Poly | 2-7-A2-Poly | A2_Poly | 2 |
| ... | ... | ... | ... |

- **sample**: count file name without extension (`1-7-A2-NP.tsv` here);
- **repetition**: the batch variable to correct for;
- **condition**: the biological covariate to preserve;
- any identifier-like column (one unique value per sample, e.g. `prefix`) is
  automatically excluded from the covariates.

Constraint: the batch must not be confounded with a covariate (each batch
should contain several levels of every covariate); otherwise the script prints
a warning.

## 3. Usage

```bash
docker run --rm -v "$(pwd)":/data xgrand/pycombatseq:1.0.0 \
    python /app/batch_correction.py \
    --design ./design_Huh75.csv \
    --counts-dir ./09_htseq_count \
    --batch-col repetition \
    --covar-cols condition \
    --out-matrix 11_combatseq/Huh75_counts_corrected.csv \
    --htseq-dir 11_combatseq/Huh75_htseq
```

Count files and design are mounted through the `/data` volume. Relative paths
are resolved against this folder inside the container.

### Options

| Option | Default | Description |
|---|---|---|
| `--design` | (required) | Design file (one row per sample) |
| `--counts-dir` | `.` | Directory containing the htseq-count files |
| `--sample-col` | `sample` | Design column with sample names |
| `--batch-col` | `repetition` | Design column defining the batches |
| `--covar-cols` | (auto) | Biological covariates to preserve; default: all columns except sample and batch (identifier-like columns excluded) |
| `--file-ext` | `.tsv` | Extension of the count files |
| `--counts-sep` | tab | Separator of the count files |
| `--design-sep` | auto | Design separator (auto-detected by default) |
| `--min-count` | `10` | Minimum count for a gene to be considered expressed |
| `--min-samples` | `3` | Minimum number of samples in which a gene must be expressed |
| `--out-matrix` | `counts_corrected.tsv` | Corrected genes × samples matrix (TSV) |
| `--out-raw` | (disabled) | Raw counts after the expression filter |
| `--out-metadata` | (disabled) | Design with the batch column renamed (DESeq2 colData) |
| `--htseq-dir` | `counts_corrected_htseq` | Output folder for htseq-format files, one per sample; empty string to disable |
| `--no-htseq-stats` | off | Do not copy the `__` summary lines into the per-sample files |

Full help is available with:

```bash
docker run --rm -v "$(pwd)":/data xgrand/pycombatseq:1.0.0 \
    python /app/batch_correction.py --help
```

## 4. Outputs

- **`--out-matrix`**: corrected count matrix (genes × samples, integers),
  suitable for data exploration (PCA, clustering).
- **`--htseq-dir`**: one file per sample, in the exact htseq-count output
  format (2 columns, no header, with the `__no_feature`, `__ambiguous`...
  summary lines copied back at the end of the file), directly consumable by a
  DESeq2 pipeline. All files share the same gene set.

The script prints the matrix size, the batch composition and the preserved
covariates to the console.

## 5. Quality checks and best practices

- **Check the result** with a PCA on the corrected counts: samples should now
  cluster by condition rather than by batch.
- **DESeq2**: for differential expression, the statistical gold standard is
  raw counts with the batch included in the model
  (`design = ~ batch + condition`). Corrected counts can be fed to DESeq2, but
  p-values are no longer strictly controlled; prefer them for exploration,
  clustering or gene prioritization.
- The expression filter (`--min-count`, `--min-samples`) prevents near-zero
  genes from destabilizing the ComBat-seq parameter estimation. Adjust it to
  your sample count.

## 6. Technical notes

- Compatible with recent inmoose versions (`covar_mod`) and older ones
  (`covar`); a gene missing from an htseq output is set to 0; htseq-count
  summary lines are excluded from the correction and copied back verbatim.
- Outputs are written inside the mounted folder (`/data`), so they appear
  directly on the host machine.