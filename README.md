# RGS6 knockout bulk RNA-seq analysis

[![DOI](https://zenodo.org/badge/1349394142.svg)](https://doi.org/10.5281/zenodo.22142601)

Reproducible preprocessing and differential-expression analysis for the RGS6
knockout bulk RNA-seq study.

## Workflow

```text
Paired-end FASTQ
  → FastQC / MultiQC
  → fastp trimming
  → HISAT2 alignment
  → SAMtools duplicate marking and indexing
  → featureCounts
  → edgeR TMM normalization and filtering
  → limma-voom differential expression
  → GSEA, Table 1, Figure 6, and GEO supplementary files
```

## Repository layout

```text
.
├── README.md
├── config_template.R
├── RNAseq_clean.R
├── DGE_clean.R
└── reference/
    └── README.md
```

## Requirements

### Software

- R and the packages checked by each analysis script
- Conda environment containing `fastp`, `hisat2`, `samtools`, and `multiqc`
- Java, required by FastQC

### Reference resources

The original workflow used mouse GRCm38/mm10 and the Ensembl GRCm38.102 GTF.
See [`reference/README.md`](reference/README.md) for the required files and
expected directory layout.

## Configuration

1. Copy `config_template.R` to `config.R`.
2. Replace every `path/to/...` placeholder with local paths.
3. Copy the relevant configuration values into Section 1 of
   `RNAseq_clean.R` and `DGE_clean.R`.
4. Do not commit `config.R` if it contains machine-specific paths.

The differential-expression metadata must contain these columns:

| Column | Description |
|---|---|
| `Mouse_ID` | Unique biological-sample identifier. |
| `Sample_code` | BAM filename stem; e.g., `sample_01` for `sample_01.bam`. |
| `Sex` | `Female` or `Male`. |
| `Mouseline` | `RGS6KO` or `WT`. |

The metadata row order is meaningful: the DGE script creates the BAM list in
`Sample_code` order, exactly as in the original analysis.

## Run order

Run the preprocessing script first:

```r
source("RNAseq_clean.R")
```

Then point `bam_dir` in `DGE_clean.R` to the resulting `Aligned/` directory
and run:

```r
source("DGE_clean.R")
```

## Analysis methods retained from the original workflow

- paired-end `featureCounts` quantification using GRCm38.102 annotation;
- edgeR `DGEList`, TMM normalization, and `filterByExpr` filtering;
- `voomWithQualityWeights` and a `~ 0 + Sex * Mouseline` model;
- limma empirical-Bayes contrasts, including the overall WT vs RGS6KO effect;
- GO Biological Process and MSigDB Hallmark GSEA;
- Table 1, Figure 6, and GEO raw-count, voom log2CPM, and metadata exports.

## Reproducibility notes

Use the same FASTQ files, BAM-generation tools, reference build, annotation,
and package versions to reproduce the original results. The original GSEA code
uses random noise to break tied ranking scores but did not record a seed;
therefore exact GSEA tie resolution from the historical run cannot be
reconstructed solely from the scripts. Set `gsea_seed` in `DGE_clean.R` for
deterministic new runs.
