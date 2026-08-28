# Reference files

This directory is reserved for external reference resources required to run the
RNA-seq preprocessing and differential-expression workflows.

## Required reference resources

| Resource | Used by | Purpose |
|---|---|---|
| HISAT2 GRCm38/mm10 index | `RNAseq_clean.R` | Paired-end alignment of trimmed reads. |
| `Mus_musculus.GRCm38.102.gtf` | `DGE_clean.R` | Gene-level quantification with `featureCounts`. |
| `TruSeq3-PE-2.fa` | `RNAseq_clean.R` | Adapter trimming with `fastp`. |

The alignment reference and the GTF annotation must refer to the same genome
build. The original analysis used the mouse GRCm38/mm10 reference and the
Ensembl GRCm38.102 GTF annotation.

## Expected layout

```text
reference/
├── README.md
├── TruSeq3-PE-2.fa
├── Mus_musculus.GRCm38.102.gtf
└── GRCm38_hisat2_index/
    ├── grcm38.1.ht2
    ├── grcm38.2.ht2
    └── ...
```

Reference files are commonly too large to commit directly to Git. Store them
locally or provide them through an institutional archive, then set their paths
in `config.R` (created from `config_template.R`).

## Reproducibility note

Changing the genome build, GTF release, or HISAT2 index can change alignment,
counting, annotation, filtering, and downstream differential-expression
results. Use the GRCm38/mm10 resources specified above when reproducing the
original analysis.
