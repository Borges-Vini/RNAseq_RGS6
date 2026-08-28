# =============================================================================
# Configuration template for the RGS6 RNA-seq repository
# =============================================================================
# Copy this file to config.R, replace the placeholder paths, and keep config.R
# out of version control (it may contain machine-specific paths). The analysis
# scripts retain their defaults; copy the relevant values below into Section 1
# of RNAseq_clean.R and DGE_clean.R before running.

# RNAseq_clean.R ---------------------------------------------------------------
rnaseq_config <- list(
  input_dir = "path/to/raw_fastq",
  output_dir = "path/to/rnaseq_results",
  hisat2_index = "path/to/GRCm38_hisat2_index",
  adapter_fasta = "path/to/TruSeq3-PE-2.fa",
  conda_env = "RNAseq_env",
  run_multiqc = TRUE
)

# Parameters used in the original paired-end preprocessing workflow.
rnaseq_parameters <- list(
  fastp_window_size = 4L,
  fastp_mean_quality = 15L,
  fastp_minimum_length = 25L,
  fastp_threads = as.integer(parallel::detectCores()),
  hisat2_threads = as.integer(parallel::detectCores() * 0.95),
  fastqc_threads = as.integer(parallel::detectCores() * 2 / 3)
)

# DGE_clean.R ------------------------------------------------------------------
dge_config <- list(
  bam_dir = "path/to/rnaseq_results/Aligned",
  results_dir = "path/to/dge_results",
  metadata_file = "path/to/dge_results/Rorabaugh_FASTQ_Sample_Guide.csv",
  gtf_file = "path/to/Mus_musculus.GRCm38.102.gtf",

  # Keep NULL to reproduce the original GSEA behavior. The historical script
  # did not set a seed before its random tie-breaking step.
  gsea_seed = NULL
)

# Required metadata schema for DGE_clean.R.
# Mouse_ID: unique biological-sample identifier used after counting.
# Sample_code: BAM filename stem (for example, SAMPLE_01 for SAMPLE_01.bam).
# Sex: exactly "Female" or "Male".
# Mouseline: exactly "RGS6KO" or "WT".
required_metadata_columns <- c("Mouse_ID", "Sample_code", "Sex", "Mouseline")
