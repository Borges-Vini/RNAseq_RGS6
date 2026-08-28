# =============================================================================
# Paired-end RNA-seq preprocessing
# =============================================================================
# Faithful clean version of RNAseq.R for the paired-end dataset used by the
# downstream DGE analysis. The operative fastp, HISAT2, and SAMtools parameters
# are retained; only configuration, validation, and file handling are organized.
# Input files must be named SAMPLE_R1.fastq.gz and SAMPLE_R2.fastq.gz.
# =============================================================================

# 1. Configuration --------------------------------------------------------------
input_dir <- "path/to/raw_fastq"
output_dir <- "path/to/rnaseq_results"
hisat2_index <- "path/to/GRCm38_hisat2_index"
adapter_fasta <- "path/to/TruSeq3-PE-2.fa"
conda_env <- "RNAseq_env"
run_multiqc <- TRUE

# Original parameters
fastp_window_size <- 4L
fastp_mean_quality <- 15L
fastp_minimum_length <- 25L
fastp_threads <- as.integer(parallel::detectCores())
hisat2_threads <- as.integer(parallel::detectCores() * 0.95)
fastqc_threads <- as.integer(parallel::detectCores() * 2 / 3)

# 2. Requirements and folders --------------------------------------------------
if (!requireNamespace("fastqcr", quietly = TRUE)) stop("Install the R package 'fastqcr' before running this script.")
if (!nzchar(Sys.which("conda"))) stop("The conda executable is not available on PATH.")
for (path in c(input_dir, hisat2_index, adapter_fasta)) if (!file.exists(path)) stop("Required path does not exist: ", path)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
qc_dir <- file.path(output_dir, "QC_results")
trimmed_dir <- file.path(output_dir, "Trimmed")
unpaired_dir <- file.path(trimmed_dir, "Unpaired")
aligned_dir <- file.path(output_dir, "Aligned")
logs_dir <- file.path(aligned_dir, "Logs")
for (dir in c(qc_dir, trimmed_dir, unpaired_dir, aligned_dir, logs_dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

# 3. Discover and validate paired-end input ------------------------------------
r1_files <- sort(list.files(input_dir, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE))
if (!length(r1_files)) stop("No *_R1.fastq.gz files found in: ", input_dir)
sample_ids <- sub("_R1\\.fastq\\.gz$", "", basename(r1_files))
r2_files <- file.path(dirname(r1_files), paste0(sample_ids, "_R2.fastq.gz"))
if (any(!file.exists(r2_files))) stop("Missing R2 file(s) for: ", paste(sample_ids[!file.exists(r2_files)], collapse = ", "))
if (anyDuplicated(sample_ids)) stop("Sample identifiers are not unique.")
fastq_pairs <- data.frame(sample = sample_ids, r1 = r1_files, r2 = r2_files, stringsAsFactors = FALSE)

# 4. Raw-read QC ---------------------------------------------------------------
fastqcr::fastqc(fq.dir = input_dir, qc.dir = qc_dir, threads = fastqc_threads)
if (run_multiqc) {
  status <- system2("conda", c("run", "-n", conda_env, "multiqc", qc_dir, "-o", qc_dir))
  if (status != 0) warning("MultiQC returned a non-zero exit status.")
}

run_checked <- function(command, args, label) {
  status <- system2(command, args)
  if (status != 0) stop(label, " failed (exit status ", status, ").")
}

# 5. Adapter and quality trimming ----------------------------------------------
# Original fastp options: --adapter_fasta; --cut_front; --cut_tail;
# --cut_window_size 4; --cut_mean_quality 15; --length_required 25.
for (i in seq_len(nrow(fastq_pairs))) {
  sample <- fastq_pairs$sample[i]
  run_checked("conda", c(
    "run", "-n", conda_env, "fastp",
    "--in1", fastq_pairs$r1[i], "--in2", fastq_pairs$r2[i],
    "--out1", file.path(trimmed_dir, paste0(sample, "_R1.fastq.gz")),
    "--out2", file.path(trimmed_dir, paste0(sample, "_R2.fastq.gz")),
    "--unpaired1", file.path(unpaired_dir, paste0(sample, "_R1.fastq.gz")),
    "--unpaired2", file.path(unpaired_dir, paste0(sample, "_R2.fastq.gz")),
    "--failed_out", file.path(trimmed_dir, paste0(sample, "_failed.out")),
    "--adapter_fasta", adapter_fasta, "--cut_front", "--cut_tail",
    "--cut_window_size", fastp_window_size, "--cut_mean_quality", fastp_mean_quality,
    "--length_required", fastp_minimum_length, "--thread", fastp_threads,
    "--html", file.path(trimmed_dir, paste0(sample, "_fastp_report.html")),
    "--json", file.path(trimmed_dir, paste0(sample, "_fastp_report.json"))
  ), paste("fastp for", sample))
}

# 6. Alignment and BAM processing ---------------------------------------------
# HISAT2 read-group fields and SAMtools order are preserved from the original.
for (sample in fastq_pairs$sample) {
  r1 <- file.path(trimmed_dir, paste0(sample, "_R1.fastq.gz"))
  r2 <- file.path(trimmed_dir, paste0(sample, "_R2.fastq.gz"))
  tmp_dir <- file.path(output_dir, "tmp", sample)
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  sam <- file.path(tmp_dir, paste0(sample, ".sam"))
  name_bam <- file.path(tmp_dir, paste0(sample, "_name_sorted.bam"))
  fixmate_bam <- file.path(tmp_dir, paste0(sample, "_fixmate.bam"))
  coord_bam <- file.path(tmp_dir, paste0(sample, "_coord_sorted.bam"))
  final_bam <- file.path(aligned_dir, paste0(sample, ".bam"))
  pipeline <- paste(
    "set -euo pipefail;",
    "hisat2 --dta --phred33", "-p", hisat2_threads,
    "--rg-id", shQuote(sample), "--rg", shQuote(paste0("SM:", sample)),
    "--rg", shQuote("CN:Genomics_Core"), "--rg", shQuote("PL:Illumina"), "--rg", shQuote("PM:NextSeq2000"),
    "-x", shQuote(hisat2_index), "-1", shQuote(r1), "-2", shQuote(r2), "-S", shQuote(sam),
    "2>", shQuote(file.path(logs_dir, paste0(sample, "_align.log"))), "&&",
    "samtools sort -n -@", hisat2_threads, "-o", shQuote(name_bam), shQuote(sam), "2>", shQuote(file.path(logs_dir, paste0(sample, "_name_sort.log"))), "&&",
    "samtools fixmate -m", shQuote(name_bam), shQuote(fixmate_bam), "2>", shQuote(file.path(logs_dir, paste0(sample, "_fixmate.log"))), "&&",
    "samtools sort -@", hisat2_threads, "-o", shQuote(coord_bam), shQuote(fixmate_bam), "2>", shQuote(file.path(logs_dir, paste0(sample, "_pos_sort.log"))), "&&",
    "samtools quickcheck -v", shQuote(coord_bam), "&&",
    "samtools markdup", shQuote(coord_bam), shQuote(final_bam), "2>", shQuote(file.path(logs_dir, paste0(sample, "_markdup.log"))), "&&",
    "samtools index", shQuote(final_bam), "2>", shQuote(file.path(logs_dir, paste0(sample, "_index.log")))
  )
  run_checked("conda", c("run", "-n", conda_env, "bash", "-c", pipeline), paste("alignment for", sample))
  unlink(tmp_dir, recursive = TRUE)
}

# 7. Alignment-rate summary ----------------------------------------------------
alignment_logs <- list.files(logs_dir, pattern = "_align\\.log$", full.names = TRUE)
alignment_rate <- vapply(alignment_logs, function(log_file) {
  lines <- readLines(log_file, warn = FALSE)
  if (length(lines) < 9) return(NA_real_)
  as.numeric(sub("([0-9.]+)%.*", "\\1", lines[9]))
}, numeric(1))
alignment_summary <- data.frame(sample = sub("_align\\.log$", "", basename(alignment_logs)), rate = alignment_rate)
alignment_summary <- alignment_summary[order(alignment_summary$rate), , drop = FALSE]
write.csv(alignment_summary, file.path(output_dir, "alignment_rates.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
