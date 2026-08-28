# =============================================================================
# Bulk RNA-seq differential-expression analysis: RGS6 knockout study
# =============================================================================
# Required metadata columns: Mouse_ID, Sample_code, Sex, Mouseline
# Expected factor values: Sex = Female/Male; Mouseline = RGS6KO/WT
# =============================================================================

# 1. Configuration --------------------------------------------------------------
bam_dir <- "path/to/aligned_bam"
results_dir <- "path/to/dge_results"
metadata_file <- file.path(results_dir, "Rorabaugh_FASTQ_Sample_Guide.csv")
gtf_file <- "path/to/Mus_musculus.GRCm38.102.gtf"

# The historical GSEA code breaks equal ranking scores with rnorm(). Leave NULL
# to retain that behavior. Set a value only for deterministic *new* reruns; the
# original run did not record a seed, so its exact GSEA tie resolution cannot be
# reconstructed from code alone.
gsea_seed <- NULL

# 2. Packages -------------------------------------------------------------------
required_packages <- c(
  "Rsubread", "edgeR", "limma", "clusterProfiler", "org.Mm.eg.db",
  "msigdbr", "AnnotationDbi", "dplyr", "tibble", "ggplot2", "ggrepel",
  "enrichplot"
)
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Rsubread); library(edgeR); library(limma); library(clusterProfiler)
  library(org.Mm.eg.db); library(msigdbr); library(AnnotationDbi)
  library(dplyr); library(tibble); library(ggplot2); library(ggrepel); library(enrichplot)
})

if (!dir.exists(bam_dir)) stop("BAM directory does not exist: ", bam_dir)
if (!file.exists(metadata_file)) stop("Metadata file does not exist: ", metadata_file)
if (!file.exists(gtf_file)) stop("GTF file does not exist: ", gtf_file)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Metadata, BAM order, and counting -----------------------------------------
sample_info <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("Mouse_ID", "Sample_code", "Sex", "Mouseline")
if (!all(required_columns %in% names(sample_info))) {
  stop("Metadata is missing: ", paste(setdiff(required_columns, names(sample_info)), collapse = ", "))
}
if (anyDuplicated(sample_info$Sample_code)) stop("Sample_code must be unique.")
if (anyDuplicated(sample_info$Mouse_ID)) stop("Mouse_ID must be unique.")

# Preserved from the original: featureCounts receives BAMs in metadata order.
bam_files <- file.path(bam_dir, paste0(sample_info$Sample_code, ".bam"))
if (any(!file.exists(bam_files))) {
  stop("Missing BAM file(s): ", paste(bam_files[!file.exists(bam_files)], collapse = ", "))
}

fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE
)
save(fc, file = file.path(results_dir, "fc.RData"))

counts <- fc$counts
clean_names <- sub("\\.bam$", "", colnames(counts))
idx <- match(clean_names, sample_info$Sample_code)
stopifnot(!anyNA(idx), !anyDuplicated(sample_info$Sample_code))
stopifnot(identical(clean_names, sample_info$Sample_code[idx]))
colnames(counts) <- sample_info[idx, "Mouse_ID"]

# Do not reorder sample_info here: original BAM construction already used its
# order, and preserving that order is necessary for numerical equivalence.
rownames(sample_info) <- sample_info$Mouse_ID
stopifnot(identical(colnames(counts), rownames(sample_info)))
sample_info$Mouseline <- factor(sample_info$Mouseline, levels = c("RGS6KO", "WT"))
sample_info$Sex <- factor(sample_info$Sex, levels = c("Female", "Male"))
if (anyNA(sample_info$Mouseline) || anyNA(sample_info$Sex)) {
  stop("Sex must be Female/Male and Mouseline must be RGS6KO/WT.")
}
sample_info <- droplevels(sample_info)
write.table(sample_info, file.path(results_dir, "samples_info.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
writeLines(capture.output(table(sample_info$Mouseline, sample_info$Sex)), file.path(results_dir, "samples_info_tab.txt"))

# 4. edgeR/limma-voom model ----------------------------------------------------
design <- model.matrix(~ 0 + Sex * Mouseline, data = sample_info)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

png(file.path(results_dir, "boxplot_before_filtering.png"), width = 2000, height = 1500, res = 300)
boxplot(log2(cpm(dge) + 1), main = "Before Filtering")
dev.off()

keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

png(file.path(results_dir, "voom.png"), width = 2000, height = 1500, res = 300)
v <- voomWithQualityWeights(dge, design, plot = TRUE)
dev.off()
save(v, file = file.path(results_dir, "voom_object.RData"))

png(file.path(results_dir, "boxplot_after_filtering.png"), width = 2000, height = 1500, res = 300)
boxplot(log2(cpm(dge) + 1), main = "After Filtering")
dev.off()

fit <- lmFit(v, design)
# These name normalizations occur in the original before the contrasts.
colnames(fit$coefficients) <- make.names(colnames(fit$coefficients))
colnames(fit$design) <- make.names(colnames(fit$design))
colnames(fit$stdev.unscaled) <- make.names(colnames(fit$stdev.unscaled))
colnames(fit$cov.coefficients) <- make.names(colnames(fit$cov.coefficients))

make_contrasts_final <- function(fit) {
  contrasts <- list(
    WT_vs_RGS6KO_Female = "MouselineWT",
    WT_vs_RGS6KO_Male = "MouselineWT + SexMale.MouselineWT",
    WT_vs_RGS6KO_overall = "MouselineWT + 0.5*SexMale.MouselineWT",
    Male_vs_Female_RGS6KO = "SexMale - SexFemale",
    Male_vs_Female_WT = "SexMale - SexFemale + SexMale.MouselineWT",
    Male_vs_Female_overall = "(SexMale - SexFemale) + 0.5*SexMale.MouselineWT",
    Sex_by_Mouseline_Interaction = "SexMale.MouselineWT",
    Female_RGS6KO = "SexFemale",
    Male_RGS6KO = "SexMale",
    Female_WT = "SexFemale + MouselineWT",
    Male_WT = "SexMale + MouselineWT + SexMale.MouselineWT"
  )
  out <- makeContrasts(contrasts = contrasts, levels = colnames(coef(fit)))
  colnames(out) <- names(contrasts)
  out
}

contrast_matrix <- make_contrasts_final(fit)
write.table(contrast_matrix, file.path(results_dir, "contrast_matrix.tsv"), sep = "\t", quote = FALSE)
stopifnot(qr(design)$rank == ncol(design), is.null(nonEstimable(design)))
stopifnot(identical(rownames(contrast_matrix), colnames(coef(fit))))

fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
top_tables_all <- lapply(colnames(contrast_matrix), function(cn) topTable(fit2, coef = cn, number = Inf, adjust = "BH"))
names(top_tables_all) <- colnames(contrast_matrix)
save(top_tables_all, file = file.path(results_dir, "top_tables_all_strains.RData"))

# Original GSEA uses all non-group-mean contrasts, with the same signed P-value
# ranking and ID conversion rules.
top_tables <- top_tables_all[!grepl("^(Female_RGS6KO|Female_WT|Male_RGS6KO|Male_WT)", names(top_tables_all))]

# 5. GSEA -----------------------------------------------------------------------
prepare_gsea_rank <- function(tbl) {
  tbl <- tbl |>
    mutate(rank_metric = -log10(P.Value) * sign(logFC)) |>
    arrange(desc(rank_metric))
  tbl <- tbl[!duplicated(rownames(tbl)), , drop = FALSE]
  gene_list <- tbl$rank_metric
  names(gene_list) <- rownames(tbl)
  gene_list <- gene_list + rnorm(length(gene_list), sd = 1e-6) # preserved historical tie-breaking
  sort(gene_list, decreasing = TRUE)
}

convert_ids <- function(gene_list) {
  mapping <- bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  mapping <- mapping[!duplicated(mapping$ENSEMBL), , drop = FALSE]
  ranked <- data.frame(ENSEMBL = names(gene_list), rank_metric = unname(gene_list)) |>
    inner_join(mapping, by = "ENSEMBL") |>
    arrange(desc(abs(rank_metric))) |>
    distinct(ENTREZID, .keep_all = TRUE)
  out <- ranked$rank_metric
  names(out) <- ranked$ENTREZID
  sort(out, decreasing = TRUE)
}

run_gsea_go <- function(gene_list) gseGO(
  geneList = gene_list, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "ENTREZID",
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE
)
run_gsea_msigdb <- function(gene_list) {
  term2gene <- msigdbr(species = "Mus musculus", collection = "H") |>
    select(gs_name, ensembl_gene) |>
    rename(term = gs_name, gene = ensembl_gene)
  GSEA(geneList = gene_list, TERM2GENE = term2gene, minGSSize = 10, maxGSSize = 500,
       pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
}

if (!is.null(gsea_seed)) set.seed(gsea_seed)
gsea_results <- lapply(top_tables, function(tbl) {
  gene_list <- prepare_gsea_rank(tbl)
  list(GO = run_gsea_go(convert_ids(gene_list)), MSIG = run_gsea_msigdb(gene_list))
})
save(gsea_results, file = file.path(results_dir, "gsea_results.RData"))

# The original retains the best ten unique descriptions for this displayed plot.
overall_go <- gsea_results$WT_vs_RGS6KO_overall$GO
overall_go@result <- overall_go@result |>
  arrange(p.adjust) |>
  distinct(Description, .keep_all = TRUE) |>
  head(10)
gsea_results$WT_vs_RGS6KO_overall$GO <- overall_go
png(file.path(results_dir, "Figure_GSEA_WT_vs_RGS6KO_overall.png"), width = 10, height = 9, units = "in", res = 300)
gsea_plot <- dotplot(overall_go, showCategory = 10) +
  theme(panel.grid.major = element_line(linewidth = 0.7), panel.grid.minor = element_blank()) +
  scale_size(range = c(5, 11)) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Adjusted P value") +
  theme_classic(base_size = 19) +
  theme(axis.line = element_line(linewidth = 0.8), axis.ticks = element_line(linewidth = 0.8), legend.position = "right")
print(gsea_plot)
dev.off()

top_ids <- overall_go@result$ID[1:min(3, nrow(overall_go@result))]
if (length(top_ids)) {
  png(file.path(results_dir, "Figure_GSEAcurve_WT_vs_RGS6KO_overall.png"), width = 10, height = 7, units = "in", res = 300)
  print(gseaplot2(overall_go, geneSetID = top_ids))
  dev.off()
}

# 6. Published tables -----------------------------------------------------------
# Group means below intentionally use edgeR logCPM, exactly as the original.
expr_logcpm <- cpm(dge, log = TRUE)
gene_symbols <- mapIds(org.Mm.eg.db, keys = rownames(expr_logcpm), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
group <- interaction(sample_info$Mouseline, sample_info$Sex)
group_means <- data.frame(
  ENSEMBL = rownames(expr_logcpm), Symbol = unname(gene_symbols),
  RGS6KO_Female = rowMeans(expr_logcpm[, group == "RGS6KO.Female", drop = FALSE]),
  RGS6KO_Male = rowMeans(expr_logcpm[, group == "RGS6KO.Male", drop = FALSE]),
  WT_Female = rowMeans(expr_logcpm[, group == "WT.Female", drop = FALSE]),
  WT_Male = rowMeans(expr_logcpm[, group == "WT.Male", drop = FALSE])
)
write.csv(filter(group_means, grepl("^Rgs", Symbol)), file.path(results_dir, "RGS_family_expression.csv"), row.names = FALSE)
write.csv(group_means, file.path(results_dir, "All_detected_transcripts.csv"), row.names = FALSE)

deg <- topTable(fit2, coef = "WT_vs_RGS6KO_overall", number = Inf, sort.by = "P")
deg$ENSEMBL <- rownames(deg)
anno <- AnnotationDbi::select(org.Mm.eg.db, keys = deg$ENSEMBL, columns = "SYMBOL", keytype = "ENSEMBL") |>
  filter(!is.na(SYMBOL)) |>
  distinct(ENSEMBL, .keep_all = TRUE)
deg <- deg |> left_join(anno, by = "ENSEMBL") |> filter(!is.na(SYMBOL))
table1 <- deg |>
  left_join(group_means, by = c("SYMBOL" = "Symbol")) |>
  filter(adj.P.Val < 0.05) |>
  arrange(adj.P.Val) |>
  transmute(
    Gene = SYMBOL, `KO Female` = sprintf("%.2f", RGS6KO.Female), `WT Female` = sprintf("%.2f", WT.Female),
    `KO Male` = sprintf("%.2f", RGS6KO.Male), `WT Male` = sprintf("%.2f", WT.Male),
    `log₂FC (WT vs KO)` = sprintf("%.2f", logFC), Direction = ifelse(logFC > 0, "Higher in WT", "Higher in KO"),
    B = round(B, 2), FDR = format(adj.P.Val, digits = 3, scientific = TRUE)
  ) |>
  mutate(order = ifelse(Gene == "Rgs6", 0, 1)) |>
  arrange(order, FDR) |>
  select(-order)
write.csv(table1, file.path(results_dir, "Table1_DifferentialExpression.csv"), row.names = FALSE, quote = FALSE)

individual_genes <- c("Rgs6", "Snapc1", "Gm4742", "Churc1", "Ap4s1", "Gpr135", "Cish")
individual_indices <- which(gene_symbols %in% individual_genes)
individual_expression <- data.frame(Mouse = as.integer(colnames(v$E)), t(v$E[individual_indices, , drop = FALSE]), check.names = FALSE)
names(individual_expression)[-1] <- unname(gene_symbols[individual_indices])
individual_expression <- left_join(sample_info, individual_expression, by = c("Mouse_ID" = "Mouse"))
write.csv(individual_expression, file.path(results_dir, "Table1_individual_mouse_expression.csv"), row.names = FALSE)

# 7. Figure 6 ------------------------------------------------------------------
rgs_order <- c("Rgs9", "Rgs4", "Rgs8", "Rgs17", "Rgs2", "Rgs5", "Rgs14", "Rgs20", "Rgs7", "Rgs10", "Rgs11", "Rgs3", "Rgs12", "Rgs19", "Rgs6", "Rgs16")
expr_rgs <- v$E[gene_symbols %in% rgs_order, , drop = FALSE]
rownames(expr_rgs) <- gene_symbols[gene_symbols %in% rgs_order]
expr_rgs <- expr_rgs[rgs_order, , drop = FALSE]
fig6b_data <- data.frame(Gene = names(rowMeans(expr_rgs)), logCPM = rowMeans(expr_rgs)) |>
  arrange(logCPM) |>
  mutate(Gene = factor(Gene, levels = Gene))
fig6b <- ggplot(fig6b_data, aes(logCPM, Gene)) +
  geom_segment(aes(x = 0, xend = logCPM, yend = Gene), linewidth = 0.7, color = "grey70") +
  geom_point(size = 3.8, color = "#3B4CC0") +
  labs(x = "Average normalized expression (logCPM)", y = NULL) + theme_classic(base_size = 13)
ggsave(file.path(results_dir, "Figure6_RGS_lollipop.png"), fig6b, width = 5, height = 5, dpi = 300)

tab <- topTable(fit2, coef = "WT_vs_RGS6KO_overall", number = Inf, sort.by = "none")
tab$ENSEMBL <- rownames(tab)
tab <- left_join(tab, anno, by = "ENSEMBL")
tab <- tab |> mutate(
  negLogFDR = -log10(adj.P.Val),
  Group = case_when(SYMBOL == "Rgs6" ~ "Rgs6", SYMBOL %in% rgs_order ~ "RGS family", TRUE ~ "Other")
)
labels <- filter(tab, SYMBOL == "Rgs6" | (adj.P.Val < 0.05 & abs(logFC) > 1))
fig6a <- ggplot(tab, aes(logFC, negLogFDR)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey50") +
  geom_point(aes(colour = Group), alpha = .75, size = 2.2) +
  geom_text_repel(data = labels, aes(label = SYMBOL), size = 4, max.overlaps = Inf, box.padding = .4, point.padding = .2, seed = 123) +
  scale_colour_manual(values = c(Other = "grey75", `RGS family` = "#2C7FB8", Rgs6 = "#C00000")) +
  labs(x = expression(log[2] * " fold-change (WT relative to RGS6KO)"), y = expression(-log[10] * "(FDR)")) +
  theme_classic(base_size = 16) + theme(legend.position = "top", legend.title = element_blank())
ggsave(file.path(results_dir, "Figure6A_volcano.pdf"), fig6a, width = 7, height = 6)

extract_limma_effect <- function(coef_name) {
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "none")
  tt$ENSEMBL <- rownames(tt)
  tt <- left_join(tt, anno, by = "ENSEMBL") |> filter(SYMBOL %in% rgs_order)
  i <- match(tt$ENSEMBL, rownames(fit2))
  se <- fit2$sigma[i] * fit2$stdev.unscaled[i, coef_name]
  mutate(tt, Lower = logFC - se, Upper = logFC + se)
}
forest_data <- bind_rows(
  mutate(extract_limma_effect("WT_vs_RGS6KO_Female"), Sex = "Female"),
  mutate(extract_limma_effect("WT_vs_RGS6KO_Male"), Sex = "Male")
) |> transmute(Gene = SYMBOL, Sex, Effect = logFC, Lower, Upper)
forest_data$Gene <- factor(forest_data$Gene, levels = rev(rgs_order))
forest_data$Highlight <- ifelse(forest_data$Gene == "Rgs6", "Rgs6", "Other")
xmax <- 1.10 * max(abs(forest_data$Lower), abs(forest_data$Upper))
fig6c <- ggplot(forest_data, aes(Effect, Gene)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, colour = Highlight), height = .12, linewidth = .6) +
  geom_point(aes(colour = Highlight, shape = Sex), size = 3) +
  scale_colour_manual(values = c(Other = "grey55", Rgs6 = "#B30000")) +
  scale_shape_manual(values = c(Female = 16, Male = 17)) +
  scale_x_continuous(limits = c(-xmax, xmax), expand = expansion(mult = .02)) +
  labs(x = "Estimated genotype effect (WT − RGS6KO, log2FC)", y = NULL) +
  theme_classic(base_size = 15) + theme(legend.position = "top", legend.title = element_blank(), axis.text.y = element_text(face = "italic", size = 14))
ggsave(file.path(results_dir, "Figure6C_limma_forest.png"), fig6c, width = 7, height = 7, dpi = 300)

# 8. GEO files ------------------------------------------------------------------
write.csv(v$E, file.path(results_dir, "RNAseq_voom_normalized_expression_logCPM.csv"), quote = FALSE)
write.csv(counts, file.path(results_dir, "RNAseq_raw_counts.csv"), quote = FALSE)
write.csv(sample_info, file.path(results_dir, "RNAseq_sample_metadata.csv"), row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
