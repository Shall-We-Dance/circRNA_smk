suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

dir.create(dirname(snakemake@log[[1]]), recursive = TRUE, showWarnings = FALSE)
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

message("Starting BSJ DEG analysis.")

counts_path <- snakemake@input[["bsj"]]
groups <- snakemake@params[["groups"]]

counts <- read.table(
  counts_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(counts) < 3) {
  stop("BSJ matrix must contain circRNA id column and at least two sample columns.")
}

id_col <- colnames(counts)[1]
rownames(counts) <- counts[[id_col]]
counts[[id_col]] <- NULL
counts[] <- lapply(counts, function(x) as.integer(round(as.numeric(x))))
counts[is.na(counts)] <- 0L

sample_to_group <- unlist(
  lapply(names(groups), function(group_name) {
    members <- unlist(groups[[group_name]])
    setNames(rep(group_name, length(members)), members)
  })
)

matched_samples <- intersect(names(sample_to_group), colnames(counts))
if (length(matched_samples) < 2) {
  stop("At least two configured samples must be present in BSJ matrix.")
}

design <- data.frame(
  sample = matched_samples,
  group = sample_to_group[matched_samples],
  stringsAsFactors = FALSE
)
design$group <- factor(design$group, levels = sort(unique(design$group)))
rownames(design) <- design$sample

counts <- counts[, design$sample, drop = FALSE]

dir.create(dirname(snakemake@output[["metadata"]]), recursive = TRUE, showWarnings = FALSE)
write.table(
  design,
  file = snakemake@output[["metadata"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts),
  colData = design,
  design = ~ group
)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

norm_means <- rowMeans(counts(dds, normalized = TRUE))
nsub <- min(1000, length(norm_means))
high_count_rows <- sum(norm_means > 5)

if (high_count_rows < nsub) {
  message(
    "Detected only ", high_count_rows,
    " rows with mean normalized count > 5 (< nsub=", nsub, "). ",
    "Using varianceStabilizingTransformation() instead of vst()."
  )
  vst_obj <- varianceStabilizingTransformation(dds, blind = FALSE)
} else {
  vst_obj <- vst(dds, blind = FALSE)
}

vst_mat <- assay(vst_obj)
dir.create(dirname(snakemake@output[["vst_counts"]]), recursive = TRUE, showWarnings = FALSE)
write.table(
  data.frame(circRNA = rownames(vst_mat), vst_mat, check.names = FALSE),
  file = snakemake@output[["vst_counts"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_volcano <- function(res_df, title, out_file, alpha = 0.05) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  res_df$padj[is.na(res_df$padj)] <- 1
  res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0
  res_df$significant <- res_df$padj < alpha
  res_df$minus_log10_padj <- -log10(pmax(res_df$padj, .Machine$double.xmin))

  p <- ggplot(res_df, aes(x = log2FoldChange, y = minus_log10_padj, color = significant)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#d7301f")) +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value",
      color = "padj < 0.05"
    )
  ggsave(out_file, p, width = 7, height = 5)
}

plot_heatmap <- function(vst_values, sig_ids, out_file, title, top_n = 50) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (length(sig_ids) < 2) {
    pdf(out_file, width = 8, height = 6)
    plot.new()
    title(main = sprintf("%s\nNot enough significant BSJs for heatmap.", title))
    dev.off()
    return()
  }
  keep_ids <- sig_ids[seq_len(min(length(sig_ids), top_n))]
  mat <- vst_values[keep_ids, , drop = FALSE]
  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    fontsize_col = 9,
    main = title,
    filename = out_file,
    width = 8,
    height = 10
  )
}

res_all <- results(dds)
res_all_df <- data.frame(circRNA = rownames(res_all), as.data.frame(res_all), check.names = FALSE)
dir.create(dirname(snakemake@output[["all_results"]]), recursive = TRUE, showWarnings = FALSE)
write.table(res_all_df, file = snakemake@output[["all_results"]], sep = "\t", quote = FALSE, row.names = FALSE)

plot_volcano(res_all_df, "BSJ DEG: all groups", snakemake@output[["all_volcano"]])
sig_all <- res_all_df[which(!is.na(res_all_df$padj) & res_all_df$padj < 0.05), ]
sig_all <- sig_all[order(sig_all$padj, -abs(sig_all$log2FoldChange)), ]
plot_heatmap(vst_mat, sig_all$circRNA, snakemake@output[["all_heatmap"]], "Top significant BSJs (all groups)")

pdf(snakemake@output[["pca"]], width = 7, height = 5)
plotPCA(vst_obj, intgroup = "group")
dev.off()

group_levels <- levels(design$group)
pairwise <- combn(group_levels, 2, simplify = FALSE)

pairwise_outputs <- function(output_paths, suffix) {
  matched <- output_paths[grepl(paste0("/", suffix, "$"), output_paths)]
  setNames(matched, basename(dirname(matched)))
}

pair_results <- pairwise_outputs(unlist(snakemake@output[["pairwise_results"]]), "deseq2_results.tsv")
pair_volcano <- pairwise_outputs(unlist(snakemake@output[["pairwise_volcano"]]), "volcano.pdf")
pair_heatmap <- pairwise_outputs(unlist(snakemake@output[["pairwise_heatmap"]]), "heatmap_top50.pdf")

for (pair in pairwise) {
  g1 <- pair[[1]]
  g2 <- pair[[2]]
  pair_name <- paste0(g1, "_vs_", g2)
  message("Running pairwise comparison: ", pair_name)
  res <- results(dds, contrast = c("group", g1, g2))
  res_df <- data.frame(circRNA = rownames(res), as.data.frame(res), check.names = FALSE)
  write.table(res_df, file = pair_results[[pair_name]], sep = "\t", quote = FALSE, row.names = FALSE)
  plot_volcano(res_df, paste0("BSJ DEG: ", g1, " vs ", g2), pair_volcano[[pair_name]])
  sig <- res_df[which(!is.na(res_df$padj) & res_df$padj < 0.05), ]
  sig <- sig[order(sig$padj, -abs(sig$log2FoldChange)), ]
  pair_samples <- design$sample[design$group %in% c(g1, g2)]
  plot_heatmap(
    vst_mat[, pair_samples, drop = FALSE],
    sig$circRNA,
    pair_heatmap[[pair_name]],
    paste0("Top significant BSJs: ", g1, " vs ", g2)
  )
}

message("BSJ DEG analysis finished.")
