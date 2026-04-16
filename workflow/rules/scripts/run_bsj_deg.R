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
ciri3_path <- snakemake@input[["ciri3"]]
groups <- snakemake@params[["groups"]]
padj_cutoff <- as.numeric(snakemake@params[["padj_cutoff"]])
lfc_cutoff <- as.numeric(snakemake@params[["lfc_cutoff"]])
min_total_count <- as.numeric(snakemake@params[["min_total_count"]])
min_samples_detected <- as.numeric(snakemake@params[["min_samples_detected"]])

if (!is.finite(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff >= 1) {
  stop("deg.padj_cutoff must be a numeric value between 0 and 1.")
}
if (!is.finite(lfc_cutoff) || lfc_cutoff < 0) {
  stop("deg.lfc_cutoff must be a non-negative numeric value.")
}
if (!is.finite(min_total_count) || min_total_count < 0) {
  stop("deg.min_total_count must be a non-negative numeric value.")
}
if (!is.finite(min_samples_detected) || min_samples_detected < 1) {
  stop("deg.min_samples_detected must be at least 1.")
}
min_total_count <- as.integer(round(min_total_count))
min_samples_detected <- as.integer(round(min_samples_detected))

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
if (id_col != "circRNA_ID") {
  warning(
    "The first column is '", id_col,
    "' (expected 'circRNA_ID'). It will still be treated as BSJ identifier."
  )
}
rownames(counts) <- counts[[id_col]]
counts[[id_col]] <- NULL
counts[] <- lapply(counts, function(x) as.integer(round(as.numeric(x))))
counts[is.na(counts)] <- 0L

parse_circrna_id <- function(ids) {
  seqname <- sub(":.*$", "", ids)
  coord_part <- sub("^[^:]+:", "", ids)
  start <- suppressWarnings(as.integer(sub("\\|.*$", "", coord_part)))
  end <- suppressWarnings(as.integer(sub("^.*\\|", "", coord_part)))
  data.frame(
    circRNA = ids,
    seqname = seqname,
    start = start,
    end = end,
    stringsAsFactors = FALSE
  )
}
bsj_annot <- parse_circrna_id(rownames(counts))

read_ciri3_annotation <- function(path) {
  annot <- tryCatch(
    read.table(
      path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = ""
    ),
    error = function(e) {
      warning("Unable to parse CIRI3 table for annotation: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(annot) || nrow(annot) == 0 || !("circRNA_ID" %in% colnames(annot))) {
    return(data.frame(
      circRNA = rownames(counts),
      gene_id = NA_character_,
      strand = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  map_column <- function(df, candidates) {
    for (nm in candidates) {
      if (nm %in% colnames(df)) return(nm)
    }
    NULL
  }

  gene_col <- map_column(annot, c("gene_id", "gene", "gene_name"))
  strand_col <- map_column(annot, c("strand", "Strand", "bsj_strand"))

  one_value <- function(v) {
    vals <- unique(v[!is.na(v) & nzchar(trimws(v))])
    if (length(vals) == 0) return(NA_character_)
    paste(vals, collapse = ";")
  }

  split_idx <- split(seq_len(nrow(annot)), annot[["circRNA_ID"]])
  out <- data.frame(
    circRNA = names(split_idx),
    gene_id = vapply(split_idx, function(idx) {
      if (is.null(gene_col)) return(NA_character_)
      one_value(as.character(annot[idx, gene_col]))
    }, character(1)),
    strand = vapply(split_idx, function(idx) {
      if (is.null(strand_col)) return(NA_character_)
      one_value(as.character(annot[idx, strand_col]))
    }, character(1)),
    stringsAsFactors = FALSE
  )
  out
}
ciri3_annot <- read_ciri3_annotation(ciri3_path)

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

keep <- rowSums(counts) >= min_total_count &
  rowSums(counts > 0) >= min_samples_detected
counts <- counts[keep, , drop = FALSE]
bsj_annot <- bsj_annot[match(rownames(counts), bsj_annot$circRNA), , drop = FALSE]

if (nrow(counts) < 2) {
  stop(
    "Not enough BSJs left after filtering (n=", nrow(counts), "). ",
    "Adjust deg.min_total_count/min_samples_detected."
  )
}
message("BSJs after filtering: ", nrow(counts))

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

plot_volcano <- function(res_df, title, out_file, alpha = 0.05, lfc_threshold = 1) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  res_df$padj[is.na(res_df$padj)] <- 1
  res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0
  res_df$reg_direction <- "Not_significant"
  res_df$reg_direction[res_df$padj < alpha & res_df$log2FoldChange >= lfc_threshold] <- "Up"
  res_df$reg_direction[res_df$padj < alpha & res_df$log2FoldChange <= -lfc_threshold] <- "Down"
  res_df$minus_log10_padj <- -log10(pmax(res_df$padj, .Machine$double.xmin))

  p <- ggplot(res_df, aes(x = log2FoldChange, y = minus_log10_padj, color = reg_direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Up" = "#d7301f", "Down" = "#2c7fb8", "Not_significant" = "grey70")) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value",
      color = paste0("padj < ", alpha, "\n& |log2FC| >= ", lfc_threshold)
    )
  ggsave(out_file, p, width = 7, height = 5)
}

plot_volcano_labeled <- function(res_df, title, out_file, alpha = 0.05, lfc_threshold = 1, label_n = 15) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  res_df$padj[is.na(res_df$padj)] <- 1
  res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0
  res_df$reg_direction <- "Not_significant"
  res_df$reg_direction[res_df$padj < alpha & res_df$log2FoldChange >= lfc_threshold] <- "Up"
  res_df$reg_direction[res_df$padj < alpha & res_df$log2FoldChange <= -lfc_threshold] <- "Down"
  res_df$minus_log10_padj <- -log10(pmax(res_df$padj, .Machine$double.xmin))
  res_df$gene_label <- ifelse(
    is.na(res_df$gene_id) | res_df$gene_id == "",
    res_df$circRNA,
    res_df$gene_id
  )
  sig_df <- res_df[res_df$reg_direction != "Not_significant", , drop = FALSE]
  sig_df <- sig_df[order(sig_df$padj, -abs(sig_df$log2FoldChange)), , drop = FALSE]
  if (nrow(sig_df) > label_n) {
    sig_df <- sig_df[seq_len(label_n), , drop = FALSE]
  }

  p <- ggplot(res_df, aes(x = log2FoldChange, y = minus_log10_padj, color = reg_direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_text(
      data = sig_df,
      aes(label = gene_label),
      size = 2.8,
      check_overlap = TRUE,
      vjust = -0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("Up" = "#d7301f", "Down" = "#2c7fb8", "Not_significant" = "grey70")) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value",
      color = paste0("padj < ", alpha, "\n& |log2FC| >= ", lfc_threshold)
    )
  ggsave(out_file, p, width = 7, height = 5)
}

plot_heatmap <- function(vst_values, sig_ids, feature_labels, out_file, title, col_groups, top_n = 50) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (length(sig_ids) < 2) {
    pdf(out_file, width = 8, height = 6)
    plot.new()
    title(main = sprintf("%s\nNot enough significant BSJs for heatmap.", title))
    dev.off()
    return()
  }
  keep_ids <- sig_ids[seq_len(min(length(sig_ids), top_n))]
  sample_ids <- intersect(names(col_groups), colnames(vst_values))
  if (length(sample_ids) == 0) {
    stop("No overlapping samples between heatmap matrix and group annotations.")
  }
  sample_order <- sample_ids[order(col_groups[sample_ids], sample_ids)]
  mat <- vst_values[keep_ids, sample_order, drop = FALSE]
  row_labels <- feature_labels[keep_ids]
  row_labels[is.na(row_labels) | row_labels == ""] <- keep_ids[is.na(row_labels) | row_labels == ""]
  row_labels <- make.unique(row_labels, sep = "_")
  annotation_col <- data.frame(
    group = factor(col_groups[sample_order], levels = unique(col_groups[sample_order])),
    row.names = sample_order,
    stringsAsFactors = FALSE
  )
  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    labels_row = row_labels,
    show_rownames = TRUE,
    fontsize_col = 9,
    main = title,
    color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(100),
    filename = out_file,
    width = 8,
    height = 10
  )
}

res_all <- results(dds)
res_all_df <- data.frame(circRNA = rownames(res_all), as.data.frame(res_all), check.names = FALSE)
res_all_df <- merge(bsj_annot, res_all_df, by = "circRNA", all.y = TRUE, sort = FALSE)
res_all_df <- merge(res_all_df, ciri3_annot, by = "circRNA", all.x = TRUE, sort = FALSE)
res_all_df$regulation <- ifelse(
  !is.na(res_all_df$padj) & res_all_df$padj < padj_cutoff & res_all_df$log2FoldChange >= lfc_cutoff,
  "Up",
  ifelse(
    !is.na(res_all_df$padj) & res_all_df$padj < padj_cutoff & res_all_df$log2FoldChange <= -lfc_cutoff,
    "Down",
    "Not_significant"
  )
)
dir.create(dirname(snakemake@output[["all_results"]]), recursive = TRUE, showWarnings = FALSE)
write.table(res_all_df, file = snakemake@output[["all_results"]], sep = "\t", quote = FALSE, row.names = FALSE)

feature_label_all <- setNames(
  ifelse(
    is.na(res_all_df$gene_id) | res_all_df$gene_id == "",
    res_all_df$circRNA,
    paste0(res_all_df$gene_id, " (", res_all_df$circRNA, ")")
  ),
  res_all_df$circRNA
)

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
pair_volcano_labeled <- pairwise_outputs(unlist(snakemake@output[["pairwise_volcano_labeled"]]), "volcano_labeled.pdf")
pair_heatmap <- pairwise_outputs(unlist(snakemake@output[["pairwise_heatmap"]]), "heatmap_top50.pdf")
pair_pca <- pairwise_outputs(unlist(snakemake@output[["pairwise_pca"]]), "pca.pdf")
pairwise_sig_rank <- data.frame(circRNA = character(0), padj = numeric(0), stringsAsFactors = FALSE)

for (pair in pairwise) {
  g1 <- pair[[1]]
  g2 <- pair[[2]]
  pair_name <- paste0(g1, "_vs_", g2)
  message("Running pairwise comparison: ", pair_name)
  res <- results(dds, contrast = c("group", g2, g1))
  res_df <- data.frame(circRNA = rownames(res), as.data.frame(res), check.names = FALSE)
  res_df <- merge(bsj_annot, res_df, by = "circRNA", all.y = TRUE, sort = FALSE)
  res_df <- merge(res_df, ciri3_annot, by = "circRNA", all.x = TRUE, sort = FALSE)
  res_df$regulation <- ifelse(
    !is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >= lfc_cutoff,
    paste0("Up_in_", g2),
    ifelse(
      !is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange <= -lfc_cutoff,
      paste0("Up_in_", g1),
      "Not_significant"
    )
  )
  write.table(res_df, file = pair_results[[pair_name]], sep = "\t", quote = FALSE, row.names = FALSE)
  plot_volcano(
    res_df,
    paste0("BSJ DEG: ", g2, " vs ", g1),
    pair_volcano[[pair_name]],
    alpha = padj_cutoff,
    lfc_threshold = lfc_cutoff
  )
  plot_volcano_labeled(
    res_df,
    paste0("BSJ DEG (labeled): ", g2, " vs ", g1),
    pair_volcano_labeled[[pair_name]],
    alpha = padj_cutoff,
    lfc_threshold = lfc_cutoff
  )
  sig <- res_df[which(!is.na(res_df$padj) & res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) >= lfc_cutoff), ]
  sig <- sig[order(sig$padj, -abs(sig$log2FoldChange)), ]
  if (nrow(sig) > 0) {
    pairwise_sig_rank <- rbind(pairwise_sig_rank, sig[, c("circRNA", "padj"), drop = FALSE])
  }
  pair_samples <- design$sample[design$group %in% c(g1, g2)]
  pair_vst <- vst_mat[, pair_samples, drop = FALSE]
  pair_design <- design[pair_samples, , drop = FALSE]
  feature_labels <- setNames(
    ifelse(
      is.na(res_df$gene_id) | res_df$gene_id == "",
      res_df$circRNA,
      paste0(res_df$gene_id, " (", res_df$circRNA, ")")
    ),
    res_df$circRNA
  )
  plot_heatmap(
    pair_vst,
    sig$circRNA,
    feature_labels,
    pair_heatmap[[pair_name]],
    paste0("Top significant BSJs: ", g1, " vs ", g2),
    setNames(as.character(pair_design$group), pair_design$sample)
  )

  pca_obj <- prcomp(t(pair_vst), center = TRUE, scale. = FALSE)
  pca_samples <- rownames(pca_obj$x)
  pca_groups <- as.character(pair_design[pca_samples, "group"])
  pca_df <- data.frame(
    sample = pca_samples,
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    group = pca_groups,
    stringsAsFactors = FALSE
  )
  ve <- pca_obj$sdev^2 / sum(pca_obj$sdev^2)
  pp <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(vjust = -0.5, size = 3, show.legend = FALSE) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0("Pairwise PCA: ", g1, " vs ", g2),
      x = sprintf("PC1 (%.1f%%)", ve[1] * 100),
      y = sprintf("PC2 (%.1f%%)", ve[2] * 100),
      color = "Group"
    )
  ggsave(pair_pca[[pair_name]], pp, width = 7, height = 5)
}

if (nrow(pairwise_sig_rank) > 0) {
  min_padj <- aggregate(padj ~ circRNA, data = pairwise_sig_rank, FUN = min)
  min_padj <- min_padj[order(min_padj$padj), , drop = FALSE]
  plot_heatmap(
    vst_mat,
    min_padj$circRNA,
    feature_label_all,
    snakemake@output[["all_heatmap"]],
    "Top significant BSJs (union of pairwise DEGs)",
    setNames(as.character(design$group), design$sample)
  )
} else {
  plot_heatmap(
    vst_mat,
    character(0),
    feature_label_all,
    snakemake@output[["all_heatmap"]],
    "Top significant BSJs (union of pairwise DEGs)",
    setNames(as.character(design$group), design$sample)
  )
}

message("BSJ DEG analysis finished.")
