suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
})

get_named <- function(x, name, default = NULL) {
  out <- tryCatch(x[[name]], error = function(e) default)
  if (is.null(out) || length(out) == 0) default else out
}

input_path <- function(name, default = NA_character_) {
  value <- get_named(snakemake@input, name, default)
  if (is.null(value) || length(value) == 0) return(default)
  as.character(value[[1]])
}

input_paths <- function(name) {
  value <- get_named(snakemake@input, name, character(0))
  as.character(unlist(value))
}

output_path <- function(name, default = NA_character_) {
  value <- get_named(snakemake@output, name, default)
  if (is.null(value) || length(value) == 0) return(default)
  as.character(value[[1]])
}

log_path <- if (length(snakemake@log) > 0) as.character(snakemake@log[[1]]) else NA_character_
if (!is.na(log_path) && nzchar(log_path)) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(log_path, open = "wt")
  sink(log_con)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)
}

method <- as.character(snakemake@params[["method"]])
mode <- as.character(snakemake@params[["mode"]])
groups <- snakemake@params[["groups"]]
comparisons <- snakemake@params[["comparisons"]]
padj_cutoff <- as.numeric(snakemake@params[["padj_cutoff"]])
lfc_cutoff <- as.numeric(snakemake@params[["lfc_cutoff"]])
comparison_name <- as.character(get_named(snakemake@params, "comparison", ""))

message("Starting CIRI3 DEG plotting: method=", method, ", mode=", mode)

if (!method %in% c("de_bsj", "de_ratio", "de_relative")) {
  stop("Unsupported CIRI3 DEG method: ", method)
}
if (!mode %in% c("pairwise", "all_groups")) {
  stop("Unsupported CIRI3 DEG plotting mode: ", mode)
}
if (!is.finite(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff >= 1) {
  stop("deg.padj_cutoff must be a numeric value between 0 and 1.")
}
if (!is.finite(lfc_cutoff) || lfc_cutoff < 0) {
  stop("deg.lfc_cutoff must be a non-negative numeric value.")
}

method_label <- switch(
  method,
  de_bsj = "CIRI3 DE_BSJ",
  de_ratio = "CIRI3 DE_Ratio",
  de_relative = "CIRI3 DE_Relative"
)

first_existing_column <- function(df, candidates) {
  for (candidate in candidates) {
    if (candidate %in% colnames(df)) return(candidate)
  }
  NULL
}

read_tabular_loose <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("Empty CIRI3 DEG result: ", path)
  }

  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  header <- split_lines[[1]]
  rows <- split_lines[-1]
  if (length(rows) == 0) {
    df <- as.data.frame(matrix(nrow = 0, ncol = length(header)), stringsAsFactors = FALSE)
    colnames(df) <- make.unique(header, sep = "_")
    return(df)
  }

  max_cols <- max(length(header), vapply(rows, length, integer(1)))
  if (length(header) == max_cols - 1) {
    header <- c("row_id", header)
  } else if (length(header) < max_cols) {
    header <- c(header, paste0("extra_", seq_len(max_cols - length(header))))
  } else if (length(header) > max_cols) {
    header <- header[seq_len(max_cols)]
  }
  header <- make.unique(header, sep = "_")

  mat <- do.call(rbind, lapply(rows, function(row) {
    length(row) <- max_cols
    row[is.na(row)] <- ""
    row
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- header
  df
}

to_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

looks_like_circ <- function(x) {
  x <- as.character(x)
  grepl("^[^:]+:[0-9]+\\|[0-9]+", x)
}

read_ciri3_result <- function(path) {
  df <- read_tabular_loose(path)
  if (nrow(df) == 0) {
    df$circRNA_raw <- character(0)
    df$circRNA <- character(0)
    df$gene_id <- character(0)
    df$pvalue <- numeric(0)
    df$padj <- numeric(0)
    df$effect <- numeric(0)
    df$effect_from_result <- logical(0)
    return(df)
  }

  id_col <- first_existing_column(df, c("circRNA_ID", "circRNA", "CircRNA", "ID"))
  raw_id <- NULL
  if (!is.null(id_col) && any(looks_like_circ(df[[id_col]]))) {
    raw_id <- as.character(df[[id_col]])
  } else if ("row_id" %in% colnames(df) && any(looks_like_circ(df[["row_id"]]))) {
    raw_id <- as.character(df[["row_id"]])
  } else {
    for (candidate in colnames(df)) {
      if (any(looks_like_circ(df[[candidate]]))) {
        raw_id <- as.character(df[[candidate]])
        break
      }
    }
  }
  if (is.null(raw_id)) {
    raw_id <- as.character(df[[1]])
    warning("Could not identify a circRNA ID column in ", path, "; using the first column.")
  }

  df$circRNA_raw <- raw_id
  df$circRNA <- sub(";.*$", "", raw_id)
  gene_from_id <- ifelse(grepl(";", raw_id), sub("^[^;]*;", "", raw_id), NA_character_)
  gene_from_id[is.na(gene_from_id) | !nzchar(gene_from_id)] <- NA_character_
  df$gene_id <- gene_from_id

  p_col <- first_existing_column(
    df,
    c("PValue", "P.Value", "pvalue", "p_value", "Pvalue", "PVAL", "P.Val", "pval")
  )
  fdr_col <- first_existing_column(
    df,
    c("FDR", "padj", "adj.P.Val", "qvalue", "q_value", "QValue", "BH")
  )
  effect_col <- first_existing_column(
    df,
    c("logFC", "log2FoldChange", "log2FC", "log2_fold_change", "Log2FC")
  )

  df$pvalue <- if (!is.null(p_col)) to_numeric(df[[p_col]]) else NA_real_
  if (!is.null(fdr_col)) {
    df$padj <- to_numeric(df[[fdr_col]])
  } else if (any(!is.na(df$pvalue))) {
    df$padj <- p.adjust(df$pvalue, method = "BH")
  } else {
    df$padj <- NA_real_
  }
  df$effect <- if (!is.null(effect_col)) to_numeric(df[[effect_col]]) else NA_real_
  df$effect_from_result <- !is.null(effect_col) && any(!is.na(df$effect))
  df
}

one_value <- function(v) {
  vals <- unique(v[!is.na(v) & nzchar(trimws(v))])
  if (length(vals) == 0) return(NA_character_)
  paste(vals, collapse = ";")
}

read_ciri3_annotation <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(data.frame(circRNA = character(0), gene_id_annot = character(0), stringsAsFactors = FALSE))
  }
  annot <- tryCatch(
    read.table(
      path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      fill = TRUE
    ),
    error = function(e) {
      warning("Unable to parse CIRI3 annotation table: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(annot) || nrow(annot) == 0 || !("circRNA_ID" %in% colnames(annot))) {
    return(data.frame(circRNA = character(0), gene_id_annot = character(0), stringsAsFactors = FALSE))
  }

  gene_col <- first_existing_column(annot, c("gene_id", "gene", "gene_name"))
  if (is.null(gene_col)) {
    return(data.frame(circRNA = character(0), gene_id_annot = character(0), stringsAsFactors = FALSE))
  }

  split_idx <- split(seq_len(nrow(annot)), annot[["circRNA_ID"]])
  data.frame(
    circRNA = names(split_idx),
    gene_id_annot = vapply(split_idx, function(idx) one_value(as.character(annot[idx, gene_col])), character(1)),
    stringsAsFactors = FALSE
  )
}

add_annotation <- function(df, annot) {
  if (nrow(df) == 0 || nrow(annot) == 0) return(df)
  merged <- merge(df, annot, by = "circRNA", all.x = TRUE, sort = FALSE)
  missing_gene <- is.na(merged$gene_id) | !nzchar(merged$gene_id)
  merged$gene_id[missing_gene] <- merged$gene_id_annot[missing_gene]
  merged$gene_id_annot <- NULL
  merged
}

read_matrix <- function(path) {
  if (is.na(path) || !nzchar(path)) stop("Missing matrix input path.")
  df <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
  if (ncol(df) < 2) {
    stop("Matrix must contain a circRNA ID column and at least one sample column: ", path)
  }
  ids <- as.character(df[[1]])
  mat_df <- df[, -1, drop = FALSE]
  mat <- as.matrix(data.frame(lapply(mat_df, to_numeric), check.names = FALSE))
  colnames(mat) <- colnames(mat_df)
  rownames(mat) <- ids
  mat[is.na(mat)] <- 0
  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }
  mat
}

build_value_matrix <- function(method, bsj_path, fsj_path = NA_character_) {
  bsj <- read_matrix(bsj_path)
  if (method == "de_ratio") {
    fsj <- read_matrix(fsj_path)
    common_rows <- intersect(rownames(bsj), rownames(fsj))
    common_cols <- intersect(colnames(bsj), colnames(fsj))
    if (length(common_rows) == 0 || length(common_cols) == 0) {
      stop("BSJ and FSJ matrices have no overlapping rows or samples.")
    }
    bsj <- bsj[common_rows, common_cols, drop = FALSE]
    fsj <- fsj[common_rows, common_cols, drop = FALSE]
    denom <- 2 * bsj + fsj
    ratio <- (2 * bsj) / denom
    ratio[!is.finite(ratio)] <- 0
    return(ratio)
  }
  log2(bsj + 1)
}

sample_group_map <- function(groups) {
  out <- unlist(lapply(names(groups), function(group_name) {
    members <- as.character(unlist(groups[[group_name]]))
    setNames(rep(group_name, length(members)), members)
  }))
  out[!duplicated(names(out))]
}

plot_blank <- function(out_file, title, message, width = 7, height = 5) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  pdf(out_file, width = width, height = height)
  plot.new()
  title(main = paste(title, message, sep = "\n"))
  dev.off()
}

label_for_rows <- function(result_df, ids) {
  label_df <- result_df[!duplicated(result_df$circRNA), , drop = FALSE]
  labels <- setNames(label_df$circRNA, label_df$circRNA)
  has_gene <- !is.na(label_df$gene_id) & nzchar(label_df$gene_id)
  labels[label_df$circRNA[has_gene]] <- paste0(label_df$gene_id[has_gene], " (", label_df$circRNA[has_gene], ")")
  out <- labels[ids]
  out[is.na(out) | !nzchar(out)] <- ids[is.na(out) | !nzchar(out)]
  make.unique(out, sep = "_")
}

rank_significant <- function(result_df, require_effect = FALSE, effect_threshold = 0) {
  if (nrow(result_df) == 0) return(result_df)
  sig <- result_df[!is.na(result_df$padj) & result_df$padj < padj_cutoff, , drop = FALSE]
  if (require_effect) {
    sig <- sig[!is.na(sig$effect) & abs(sig$effect) >= effect_threshold, , drop = FALSE]
  }
  if (nrow(sig) == 0) return(sig)
  sig[order(sig$padj, -abs(sig$effect)), , drop = FALSE]
}

plot_heatmap <- function(value_mat, result_df, feature_ids, out_file, title, group_map, top_n = 50) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  feature_ids <- unique(feature_ids)
  feature_ids <- feature_ids[feature_ids %in% rownames(value_mat)]
  if (length(feature_ids) < 2) {
    plot_blank(out_file, title, "Not enough significant circRNAs for heatmap.", width = 8, height = 6)
    return()
  }

  sample_ids <- intersect(names(group_map), colnames(value_mat))
  if (length(sample_ids) < 2) {
    plot_blank(out_file, title, "Not enough samples for heatmap.", width = 8, height = 6)
    return()
  }
  sample_order <- sample_ids[order(group_map[sample_ids], sample_ids)]
  keep_ids <- feature_ids[seq_len(min(length(feature_ids), top_n))]
  mat <- value_mat[keep_ids, sample_order, drop = FALSE]
  mat[!is.finite(mat)] <- 0

  annotation_col <- data.frame(
    group = factor(group_map[sample_order], levels = unique(group_map[sample_order])),
    row.names = sample_order,
    stringsAsFactors = FALSE
  )
  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    labels_row = label_for_rows(result_df, keep_ids),
    show_rownames = TRUE,
    fontsize_col = 9,
    main = title,
    color = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(100),
    filename = out_file,
    width = 8,
    height = 10
  )
}

plot_pca <- function(value_mat, out_file, title, group_map) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  sample_ids <- intersect(names(group_map), colnames(value_mat))
  if (length(sample_ids) < 2) {
    plot_blank(out_file, title, "Not enough samples for PCA.")
    return()
  }

  mat <- value_mat[, sample_ids, drop = FALSE]
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x)) && stats::var(x) > 0), , drop = FALSE]
  if (nrow(mat) < 2) {
    plot_blank(out_file, title, "Not enough variable circRNAs for PCA.")
    return()
  }

  pca_obj <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pca_samples <- rownames(pca_obj$x)
  pca_df <- data.frame(
    sample = pca_samples,
    PC1 = pca_obj$x[, 1],
    PC2 = if (ncol(pca_obj$x) >= 2) pca_obj$x[, 2] else 0,
    group = group_map[pca_samples],
    stringsAsFactors = FALSE
  )
  ve <- pca_obj$sdev^2
  ve <- if (sum(ve) > 0) ve / sum(ve) else rep(NA_real_, length(ve))
  pc1_label <- if (length(ve) >= 1 && is.finite(ve[1])) sprintf("PC1 (%.1f%%)", ve[1] * 100) else "PC1"
  pc2_label <- if (length(ve) >= 2 && is.finite(ve[2])) sprintf("PC2 (%.1f%%)", ve[2] * 100) else "PC2"

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(vjust = -0.5, size = 3, show.legend = FALSE) +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      x = pc1_label,
      y = pc2_label,
      color = "Group"
    )
  ggsave(out_file, p, width = 7, height = 5)
}

matrix_effect <- function(result_df, value_mat, group_map, case_group, control_group) {
  case_samples <- names(group_map)[group_map == case_group]
  control_samples <- names(group_map)[group_map == control_group]
  case_samples <- intersect(case_samples, colnames(value_mat))
  control_samples <- intersect(control_samples, colnames(value_mat))
  if (length(case_samples) == 0 || length(control_samples) == 0) {
    return(rep(NA_real_, nrow(result_df)))
  }
  out <- rep(NA_real_, nrow(result_df))
  matched <- result_df$circRNA %in% rownames(value_mat)
  if (any(matched)) {
    mat <- value_mat[result_df$circRNA[matched], , drop = FALSE]
    out[matched] <- rowMeans(mat[, case_samples, drop = FALSE]) -
      rowMeans(mat[, control_samples, drop = FALSE])
  }
  out
}

plot_volcano <- function(result_df, out_file, title, x_label, effect_threshold, use_effect_threshold = TRUE) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (nrow(result_df) == 0 || all(is.na(result_df$padj))) {
    plot_blank(out_file, title, "No plottable CIRI3 DEG results.")
    return()
  }

  plot_df <- result_df
  plot_df$padj[is.na(plot_df$padj)] <- 1
  plot_df$effect[is.na(plot_df$effect)] <- 0
  plot_df$minus_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
  plot_df$reg_direction <- "Not_significant"
  sig <- plot_df$padj < padj_cutoff
  if (use_effect_threshold) {
    sig <- sig & abs(plot_df$effect) >= effect_threshold
  }
  plot_df$reg_direction[sig & plot_df$effect >= 0] <- "Up"
  plot_df$reg_direction[sig & plot_df$effect < 0] <- "Down"

  legend_label <- paste0("FDR < ", padj_cutoff)
  if (use_effect_threshold && effect_threshold > 0) {
    legend_label <- paste0(legend_label, "\n& |effect| >= ", effect_threshold)
  }

  p <- ggplot(plot_df, aes(x = effect, y = minus_log10_padj, color = reg_direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Up" = "#d7301f", "Down" = "#2c7fb8", "Not_significant" = "grey70")) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      x = x_label,
      y = "-log10 FDR",
      color = legend_label
    )
  if (use_effect_threshold && effect_threshold > 0) {
    p <- p + geom_vline(xintercept = c(-effect_threshold, effect_threshold), linetype = "dashed", color = "grey60", linewidth = 0.3)
  } else {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3)
  }
  ggsave(out_file, p, width = 7, height = 5)
}

plot_volcano_labeled <- function(result_df, out_file, title, x_label, effect_threshold, use_effect_threshold = TRUE, label_n = 15) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (nrow(result_df) == 0 || all(is.na(result_df$padj))) {
    plot_blank(out_file, title, "No plottable CIRI3 DEG results.")
    return()
  }

  plot_df <- result_df
  plot_df$padj[is.na(plot_df$padj)] <- 1
  plot_df$effect[is.na(plot_df$effect)] <- 0
  plot_df$minus_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
  plot_df$reg_direction <- "Not_significant"
  sig <- plot_df$padj < padj_cutoff
  if (use_effect_threshold) {
    sig <- sig & abs(plot_df$effect) >= effect_threshold
  }
  plot_df$reg_direction[sig & plot_df$effect >= 0] <- "Up"
  plot_df$reg_direction[sig & plot_df$effect < 0] <- "Down"
  plot_df$gene_label <- ifelse(
    is.na(plot_df$gene_id) | plot_df$gene_id == "",
    plot_df$circRNA,
    plot_df$gene_id
  )

  label_df <- plot_df[plot_df$reg_direction != "Not_significant", , drop = FALSE]
  label_df <- label_df[order(label_df$padj, -abs(label_df$effect)), , drop = FALSE]
  if (nrow(label_df) > label_n) label_df <- label_df[seq_len(label_n), , drop = FALSE]

  legend_label <- paste0("FDR < ", padj_cutoff)
  if (use_effect_threshold && effect_threshold > 0) {
    legend_label <- paste0(legend_label, "\n& |effect| >= ", effect_threshold)
  }

  p <- ggplot(plot_df, aes(x = effect, y = minus_log10_padj, color = reg_direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_text(
      data = label_df,
      aes(label = gene_label),
      size = 2.8,
      check_overlap = TRUE,
      vjust = -0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("Up" = "#d7301f", "Down" = "#2c7fb8", "Not_significant" = "grey70")) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      x = x_label,
      y = "-log10 FDR",
      color = legend_label
    )
  if (use_effect_threshold && effect_threshold > 0) {
    p <- p + geom_vline(xintercept = c(-effect_threshold, effect_threshold), linetype = "dashed", color = "grey60", linewidth = 0.3)
  } else {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3)
  }
  ggsave(out_file, p, width = 7, height = 5)
}

bsj_path <- input_path("bsj_matrix", input_path("matrix"))
fsj_path <- input_path("fsj_matrix")
ciri3_annotation <- read_ciri3_annotation(input_path("ciri3"))
value_mat <- build_value_matrix(method, bsj_path, fsj_path)
all_group_map <- sample_group_map(groups)

if (mode == "pairwise") {
  result_df <- add_annotation(read_ciri3_result(input_path("result")), ciri3_annotation)
  if (!comparison_name %in% names(comparisons)) {
    stop("Unknown CIRI3 comparison for plotting: ", comparison_name)
  }
  comparison <- comparisons[[comparison_name]]
  case_group <- as.character(comparison[["case_group"]])
  control_group <- as.character(comparison[["control_group"]])
  pair_samples <- c(as.character(unlist(comparison[["case"]])), as.character(unlist(comparison[["control"]])))
  group_map <- all_group_map[pair_samples]
  group_map <- group_map[!is.na(group_map)]

  if (!isTRUE(result_df$effect_from_result[1])) {
    result_df$effect <- matrix_effect(result_df, value_mat, group_map, case_group, control_group)
  }
  effect_from_result <- isTRUE(result_df$effect_from_result[1])
  x_label <- if (effect_from_result) {
    "log2 fold change"
  } else if (method == "de_ratio") {
    "Case-Control mean junction ratio"
  } else {
    "Case-Control mean log2 count"
  }
  effect_threshold <- if (method == "de_ratio" && !effect_from_result) 0 else lfc_cutoff
  use_effect_threshold <- !(method == "de_ratio" && !effect_from_result)

  sig <- rank_significant(result_df, require_effect = use_effect_threshold, effect_threshold = effect_threshold)
  plot_volcano(
    result_df,
    output_path("volcano"),
    paste0(method_label, ": ", case_group, " vs ", control_group),
    x_label,
    effect_threshold,
    use_effect_threshold
  )
  plot_volcano_labeled(
    result_df,
    output_path("volcano_labeled"),
    paste0(method_label, " (labeled): ", case_group, " vs ", control_group),
    x_label,
    effect_threshold,
    use_effect_threshold
  )
  plot_heatmap(
    value_mat,
    result_df,
    sig$circRNA,
    output_path("heatmap"),
    paste0(method_label, " top significant: ", case_group, " vs ", control_group),
    group_map
  )
  plot_pca(
    value_mat,
    output_path("pca"),
    paste0(method_label, " PCA: ", case_group, " vs ", control_group),
    group_map
  )
} else {
  result_paths <- input_paths("results")
  if (length(result_paths) == 0) {
    stop("No CIRI3 pairwise results supplied for all_groups plotting.")
  }
  result_list <- lapply(result_paths, function(path) {
    df <- add_annotation(read_ciri3_result(path), ciri3_annotation)
    df$comparison <- basename(dirname(dirname(path)))
    df
  })
  all_results <- do.call(rbind, result_list)
  all_sig <- all_results[!is.na(all_results$padj) & all_results$padj < padj_cutoff, , drop = FALSE]
  if (nrow(all_sig) > 0) {
    min_padj <- aggregate(padj ~ circRNA, data = all_sig, FUN = min)
    min_padj <- min_padj[order(min_padj$padj), , drop = FALSE]
    heatmap_ids <- min_padj$circRNA
  } else {
    heatmap_ids <- character(0)
  }

  plot_heatmap(
    value_mat,
    all_results,
    heatmap_ids,
    output_path("heatmap"),
    paste0(method_label, " all groups top significant circRNAs"),
    all_group_map
  )
  plot_pca(
    value_mat,
    output_path("pca"),
    paste0(method_label, " all groups PCA"),
    all_group_map
  )
}

message("CIRI3 DEG plotting finished: method=", method, ", mode=", mode)
