suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(grid)
})

# -----------------------------
# Snakemake plumbing
# -----------------------------
get_named <- function(x, name, default = NULL) {
  out <- tryCatch(x[[name]], error = function(e) default)
  if (is.null(out) || length(out) == 0) default else out
}

input_path   <- function(name) as.character(get_named(snakemake@input,  name)[[1]])
output_path  <- function(name) as.character(get_named(snakemake@output, name)[[1]])
param <- function(name, default = NULL) {
  val <- get_named(snakemake@params, name, default)
  if (is.null(val)) return(default)
  val
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

options(stringsAsFactors = FALSE, scipen = 999)

# -----------------------------
# Inputs / params
# -----------------------------
bsj_matrix_file <- input_path("bsj_matrix")
metadata_file   <- input_path("metadata")
star_dir        <- as.character(param("star_dir"))
if (is.null(star_dir) || !nzchar(star_dir)) {
  stop("snakemake.params.star_dir is required.")
}

# Order of conditions; row order in metadata also follows this.
condition_order <- as.character(unlist(param("condition_order", character(0))))
if (length(condition_order) < 1) {
  stop("snakemake.params.condition_order must be a non-empty character vector.")
}

# Named character vector: condition -> hex color
condition_colors_in <- param("condition_colors", NULL)
if (is.null(condition_colors_in) || length(condition_colors_in) == 0) {
  stop("snakemake.params.condition_colors must be a named character vector.")
}
condition_colors <- setNames(
  as.character(unlist(condition_colors_in)),
  names(condition_colors_in)
)

# Numeric / behavioral knobs
per_million                          <- as.numeric(param("per_million", 1e6))
heatmap_log_pseudocount              <- as.numeric(param("heatmap_log_pseudocount", 1))
heatmap_top_n_values                 <- as.integer(unlist(param("heatmap_top_n_values", c(50, 25))))
heatmap_cluster_rows                 <- isTRUE(param("heatmap_cluster_rows", TRUE))
heatmap_cluster_cols                 <- isTRUE(param("heatmap_cluster_cols", FALSE))
heatmap_show_rownames_all            <- isTRUE(param("heatmap_show_rownames_all", FALSE))
heatmap_show_rownames_top            <- isTRUE(param("heatmap_show_rownames_top", TRUE))
heatmap_absolute_upper_clip_quantile <- as.numeric(param("heatmap_absolute_upper_clip_quantile", 0.99))
heatmap_absolute_use_log_for_colormap <- isTRUE(param("heatmap_absolute_use_log_for_colormap", FALSE))
heatmap_cellwidth                    <- as.numeric(param("heatmap_cellwidth", 18))
heatmap_cellheight_all               <- as.numeric(param("heatmap_cellheight_all", 2))
heatmap_cellheight_top               <- as.numeric(param("heatmap_cellheight_top", 14))
heatmap_fontsize_col                 <- as.numeric(param("heatmap_fontsize_col", 10))
heatmap_fontsize_row_top             <- as.numeric(param("heatmap_fontsize_row_top", 7))
heatmap_fontsize_row_all             <- as.numeric(param("heatmap_fontsize_row_all", 4))
heatmap_zscore_breaks_limit          <- as.numeric(param("heatmap_zscore_breaks_limit", 3))
plot_dpi                             <- as.numeric(param("plot_dpi", 200))

# Outputs that are fixed (single file paths)
out_sample_summary    <- output_path("sample_summary")
out_bsj_cpm_matrix    <- output_path("bsj_cpm_matrix")
out_condition_summary <- output_path("condition_summary")

out_sample_barplot_png  <- output_path("sample_barplot_png")
out_sample_barplot_pdf  <- output_path("sample_barplot_pdf")
out_conditions_box_png  <- output_path("conditions_boxplot_png")
out_conditions_box_pdf  <- output_path("conditions_boxplot_pdf")

out_all_cpm_matrix      <- output_path("all_cpm_matrix")
out_all_zscore_matrix   <- output_path("all_zscore_matrix")
out_all_cpm_heat_png    <- output_path("all_cpm_heatmap_png")
out_all_cpm_heat_pdf    <- output_path("all_cpm_heatmap_pdf")
out_all_z_heat_png      <- output_path("all_zscore_heatmap_png")
out_all_z_heat_pdf      <- output_path("all_zscore_heatmap_pdf")

# Outputs that are dynamic (one per top_n value).
# We expect the Snakemake rule to declare them via expand() in the same order
# as heatmap_top_n_values, with output names of the form
# `topN_cpm_matrix`, `topN_zscore_matrix`, `topN_cpm_heatmap_png`,
# `topN_cpm_heatmap_pdf`, `topN_zscore_heatmap_png`, `topN_zscore_heatmap_pdf`
# (one for each N). To stay tolerant we resolve these by name lookup but fall
# back to constructing the expected path from the `topN_outdir_template` param
# if needed.
topN_outdir_template <- as.character(param("topN_outdir_template",
                                           file.path(dirname(out_all_cpm_matrix), "..",
                                                     "top{n}_variable")))

resolve_topN_output <- function(top_n, suffix) {
  candidate <- paste0("top", top_n, "_", suffix)
  val <- get_named(snakemake@output, candidate, NULL)
  if (!is.null(val)) return(as.character(val[[1]]))
  # Fall back to template
  base_dir <- gsub("{n}", as.character(top_n), topN_outdir_template, fixed = TRUE)
  ext <- switch(
    suffix,
    cpm_matrix          = "absolute_CPM_matrix.tsv",
    zscore_matrix       = "log2CPM_zscore_matrix.tsv",
    cpm_heatmap_png     = "absolute_CPM_heatmap.png",
    cpm_heatmap_pdf     = "absolute_CPM_heatmap.pdf",
    zscore_heatmap_png  = "zscore_heatmap.png",
    zscore_heatmap_pdf  = "zscore_heatmap.pdf",
    stop("Unknown topN suffix: ", suffix)
  )
  file.path(base_dir, ext)
}

# Ensure all directories exist
ensure_dirs <- function(...) {
  for (p in c(...)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
}
ensure_dirs(out_sample_summary, out_bsj_cpm_matrix, out_condition_summary,
            out_sample_barplot_png, out_sample_barplot_pdf,
            out_conditions_box_png, out_conditions_box_pdf,
            out_all_cpm_matrix, out_all_zscore_matrix,
            out_all_cpm_heat_png, out_all_cpm_heat_pdf,
            out_all_z_heat_png, out_all_z_heat_pdf)

# Palettes (built once, reused everywhere)
abs_palette    <- colorRampPalette(c("#F7FBFF", "#6BAED6", "#08519C"))(101)
zscore_palette <- colorRampPalette(c("#313695", "#FFFFFF", "#A50026"))(101)
zscore_breaks  <- seq(-heatmap_zscore_breaks_limit, heatmap_zscore_breaks_limit,
                      length.out = length(zscore_palette) + 1)

# -----------------------------
# Helpers
# -----------------------------
safe_read_tsv <- function(path, header = TRUE, ...) {
  if (!file.exists(path)) stop("File not found: ", path)
  read.delim(
    path,
    header           = header,
    sep              = "\t",
    check.names      = FALSE,
    quote            = "",
    comment.char     = "",
    stringsAsFactors = FALSE,
    ...
  )
}

write_tsv_base <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

write_matrix_tsv <- function(mat, path) {
  write_tsv_base(
    data.frame(circRNA_ID = rownames(mat), mat, check.names = FALSE),
    path
  )
}

save_plot_png_pdf <- function(plot_obj, png_path, pdf_path, width, height, dpi = 200) {
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(pdf_path, plot = plot_obj, width = width, height = height)
}

extract_uniquely_mapped_reads <- function(log_file) {
  lines <- readLines(log_file, warn = FALSE)
  hit <- grep("^\\s*Uniquely mapped reads number\\s*\\|", lines, value = TRUE)
  if (length(hit) == 0) {
    stop("Could not find 'Uniquely mapped reads number' in: ", log_file)
  }
  value <- sub("^.*\\|\\s*", "", hit[1])
  value <- gsub(",", "", trimws(value))
  as.numeric(value)
}

find_star_log_file <- function(sample_name, star_dir) {
  sample_dir <- file.path(star_dir, sample_name)
  if (!dir.exists(sample_dir)) {
    stop("STAR sample directory not found for sample: ", sample_name,
         "\nExpected: ", sample_dir)
  }
  candidates <- list.files(sample_dir, pattern = "Log\\.final", full.names = TRUE)
  if (length(candidates) == 0) {
    stop("No STAR Log.final file found in: ", sample_dir)
  }
  candidates[1]
}

make_condition_summary <- function(sample_df, condition_levels) {
  if (nrow(sample_df) == 0) return(data.frame())

  tmp <- split(sample_df, sample_df$Condition)
  cond_summary <- do.call(rbind, lapply(tmp, function(d) {
    x <- d$total_BSJ_CPM_uniqueMapped
    data.frame(
      Condition = as.character(d$Condition[1]),
      n         = length(x),
      mean      = mean(x),
      median    = median(x),
      sd        = sd(x),
      se        = sd(x) / sqrt(length(x)),
      stringsAsFactors = FALSE
    )
  }))

  cond_summary <- cond_summary[match(condition_levels, cond_summary$Condition), , drop = FALSE]
  cond_summary <- cond_summary[!is.na(cond_summary$Condition), , drop = FALSE]
  rownames(cond_summary) <- NULL
  cond_summary
}

build_group_boxplot <- function(sample_df, condition_levels, title_text, palette_used) {
  sample_df$Condition <- factor(sample_df$Condition, levels = condition_levels)
  ggplot(sample_df, aes(x = Condition, y = total_BSJ_CPM_uniqueMapped, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 2, shape = 21, fill = "white") +
    scale_fill_manual(values = palette_used, drop = FALSE) +
    labs(
      x     = "Condition",
      y     = "Total circRNA BSJ abundance\n(reads per million uniquely mapped)",
      title = title_text
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 30, hjust = 1)
    )
}

row_zscore_matrix <- function(mat) {
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds   <- apply(mat, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  z <- sweep(mat, 1, row_means, "-")
  z <- sweep(z, 1, row_sds, "/")
  z[is.na(z)] <- 0
  z
}

get_top_variable_rows <- function(mat, top_n) {
  vars <- apply(mat, 1, var, na.rm = TRUE)
  vars[is.na(vars)] <- -Inf
  ord  <- order(vars, decreasing = TRUE)
  keep <- ord[seq_len(min(top_n, nrow(mat)))]
  mat[keep, , drop = FALSE]
}

clip_matrix_upper <- function(mat, q = 0.99) {
  finite_vals <- as.numeric(mat)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  if (length(finite_vals) == 0) return(mat)
  upper <- as.numeric(quantile(finite_vals, probs = q, na.rm = TRUE))
  mat[mat > upper] <- upper
  mat
}

make_col_annotation <- function(metadata_df) {
  data.frame(
    Condition = factor(as.character(metadata_df$Condition), levels = condition_order),
    row.names = metadata_df$Sample,
    stringsAsFactors = FALSE
  )
}

save_pheatmap_png_pdf <- function(mat_display,
                                  png_path,
                                  pdf_path,
                                  main_title,
                                  annotation_col,
                                  annotation_colors,
                                  show_rownames = TRUE,
                                  cellheight    = 10,
                                  cellwidth     = 18,
                                  fontsize_row  = 7,
                                  fontsize_col  = 10,
                                  color_palette = abs_palette,
                                  breaks        = NULL) {

  n_rows <- nrow(mat_display)
  n_cols <- ncol(mat_display)
  if (n_rows < 2 || n_cols < 1) {
    message("  [skip heatmap] ", basename(png_path),
            " (n_rows=", n_rows, ", n_cols=", n_cols, ")")
    return(invisible(NULL))
  }

  width_in  <- max(8,  3   + n_cols * 0.30)
  height_in <- max(6,  2.5 + n_rows * (cellheight / 50))
  height_in <- min(height_in, 60)
  width_in  <- min(width_in,  40)

  draw_one <- function(filename, device) {
    if (device == "png") {
      png(filename, width = width_in, height = height_in, units = "in", res = 200)
    } else {
      pdf(filename, width = width_in, height = height_in)
    }
    on.exit(dev.off(), add = TRUE)

    ph <- pheatmap(
      mat_display,
      color             = color_palette,
      breaks            = breaks,
      cluster_rows      = heatmap_cluster_rows,
      cluster_cols      = heatmap_cluster_cols,
      annotation_col    = annotation_col,
      annotation_colors = annotation_colors,
      show_rownames     = show_rownames,
      show_colnames     = TRUE,
      fontsize_col      = fontsize_col,
      fontsize_row      = fontsize_row,
      cellheight        = cellheight,
      cellwidth         = cellwidth,
      angle_col         = 90,
      border_color      = NA,
      main              = main_title,
      silent            = TRUE
    )
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
  }

  draw_one(png_path, "png")
  draw_one(pdf_path, "pdf")
}

# -----------------------------
# Load & preprocess
# -----------------------------
message("Reading BSJ matrix: ",  bsj_matrix_file)
message("Reading metadata: ",    metadata_file)
message("STAR log directory: ",  star_dir)

bsj      <- safe_read_tsv(bsj_matrix_file)
metadata <- safe_read_tsv(metadata_file)

if (!"circRNA_ID" %in% colnames(bsj)) {
  stop("BSJ matrix must contain a 'circRNA_ID' column: ", bsj_matrix_file)
}
required_meta_cols <- c("Sample", "Condition", "Rep")
missing_meta_cols <- setdiff(required_meta_cols, colnames(metadata))
if (length(missing_meta_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_meta_cols, collapse = ", "))
}

bsj_samples  <- setdiff(colnames(bsj), "circRNA_ID")
meta_samples <- metadata$Sample

missing_in_metadata <- setdiff(bsj_samples, meta_samples)
if (length(missing_in_metadata) > 0) {
  stop("These BSJ samples are missing from metadata:\n",
       paste(missing_in_metadata, collapse = "\n"))
}
metadata <- metadata[metadata$Sample %in% bsj_samples, , drop = FALSE]
metadata <- metadata[match(bsj_samples, metadata$Sample), , drop = FALSE]

# STAR uniquely-mapped reads per sample
metadata$star_log_file <- vapply(
  metadata$Sample, FUN = find_star_log_file,
  FUN.VALUE = character(1), star_dir = star_dir
)
metadata$uniquely_mapped_reads <- vapply(
  metadata$star_log_file, FUN = extract_uniquely_mapped_reads,
  FUN.VALUE = numeric(1)
)
if (any(is.na(metadata$uniquely_mapped_reads) | metadata$uniquely_mapped_reads <= 0)) {
  stop("Some uniquely mapped read counts are missing or non-positive.")
}

# Validate / order conditions
unknown_conds <- setdiff(unique(metadata$Condition), condition_order)
if (length(unknown_conds) > 0) {
  stop("Metadata contains conditions not in condition_order: ",
       paste(unknown_conds, collapse = ", "))
}
metadata$Condition <- factor(metadata$Condition, levels = condition_order, ordered = TRUE)

metadata$RepOrder <- suppressWarnings(as.numeric(metadata$Rep))
if (any(is.na(metadata$RepOrder))) {
  metadata$RepOrder[is.na(metadata$RepOrder)] <- seq_len(sum(is.na(metadata$RepOrder)))
}
metadata <- metadata[order(metadata$Condition, metadata$RepOrder, metadata$Sample),
                     , drop = FALSE]

# Numeric BSJ matrix with sample order matched to metadata
bsj_numeric <- bsj[, bsj_samples, drop = FALSE]
for (j in seq_along(bsj_numeric)) bsj_numeric[[j]] <- as.numeric(bsj_numeric[[j]])
rownames(bsj_numeric) <- bsj$circRNA_ID

sample_order <- metadata$Sample
bsj_numeric  <- bsj_numeric[, sample_order, drop = FALSE]

unique_map_vec <- metadata$uniquely_mapped_reads
names(unique_map_vec) <- metadata$Sample

bsj_cpm <- sweep(as.matrix(bsj_numeric), 2,
                 unique_map_vec[colnames(bsj_numeric)], "/") * per_million
bsj_cpm_df <- data.frame(circRNA_ID = rownames(bsj_cpm), bsj_cpm, check.names = FALSE)

total_bsj_counts <- colSums(bsj_numeric, na.rm = TRUE)
total_bsj_cpm    <- (total_bsj_counts / unique_map_vec[names(total_bsj_counts)]) * per_million

sample_summary <- metadata
sample_summary$total_BSJ_counts          <- total_bsj_counts[metadata$Sample]
sample_summary$total_BSJ_CPM_uniqueMapped <- total_bsj_cpm[metadata$Sample]

# -----------------------------
# Part 2: normalized BSJ tables + plots
# -----------------------------
write_tsv_base(sample_summary, out_sample_summary)
write_tsv_base(bsj_cpm_df,     out_bsj_cpm_matrix)

overall_condition_summary <- make_condition_summary(sample_summary, condition_order)
write_tsv_base(overall_condition_summary, out_condition_summary)

# All-sample barplot
sample_summary$Sample <- factor(sample_summary$Sample, levels = sample_summary$Sample)
present_conditions <- unique(as.character(sample_summary$Condition))
present_conditions <- present_conditions[!is.na(present_conditions)]
palette_used_root  <- condition_colors[present_conditions]

p_sample <- ggplot(sample_summary,
                   aes(x = Sample, y = total_BSJ_CPM_uniqueMapped, fill = Condition)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = palette_used_root, drop = FALSE) +
  labs(
    x     = NULL,
    y     = "Total circRNA BSJ abundance\n(reads per million uniquely mapped)",
    fill  = "Condition",
    title = "Sample-level circRNA BSJ abundance normalized by uniquely mapped reads"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 60, hjust = 1),
    panel.grid.major.x = element_blank()
  )

save_plot_png_pdf(p_sample, out_sample_barplot_png, out_sample_barplot_pdf,
                  width = 12, height = 6, dpi = plot_dpi)

# All-condition boxplot
p_root_box <- build_group_boxplot(
  sample_summary, condition_order,
  title_text   = "All conditions with all replicates",
  palette_used = condition_colors[condition_order]
)
save_plot_png_pdf(p_root_box, out_conditions_box_png, out_conditions_box_pdf,
                  width = 9, height = 6, dpi = plot_dpi)

# -----------------------------
# Part 3: heatmaps
# -----------------------------
metadata_for_cols <- metadata[match(colnames(bsj_cpm), metadata$Sample), , drop = FALSE]
annotation_col    <- make_col_annotation(metadata_for_cols)
annotation_colors <- list(Condition = condition_colors)

# All circRNAs matrices
abs_mat_all_raw <- bsj_cpm
abs_mat_all_display <- if (heatmap_absolute_use_log_for_colormap) {
  log2(abs_mat_all_raw + heatmap_log_pseudocount)
} else {
  abs_mat_all_raw
}
abs_mat_all_display <- clip_matrix_upper(abs_mat_all_display,
                                         q = heatmap_absolute_upper_clip_quantile)

log_mat_all <- log2(abs_mat_all_raw + heatmap_log_pseudocount)
z_mat_all   <- row_zscore_matrix(log_mat_all)

write_matrix_tsv(abs_mat_all_raw, out_all_cpm_matrix)
write_matrix_tsv(z_mat_all,       out_all_zscore_matrix)

save_pheatmap_png_pdf(
  abs_mat_all_display,
  png_path          = out_all_cpm_heat_png,
  pdf_path          = out_all_cpm_heat_pdf,
  main_title        = "All circRNAs (no z-score)",
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames     = heatmap_show_rownames_all,
  cellheight        = heatmap_cellheight_all,
  cellwidth         = heatmap_cellwidth,
  fontsize_row      = heatmap_fontsize_row_all,
  fontsize_col      = heatmap_fontsize_col,
  color_palette     = abs_palette
)
save_pheatmap_png_pdf(
  z_mat_all,
  png_path          = out_all_z_heat_png,
  pdf_path          = out_all_z_heat_pdf,
  main_title        = "All circRNAs (z-score)",
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames     = heatmap_show_rownames_all,
  cellheight        = heatmap_cellheight_all,
  cellwidth         = heatmap_cellwidth,
  fontsize_row      = heatmap_fontsize_row_all,
  fontsize_col      = heatmap_fontsize_col,
  color_palette     = zscore_palette,
  breaks            = zscore_breaks
)

# Top-N variable heatmaps
for (top_n in heatmap_top_n_values) {
  top_log        <- get_top_variable_rows(log_mat_all, top_n)
  abs_subset_raw <- abs_mat_all_raw[rownames(top_log), , drop = FALSE]

  abs_subset_display <- if (heatmap_absolute_use_log_for_colormap) {
    log2(abs_subset_raw + heatmap_log_pseudocount)
  } else {
    abs_subset_raw
  }
  abs_subset_display <- clip_matrix_upper(abs_subset_display,
                                          q = heatmap_absolute_upper_clip_quantile)

  z_subset <- row_zscore_matrix(log2(abs_subset_raw + heatmap_log_pseudocount))

  out_top_cpm_mat    <- resolve_topN_output(top_n, "cpm_matrix")
  out_top_z_mat      <- resolve_topN_output(top_n, "zscore_matrix")
  out_top_cpm_png    <- resolve_topN_output(top_n, "cpm_heatmap_png")
  out_top_cpm_pdf    <- resolve_topN_output(top_n, "cpm_heatmap_pdf")
  out_top_z_png      <- resolve_topN_output(top_n, "zscore_heatmap_png")
  out_top_z_pdf      <- resolve_topN_output(top_n, "zscore_heatmap_pdf")

  ensure_dirs(out_top_cpm_mat, out_top_z_mat,
              out_top_cpm_png, out_top_cpm_pdf,
              out_top_z_png,   out_top_z_pdf)

  write_matrix_tsv(abs_subset_raw, out_top_cpm_mat)
  write_matrix_tsv(z_subset,       out_top_z_mat)

  save_pheatmap_png_pdf(
    abs_subset_display,
    png_path          = out_top_cpm_png,
    pdf_path          = out_top_cpm_pdf,
    main_title        = paste0("Top ", top_n, " (no z-score)"),
    annotation_col    = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames     = heatmap_show_rownames_top,
    cellheight        = heatmap_cellheight_top,
    cellwidth         = heatmap_cellwidth,
    fontsize_row      = heatmap_fontsize_row_top,
    fontsize_col      = heatmap_fontsize_col,
    color_palette     = abs_palette
  )

  save_pheatmap_png_pdf(
    z_subset,
    png_path          = out_top_z_png,
    pdf_path          = out_top_z_pdf,
    main_title        = paste0("Top ", top_n, " (z-score)"),
    annotation_col    = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames     = heatmap_show_rownames_top,
    cellheight        = heatmap_cellheight_top,
    cellwidth         = heatmap_cellwidth,
    fontsize_row      = heatmap_fontsize_row_top,
    fontsize_col      = heatmap_fontsize_col,
    color_palette     = zscore_palette,
    breaks            = zscore_breaks
  )
}

message("Dataset summary complete: ",
        nrow(bsj_cpm), " circRNAs across ", ncol(bsj_cpm), " samples / ",
        length(condition_order), " conditions.")
