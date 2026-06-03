suppressPackageStartupMessages({
  library(ggplot2)
})

# -----------------------------
# Snakemake plumbing
# -----------------------------
get_named <- function(x, name, default = NULL) {
  out <- tryCatch(x[[name]], error = function(e) default)
  if (is.null(out) || length(out) == 0) default else out
}

input_path  <- function(name) as.character(get_named(snakemake@input,  name)[[1]])
output_path <- function(name) as.character(get_named(snakemake@output, name)[[1]])
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
fsj_cpm_file      <- input_path("fsj_cpm_matrix")
totalrna_cpm_file <- input_path("totalrna_cpm_matrix")
metadata_file     <- input_path("metadata")

# Outputs (FSJ branch)
out_fsj_bar_png <- output_path("fsj_barplot_png")
out_fsj_bar_pdf <- output_path("fsj_barplot_pdf")
out_fsj_box_png <- output_path("fsj_boxplot_png")
out_fsj_box_pdf <- output_path("fsj_boxplot_pdf")

# Outputs (totalRNA branch)
out_tot_bar_png <- output_path("totalrna_barplot_png")
out_tot_bar_pdf <- output_path("totalrna_barplot_pdf")
out_tot_box_png <- output_path("totalrna_boxplot_png")
out_tot_box_pdf <- output_path("totalrna_boxplot_pdf")

# Params
condition_order_in  <- as.character(unlist(param("condition_order", character(0))))
condition_colors_in <- param("condition_colors", NULL)
plot_dpi            <- as.numeric(param("plot_dpi", 200))

if (length(condition_order_in) < 1) {
  stop("snakemake.params.condition_order must be a non-empty character vector.")
}
if (is.null(condition_colors_in) || length(condition_colors_in) == 0) {
  stop("snakemake.params.condition_colors must be a named character vector.")
}
condition_order  <- condition_order_in
condition_colors <- setNames(
  as.character(unlist(condition_colors_in)),
  names(condition_colors_in)
)

ensure_dirs <- function(...) {
  for (p in c(...)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
}
ensure_dirs(out_fsj_bar_png, out_fsj_bar_pdf, out_fsj_box_png, out_fsj_box_pdf,
            out_tot_bar_png, out_tot_bar_pdf, out_tot_box_png, out_tot_box_pdf)

# -----------------------------
# Helpers
# -----------------------------
safe_read_tsv <- function(path, header = TRUE, comment.char = "") {
  if (!file.exists(path)) stop("File not found: ", path)
  read.delim(
    path, header = header, sep = "\t",
    check.names = FALSE, quote = "", comment.char = comment.char,
    stringsAsFactors = FALSE
  )
}

# Map user-friendly metadata column names to the canonical names used by the
# rest of the pipeline. No-op when canonical names are already present (as is
# the case for the auto-generated metadata from generate_circrna_metadata.py).
canonicalize_metadata_columns <- function(metadata) {
  colmap <- c(
    sample_name          = "Sample",
    condition            = "Condition",
    biological_replicate = "Rep"
  )
  for (src in names(colmap)) {
    dst <- colmap[[src]]
    if (src %in% colnames(metadata) && !(dst %in% colnames(metadata))) {
      colnames(metadata)[colnames(metadata) == src] <- dst
    }
  }
  metadata
}

prepare_metadata <- function(metadata_file, condition_order) {
  metadata <- safe_read_tsv(metadata_file)
  metadata <- canonicalize_metadata_columns(metadata)

  req <- c("Sample", "Condition")
  miss <- setdiff(req, colnames(metadata))
  if (length(miss) > 0) {
    stop("Metadata missing required columns: ", paste(miss, collapse = ", "))
  }
  metadata$Sample    <- as.character(metadata$Sample)
  metadata$Condition <- trimws(as.character(metadata$Condition))

  unrecognized <- setdiff(unique(metadata$Condition), condition_order)
  if (length(unrecognized) > 0) {
    bad_samples <- metadata$Sample[metadata$Condition %in% unrecognized]
    warning(
      "Condition labels not in condition_order will become NA and be dropped from plots: ",
      paste(unrecognized, collapse = ", "),
      "\nAffected samples: ", paste(bad_samples, collapse = ", ")
    )
  }
  metadata$ConditionLabel <- factor(metadata$Condition,
                                    levels = condition_order, ordered = TRUE)

  rep_col <- intersect(c("Rep", "rep", "Replicate", "replicate", "trial", "Trial"),
                       colnames(metadata))
  if (length(rep_col) > 0) {
    raw_rep <- metadata[[rep_col[1]]]
    rep_num <- suppressWarnings(as.numeric(as.character(raw_rep)))
    if (any(is.na(rep_num) & !is.na(raw_rep))) {
      rep_num <- suppressWarnings(
        as.numeric(gsub("[^0-9]", "", as.character(raw_rep)))
      )
    }
    if (any(is.na(rep_num))) {
      bad <- unique(as.character(raw_rep)[is.na(rep_num)])
      warning("Could not parse replicate numbers from column '", rep_col[1],
              "' for values: ", paste(bad, collapse = ", "),
              "\nFalling back to row-order positions for those samples.")
      rep_num[is.na(rep_num)] <- seq_len(sum(is.na(rep_num)))
    }
    metadata$RepOrder <- rep_num
  } else {
    metadata$RepOrder <- seq_len(nrow(metadata))
  }

  metadata
}

strict_numeric <- function(x, context = "values") {
  was_na <- is.na(x)
  out <- suppressWarnings(as.numeric(as.character(x)))
  newly_na <- is.na(out) & !was_na
  if (any(newly_na)) {
    bad <- unique(as.character(x)[newly_na])
    stop(sprintf("Non-numeric value(s) encountered while parsing %s: %s",
                 context, paste(head(bad, 5), collapse = ", ")))
  }
  out
}

build_replicate_barplot <- function(sample_df, condition_levels, value_col,
                                    ylab_text, title_text, palette_used) {
  sample_df$ConditionLabel <- factor(sample_df$ConditionLabel, levels = condition_levels)
  sample_df <- sample_df[order(sample_df$ConditionLabel, sample_df$RepOrder, sample_df$Sample),
                         , drop = FALSE]
  sample_df$Sample <- factor(sample_df$Sample, levels = sample_df$Sample)

  ggplot(sample_df, aes(x = Sample, y = .data[[value_col]], fill = ConditionLabel)) +
    geom_col(width = 0.75, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = palette_used, drop = FALSE) +
    labs(x = NULL, y = ylab_text, fill = "Condition", title = title_text) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 60, hjust = 1),
      panel.grid.major.x = element_blank()
    )
}

build_group_boxplot <- function(sample_df, condition_levels, value_col,
                                ylab_text, title_text, palette_used) {
  sample_df$ConditionLabel <- factor(sample_df$ConditionLabel, levels = condition_levels)
  ggplot(sample_df, aes(x = ConditionLabel, y = .data[[value_col]], fill = ConditionLabel)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    # height = 0 keeps dot y-positions truthful; only horizontal jitter applied.
    geom_jitter(width = 0.12, height = 0, size = 2, shape = 21, fill = "white") +
    scale_fill_manual(values = palette_used, drop = FALSE) +
    labs(x = "Condition", y = ylab_text, title = title_text) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 30, hjust = 1)
    )
}

save_plot_png_pdf <- function(plot_obj, png_path, pdf_path,
                              width, height, dpi = 200) {
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(pdf_path, plot = plot_obj, width = width, height = height)
}

# Given a CPM matrix DataFrame (first column = id, rest per-sample),
# return per-sample colSums of the numeric block, restricted to samples
# present in `keep_samples` (in keep_samples order).
per_sample_totals <- function(cpm_df, id_cols, keep_samples, context) {
  sample_cols <- setdiff(colnames(cpm_df), id_cols)
  matrix_samples <- intersect(keep_samples, sample_cols)
  missing <- setdiff(keep_samples, matrix_samples)
  if (length(missing) > 0) {
    warning(context, ": these metadata samples are missing from the matrix and will be skipped: ",
            paste(missing, collapse = ", "))
  }
  num_block <- cpm_df[, matrix_samples, drop = FALSE]
  for (nm in matrix_samples) {
    num_block[[nm]] <- strict_numeric(num_block[[nm]],
                                      context = paste0(context, " column '", nm, "'"))
  }
  totals <- colSums(as.matrix(num_block))
  totals[matrix_samples]
}

# -----------------------------
# Read inputs
# -----------------------------
message("FSJ CPM matrix:      ", fsj_cpm_file)
message("totalRNA CPM matrix: ", totalrna_cpm_file)
message("Metadata file:       ", metadata_file)

metadata <- prepare_metadata(metadata_file, condition_order)
# Drop unrecognized-Condition rows so they don't sneak into plots.
metadata <- metadata[!is.na(metadata$ConditionLabel), , drop = FALSE]
if (nrow(metadata) == 0) {
  stop("No samples remain in metadata after dropping unrecognized Condition labels.")
}

fsj_df      <- safe_read_tsv(fsj_cpm_file)
totalrna_df <- safe_read_tsv(totalrna_cpm_file)

# -----------------------------
# Build sample_summary frames and plot
# -----------------------------
plot_branch <- function(label,
                        cpm_df,
                        id_cols,
                        ylab_text,
                        bar_title,
                        box_title,
                        out_bar_png, out_bar_pdf,
                        out_box_png, out_box_pdf,
                        bar_size = c(10, 6),
                        box_size = c(8, 6)) {

  keep_samples <- metadata$Sample
  totals <- per_sample_totals(cpm_df, id_cols, keep_samples, context = label)

  sample_summary <- metadata[match(names(totals), metadata$Sample), , drop = FALSE]
  sample_summary$total_CPM_uniqueMapped <- as.numeric(totals)
  # Re-order rows by Condition x Rep so barplot follows the same convention.
  sample_summary <- sample_summary[order(sample_summary$ConditionLabel,
                                         sample_summary$RepOrder,
                                         sample_summary$Sample), , drop = FALSE]

  palette_used <- condition_colors[condition_order]

  bar <- build_replicate_barplot(sample_summary, condition_order,
                                 "total_CPM_uniqueMapped",
                                 ylab_text, bar_title, palette_used)
  box <- build_group_boxplot(sample_summary, condition_order,
                             "total_CPM_uniqueMapped",
                             ylab_text, box_title, palette_used)

  save_plot_png_pdf(bar, out_bar_png, out_bar_pdf,
                    width = bar_size[1], height = bar_size[2], dpi = plot_dpi)
  save_plot_png_pdf(box, out_box_png, out_box_pdf,
                    width = box_size[1], height = box_size[2], dpi = plot_dpi)

  message(label, ": plotted ", nrow(sample_summary), " sample(s) across ",
          length(unique(as.character(sample_summary$ConditionLabel))), " condition(s).")
}

plot_branch(
  label        = "FSJ",
  cpm_df       = fsj_df,
  id_cols      = "circRNA_ID",
  ylab_text    = "Total FSJ abundance at circRNA loci\n(reads per million STAR-uniquely-mapped)",
  bar_title    = "Sample-level FSJ abundance at circRNA loci normalized by uniquely mapped reads",
  box_title    = "FSJ abundance at circRNA loci across conditions",
  out_bar_png  = out_fsj_bar_png,
  out_bar_pdf  = out_fsj_bar_pdf,
  out_box_png  = out_fsj_box_png,
  out_box_pdf  = out_fsj_box_pdf
)

plot_branch(
  label        = "totalRNA",
  cpm_df       = totalrna_df,
  id_cols      = c("Geneid", "Chr", "Start", "End", "Strand", "Length"),
  ylab_text    = "Total linear gene counts\n(reads per million uniquely mapped)",
  bar_title    = "Sample-level total linear gene counts normalized by uniquely mapped reads",
  box_title    = "Total linear gene counts across conditions",
  out_bar_png  = out_tot_bar_png,
  out_bar_pdf  = out_tot_bar_pdf,
  out_box_png  = out_tot_box_png,
  out_box_pdf  = out_tot_box_pdf
)

message("Linear-counterpart summary plots complete.")
