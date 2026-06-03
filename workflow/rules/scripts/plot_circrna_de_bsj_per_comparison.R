suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# -----------------------------
# Snakemake plumbing
# -----------------------------
get_named <- function(x, name, default = NULL) {
  out <- tryCatch(x[[name]], error = function(e) default)
  if (is.null(out) || length(out) == 0) default else out
}

input_path <- function(name) as.character(get_named(snakemake@input,  name)[[1]])
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
result_file        <- input_path("result")
infor_file         <- input_path("infor")

comparison_name    <- as.character(param("comparison", "comparison"))
case_label         <- as.character(param("case_label", "Case"))
control_label      <- as.character(param("control_label", "Control"))
fdr_cutoff         <- as.numeric(param("fdr_cutoff", 0.05))
logfc_cutoff       <- as.numeric(param("logfc_cutoff", 1))
label_top_n        <- as.integer(param("label_top_n", 15))
reverse_lfc        <- isTRUE(param("reverse_ciri3_logfc_direction", TRUE))
fdr_floor_for_plot <- as.numeric(param("fdr_floor_for_plot", 1e-300))

# Plotting style (sensible defaults; rarely need overriding)
ma_width   <- as.numeric(param("ma_width",   12))
ma_height  <- as.numeric(param("ma_height",   9))
vol_width  <- as.numeric(param("vol_width",  12))
vol_height <- as.numeric(param("vol_height",  9))
plot_dpi   <- as.numeric(param("plot_dpi",  200))

col_case    <- as.character(param("col_case",    "#D95F02"))
col_control <- as.character(param("col_control", "#0072B2"))
col_nonsig  <- as.character(param("col_nonsig",  "#B3B3B3"))

# Outputs
out_cleaned_tsv        <- output_path("cleaned")
out_candidates_tsv     <- output_path("candidates")
out_ma_unlabeled_png   <- output_path("ma_unlabeled_png")
out_ma_unlabeled_pdf   <- output_path("ma_unlabeled_pdf")
out_ma_labeled_png     <- output_path("ma_labeled_png")
out_ma_labeled_pdf     <- output_path("ma_labeled_pdf")
out_vol_unlabeled_png  <- output_path("volcano_unlabeled_png")
out_vol_unlabeled_pdf  <- output_path("volcano_unlabeled_pdf")
out_vol_labeled_png    <- output_path("volcano_labeled_png")
out_vol_labeled_pdf    <- output_path("volcano_labeled_pdf")

for (p in c(out_cleaned_tsv, out_candidates_tsv,
            out_ma_unlabeled_png, out_ma_unlabeled_pdf,
            out_ma_labeled_png,   out_ma_labeled_pdf,
            out_vol_unlabeled_png, out_vol_unlabeled_pdf,
            out_vol_labeled_png,   out_vol_labeled_pdf)) {
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
}

message("CIRI3 DE_BSJ plotting for comparison: ", comparison_name)
message("  case_label = ", case_label, " | control_label = ", control_label)
message("  FDR cutoff = ", fdr_cutoff, " | |logFC| cutoff = ", logfc_cutoff,
        " | reverse_ciri3_logfc_direction = ", reverse_lfc)

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

save_plot_png_pdf <- function(plot_obj, png_path, pdf_path, width, height, dpi = 200) {
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(pdf_path, plot = plot_obj, width = width, height = height)
}

read_ciri3_result <- function(path) {
  raw <- safe_read_tsv(path, header = TRUE)

  if (!"circRNA_ID" %in% colnames(raw)) {
    if (ncol(raw) == 6) {
      rn <- rownames(raw)
      if (is.null(rn) || length(rn) != nrow(raw)) {
        stop("Could not recover circRNA_ID values from row names in: ", path)
      }
      default_rn <- as.character(seq_len(nrow(raw)))
      if (all(rn == default_rn)) {
        stop("result.txt has 6 columns but row names do not contain circRNA IDs: ", path)
      }
      raw <- data.frame(circRNA_ID = rn, raw, check.names = FALSE,
                        stringsAsFactors = FALSE)
      rownames(raw) <- NULL
    } else if (ncol(raw) == 7) {
      colnames(raw) <- c("circRNA_ID", "logFC", "logCPM", "LR", "PValue", "DE", "FDR")
    } else {
      stop("Unexpected number of columns in result.txt: ", ncol(raw),
           ". Expected 6 (with circRNA_ID in row names) or 7.")
    }
  }

  required <- c("circRNA_ID", "logFC", "logCPM", "LR", "PValue", "DE", "FDR")
  missing_cols <- setdiff(required, colnames(raw))
  if (length(missing_cols) > 0) {
    stop("Missing expected columns in result.txt after repair: ",
         paste(missing_cols, collapse = ", "))
  }

  raw <- raw[, required, drop = FALSE]
  num_cols <- c("logFC", "logCPM", "LR", "PValue", "DE", "FDR")
  for (nm in num_cols) raw[[nm]] <- suppressWarnings(as.numeric(raw[[nm]]))
  raw$circRNA_ID <- as.character(raw$circRNA_ID)
  raw
}

classify_points <- function(df, case_label, control_label, fdr_cutoff, logfc_cutoff) {
  df$direction_class <- "Not significant"
  sig <- !is.na(df$FDR) & !is.na(df$logFC_aligned) &
    df$FDR < fdr_cutoff & abs(df$logFC_aligned) >= logfc_cutoff

  df$direction_class[sig & df$logFC_aligned > 0] <- paste0("Higher in ", case_label)
  df$direction_class[sig & df$logFC_aligned < 0] <- paste0("Higher in ", control_label)

  df$direction_class <- factor(
    df$direction_class,
    levels = c(
      paste0("Higher in ", case_label),
      paste0("Higher in ", control_label),
      "Not significant"
    )
  )
  df
}

choose_labels <- function(df, top_n = 15) {
  sig_df <- df[!is.na(df$FDR) & !is.na(df$logFC_aligned) &
                 df$direction_class != "Not significant", , drop = FALSE]
  if (nrow(sig_df) == 0) return(df[0, , drop = FALSE])

  sig_df$abs_logFC_aligned <- abs(sig_df$logFC_aligned)
  sig_df <- sig_df[order(sig_df$FDR, -sig_df$abs_logFC_aligned, -sig_df$neg_log10_FDR_plot),
                   , drop = FALSE]
  head(sig_df, top_n)
}

de_plot_theme <- function() {
  theme_bw(base_size = 16) +
    theme(
      plot.title       = element_text(face = "bold", size = 20),
      plot.subtitle    = element_text(size = 14),
      legend.position  = "top",
      legend.text      = element_text(size = 12),
      panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.8),
      panel.grid.minor = element_line(color = "#E6E6E6", linewidth = 0.5),
      axis.title       = element_text(size = 16),
      axis.text        = element_text(size = 13)
    )
}

de_color_values <- function(case_label, control_label) {
  c(
    setNames(col_case,    paste0("Higher in ", case_label)),
    setNames(col_control, paste0("Higher in ", control_label)),
    "Not significant" = col_nonsig
  )
}

add_repel_labels <- function(p, labeled_df) {
  if (is.null(labeled_df) || nrow(labeled_df) == 0) return(p)
  p + geom_text_repel(
    data          = labeled_df,
    aes(label     = circRNA_ID),
    size          = 4.3,
    max.overlaps  = Inf,
    box.padding   = 0.35,
    point.padding = 0.25,
    segment.color = "grey50",
    show.legend   = FALSE
  )
}

build_ma_plot <- function(df, title_text, subtitle_text, case_label, control_label,
                          labeled_df = NULL) {
  p <- ggplot(df, aes(x = logCPM, y = logFC_aligned, color = direction_class)) +
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
    geom_hline(yintercept = c(-logfc_cutoff, logfc_cutoff),
               linetype = "dashed", linewidth = 0.8, color = "black") +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_color_manual(values = de_color_values(case_label, control_label), drop = FALSE) +
    labs(
      title    = title_text,
      subtitle = subtitle_text,
      x        = "Average abundance (logCPM)",
      y        = "CIRI3_DE logFC_aligned",
      color    = NULL
    ) +
    de_plot_theme()
  add_repel_labels(p, labeled_df)
}

build_volcano_plot <- function(df, title_text, subtitle_text, case_label, control_label,
                               labeled_df = NULL) {
  p <- ggplot(df, aes(x = logFC_aligned, y = neg_log10_FDR_plot, color = direction_class)) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
               linetype = "dashed", linewidth = 0.8, color = "black") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", linewidth = 0.8, color = "black") +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_color_manual(values = de_color_values(case_label, control_label), drop = FALSE) +
    labs(
      title    = title_text,
      subtitle = subtitle_text,
      x        = "CIRI3_DE logFC_aligned",
      y        = expression(-log[10](FDR)),
      color    = NULL
    ) +
    de_plot_theme()
  add_repel_labels(p, labeled_df)
}

# -----------------------------
# Main
# -----------------------------
res   <- read_ciri3_result(result_file)
infor <- safe_read_tsv(infor_file)

if (!all(c("Sample", "Class") %in% colnames(infor))) {
  stop("infor.tsv must contain Sample and Class columns: ", infor_file)
}

# Sanity: check that infor Case/Control samples align with what the pipeline thinks
infor_case_samples    <- infor$Sample[tolower(infor$Class) == "case"]
infor_control_samples <- infor$Sample[tolower(infor$Class) == "control"]
message("  infor Case samples (n=", length(infor_case_samples), "): ",
        paste(infor_case_samples, collapse = ", "))
message("  infor Control samples (n=", length(infor_control_samples), "): ",
        paste(infor_control_samples, collapse = ", "))

# Align logFC sign for plotting/classification
res$logFC_raw     <- res$logFC
res$logFC_aligned <- if (reverse_lfc) -res$logFC_raw else res$logFC_raw

# Classify points and prepare plotting columns
res <- classify_points(res, case_label, control_label, fdr_cutoff, logfc_cutoff)
res$comparison    <- comparison_name
res$case_label    <- case_label
res$control_label <- control_label
res$abs_logFC     <- abs(res$logFC_aligned)

res$FDR_plot <- res$FDR
res$FDR_plot[is.na(res$FDR_plot)] <- 1
res$FDR_plot[res$FDR_plot <= 0]   <- fdr_floor_for_plot
res$neg_log10_FDR_plot <- -log10(res$FDR_plot)

labeled_df <- choose_labels(res, top_n = label_top_n)

# Cleaned results table
cleaned_out <- res[, c("comparison", "circRNA_ID", "logFC_raw", "logFC_aligned",
                       "logCPM", "LR", "PValue", "DE", "FDR",
                       "abs_logFC", "neg_log10_FDR_plot", "direction_class",
                       "case_label", "control_label")]
write_tsv_base(cleaned_out, out_cleaned_tsv)

# Significant-candidates table (ranked)
sig_out <- cleaned_out[cleaned_out$direction_class != "Not significant", , drop = FALSE]
if (nrow(sig_out) > 0) {
  sig_out <- sig_out[order(sig_out$FDR, -abs(sig_out$logFC_aligned)), , drop = FALSE]
  sig_out$candidate_rank <- seq_len(nrow(sig_out))
} else {
  sig_out$candidate_rank <- integer(0)
}
write_tsv_base(sig_out, out_candidates_tsv)

# Plots
title_ma  <- paste0(comparison_name, " CIRI3_DE DE_BSJ MA Plot")
title_vol <- paste0(comparison_name, " CIRI3_DE DE_BSJ Volcano Plot")
subtitle  <- paste0(
  "Case group: ", case_label,
  " | Control group: ", control_label,
  " | Thresholds: FDR < ", fdr_cutoff, ", |logFC| >= ", logfc_cutoff
)

ma_unlabeled  <- build_ma_plot(res,      title_ma,  subtitle, case_label, control_label, NULL)
ma_labeled    <- build_ma_plot(res,      title_ma,  subtitle, case_label, control_label, labeled_df)
vol_unlabeled <- build_volcano_plot(res, title_vol, subtitle, case_label, control_label, NULL)
vol_labeled   <- build_volcano_plot(res, title_vol, subtitle, case_label, control_label, labeled_df)

save_plot_png_pdf(ma_unlabeled,  out_ma_unlabeled_png,  out_ma_unlabeled_pdf,
                  width = ma_width,  height = ma_height,  dpi = plot_dpi)
save_plot_png_pdf(ma_labeled,    out_ma_labeled_png,    out_ma_labeled_pdf,
                  width = ma_width,  height = ma_height,  dpi = plot_dpi)
save_plot_png_pdf(vol_unlabeled, out_vol_unlabeled_png, out_vol_unlabeled_pdf,
                  width = vol_width, height = vol_height, dpi = plot_dpi)
save_plot_png_pdf(vol_labeled,   out_vol_labeled_png,   out_vol_labeled_pdf,
                  width = vol_width, height = vol_height, dpi = plot_dpi)

message("Done. Wrote ", nrow(res), " circRNAs total, ",
        nrow(sig_out), " significant after FDR<", fdr_cutoff,
        " and |logFC|>=", logfc_cutoff, ".")
