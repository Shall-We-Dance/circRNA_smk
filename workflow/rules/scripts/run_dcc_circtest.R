suppressPackageStartupMessages({
  library(ggplot2)
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

if (!requireNamespace("CircTest", quietly = TRUE)) {
  stop(
    "The R package 'CircTest' is required for dcc.run_circtest. ",
    "Install circtools/CircTest in workflow/rules/envs/circtest.yaml or disable dcc.run_circtest."
  )
}
suppressPackageStartupMessages(library(CircTest))

list_param <- function(value) {
  as.character(unlist(value, use.names = FALSE))
}

numeric_param <- function(value, label, min_value = NULL, max_value = NULL) {
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed)) {
    stop(label, " must be numeric.")
  }
  if (!is.null(min_value) && parsed < min_value) {
    stop(label, " must be >= ", min_value, ".")
  }
  if (!is.null(max_value) && parsed > max_value) {
    stop(label, " must be <= ", max_value, ".")
  }
  parsed
}

int_param <- function(value, label, min_value = NULL) {
  as.integer(round(numeric_param(value, label, min_value = min_value)))
}

all_samples <- list_param(snakemake@params[["all_samples"]])
selected_samples <- list_param(snakemake@params[["selected_samples"]])
group_labels <- list_param(snakemake@params[["group_labels"]])
comparison <- as.character(snakemake@params[["comparison"]])
case_group <- as.character(snakemake@params[["case_group"]])
control_group <- as.character(snakemake@params[["control_group"]])
padj_cutoff <- numeric_param(snakemake@params[["padj_cutoff"]], "deg.padj_cutoff", 0, 1)
filter_sample <- int_param(snakemake@params[["filter_sample"]], "dcc.circtest_filter_sample", 1)
filter_count <- int_param(snakemake@params[["filter_count"]], "dcc.circtest_filter_count", 0)
percentage <- numeric_param(snakemake@params[["percentage"]], "dcc.circtest_percentage", 0)
max_plots <- int_param(snakemake@params[["max_plots"]], "dcc.circtest_max_plots", 0)

if (length(selected_samples) < 2) {
  stop("DCC CircTest requires at least two selected samples.")
}
if (length(selected_samples) != length(group_labels)) {
  stop("DCC CircTest selected sample and group label counts differ.")
}
if (!all(selected_samples %in% all_samples)) {
  missing <- setdiff(selected_samples, all_samples)
  stop("DCC CircTest selected samples missing from DCC sample order: ", paste(missing, collapse = ", "))
}
if (length(unique(group_labels)) < 2) {
  stop("DCC CircTest requires at least two groups.")
}

message("Starting DCC CircTest analysis for ", comparison, ".")
message("Selected samples: ", paste(selected_samples, collapse = ", "))
message("Group labels: ", paste(group_labels, collapse = ", "))

looks_like_header <- function(row) {
  if (length(row) < 3) {
    return(TRUE)
  }
  start <- suppressWarnings(as.numeric(row[[2]]))
  end <- suppressWarnings(as.numeric(row[[3]]))
  is.na(start) || is.na(end)
}

read_tsv_rows <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(list())
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    return(list())
  }
  strsplit(lines, "\t", fixed = TRUE)
}

normalize_int <- function(value) {
  parsed <- suppressWarnings(as.numeric(value))
  parsed[!is.finite(parsed) | parsed < 0] <- 0
  as.integer(round(parsed))
}

read_dcc_counts <- function(path, sample_order) {
  rows <- read_tsv_rows(path)
  if (length(rows) == 0) {
    out <- data.frame(
      circRNA = character(0),
      chrom = character(0),
      start = integer(0),
      end = integer(0),
      stringsAsFactors = FALSE
    )
    for (sample in sample_order) {
      out[[sample]] <- integer(0)
    }
    return(out)
  }

  header <- looks_like_header(rows[[1]])
  data_rows <- if (header) rows[-1] else rows
  if (length(data_rows) == 0) {
    out <- data.frame(
      circRNA = character(0),
      chrom = character(0),
      start = integer(0),
      end = integer(0),
      stringsAsFactors = FALSE
    )
    for (sample in sample_order) {
      out[[sample]] <- integer(0)
    }
    return(out)
  }

  required_cols <- 3 + length(sample_order)
  padded <- lapply(data_rows, function(row) {
    length(row) <- max(length(row), required_cols)
    row[seq_len(required_cols)]
  })
  mat <- do.call(rbind, padded)
  df <- data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- c("chrom", "start", "end", sample_order)
  df$start <- suppressWarnings(as.integer(as.numeric(df$start)))
  df$end <- suppressWarnings(as.integer(as.numeric(df$end)))
  keep <- !is.na(df$start) & !is.na(df$end) & nzchar(trimws(df$chrom))
  df <- df[keep, , drop = FALSE]
  df$chrom <- trimws(df$chrom)
  df$circRNA <- paste0(df$chrom, ":", df$start, "|", df$end)
  for (sample in sample_order) {
    df[[sample]] <- normalize_int(df[[sample]])
  }
  df[, c("circRNA", "chrom", "start", "end", sample_order), drop = FALSE]
}

read_dcc_coordinates <- function(path) {
  rows <- read_tsv_rows(path)
  if (length(rows) == 0) {
    return(data.frame(
      circRNA = character(0),
      gene_id = character(0),
      strand = character(0),
      junction_type = character(0),
      stringsAsFactors = FALSE
    ))
  }
  default_header <- c(
    "chrom", "start", "end", "gene_id", "junction_type",
    "strand", "circ_region", "overall_region"
  )
  has_header <- looks_like_header(rows[[1]])
  header <- if (has_header) rows[[1]] else default_header
  data_rows <- if (has_header) rows[-1] else rows
  normalized_header <- tolower(trimws(header))

  pick_col <- function(candidates, fallback) {
    idx <- match(tolower(candidates), normalized_header, nomatch = 0)
    idx <- idx[idx > 0]
    if (length(idx) > 0) idx[[1]] else fallback
  }

  chrom_col <- pick_col(c("chrom", "chr", "#chrom"), 1)
  start_col <- pick_col(c("start", "chromstart"), 2)
  end_col <- pick_col(c("end", "chromend"), 3)
  gene_col <- pick_col(c("gene_id", "genename", "gene", "gene_name"), 4)
  junction_col <- pick_col(c("junction_type", "junctiontype"), 5)
  strand_col <- pick_col(c("strand"), 6)

  read_cell <- function(row, idx, default = "") {
    if (!is.na(idx) && idx <= length(row)) {
      return(trimws(row[[idx]]))
    }
    default
  }

  out <- lapply(data_rows, function(row) {
    chrom <- read_cell(row, chrom_col)
    start <- suppressWarnings(as.integer(as.numeric(read_cell(row, start_col))))
    end <- suppressWarnings(as.integer(as.numeric(read_cell(row, end_col))))
    if (!nzchar(chrom) || is.na(start) || is.na(end)) {
      return(NULL)
    }
    data.frame(
      circRNA = paste0(chrom, ":", start, "|", end),
      gene_id = read_cell(row, gene_col, NA_character_),
      strand = read_cell(row, strand_col, NA_character_),
      junction_type = read_cell(row, junction_col, NA_character_),
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) {
    return(data.frame(
      circRNA = character(0),
      gene_id = character(0),
      strand = character(0),
      junction_type = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}

write_table <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_blank_pdf <- function(path, message) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pdf(path, width = 7, height = 5)
  plot.new()
  title(main = message)
  dev.off()
}

empty_results <- function() {
  data.frame(
    circRNA = character(0),
    chrom = character(0),
    start = integer(0),
    end = integer(0),
    gene_id = character(0),
    strand = character(0),
    junction_type = character(0),
    pvalue = numeric(0),
    padj = numeric(0),
    regulation = character(0),
    stringsAsFactors = FALSE
  )
}

write_summary <- function(input_n, filtered_n, tested_n, significant_n, status) {
  summary <- data.frame(
    comparison = comparison,
    case_group = case_group,
    control_group = control_group,
    selected_samples = paste(selected_samples, collapse = ","),
    group_labels = paste(group_labels, collapse = ","),
    input_circRNAs = input_n,
    filtered_circRNAs = filtered_n,
    tested_circRNAs = tested_n,
    significant_circRNAs = significant_n,
    padj_cutoff = padj_cutoff,
    filter_sample = filter_sample,
    filter_count = filter_count,
    percentage = percentage,
    status = status,
    stringsAsFactors = FALSE
  )
  write_table(summary, snakemake@output[["summary"]])
}

make_circtest_tables <- function(circ_counts, linear_counts, samples) {
  original_n <- nrow(circ_counts)
  common_ids <- intersect(circ_counts$circRNA, linear_counts$circRNA)
  circ_counts <- circ_counts[match(common_ids, circ_counts$circRNA), , drop = FALSE]
  linear_counts <- linear_counts[match(common_ids, linear_counts$circRNA), , drop = FALSE]

  if (length(common_ids) == 0) {
    Circ <- data.frame(circRNA = character(0), stringsAsFactors = FALSE)
    Linear <- data.frame(circRNA = character(0), stringsAsFactors = FALSE)
    for (sample in samples) {
      Circ[[sample]] <- integer(0)
      Linear[[sample]] <- integer(0)
    }
    return(list(Circ = Circ, Linear = Linear, dropped = original_n))
  }

  Circ <- data.frame(circRNA = circ_counts$circRNA, circ_counts[, samples, drop = FALSE], check.names = FALSE)
  Linear <- data.frame(circRNA = linear_counts$circRNA, linear_counts[, samples, drop = FALSE], check.names = FALSE)
  rownames(Circ) <- Circ$circRNA
  rownames(Linear) <- Linear$circRNA
  list(Circ = Circ, Linear = Linear, dropped = original_n - length(common_ids))
}

simple_filter <- function(Circ, Linear) {
  count_cols <- setdiff(colnames(Circ), "circRNA")
  circ_mat <- as.matrix(Circ[, count_cols, drop = FALSE])
  linear_mat <- as.matrix(Linear[, count_cols, drop = FALSE])
  storage.mode(circ_mat) <- "numeric"
  storage.mode(linear_mat) <- "numeric"
  denominator <- circ_mat + linear_mat
  ratio <- circ_mat / denominator
  ratio[!is.finite(ratio)] <- 0
  keep <- rowSums(circ_mat >= filter_count) >= filter_sample &
    apply(ratio, 1, max, na.rm = TRUE) >= percentage
  keep[is.na(keep)] <- FALSE
  Circ[keep, , drop = FALSE]
}

plot_ratios <- function(Circ, Linear, result_df, out_pdf) {
  sig_df <- result_df[!is.na(result_df$padj) & result_df$padj < padj_cutoff, , drop = FALSE]
  sig_df <- sig_df[order(sig_df$padj), , drop = FALSE]
  if (nrow(sig_df) == 0 || max_plots == 0) {
    write_blank_pdf(out_pdf, "No significant DCC CircTest circRNAs to plot.")
    return()
  }
  sig_ids <- head(sig_df$circRNA, max_plots)
  sig_ids <- sig_ids[sig_ids %in% rownames(Circ)]
  if (length(sig_ids) == 0) {
    write_blank_pdf(out_pdf, "No significant DCC CircTest circRNAs to plot.")
    return()
  }

  count_cols <- selected_samples
  circ_mat <- as.matrix(Circ[sig_ids, count_cols, drop = FALSE])
  linear_mat <- as.matrix(Linear[sig_ids, count_cols, drop = FALSE])
  storage.mode(circ_mat) <- "numeric"
  storage.mode(linear_mat) <- "numeric"
  ratio <- circ_mat / (circ_mat + linear_mat)
  ratio[!is.finite(ratio)] <- 0

  plot_df <- data.frame(
    circRNA = rep(sig_ids, times = length(count_cols)),
    sample = rep(count_cols, each = length(sig_ids)),
    group = rep(group_labels, each = length(sig_ids)),
    ratio = as.vector(ratio),
    circ_count = as.vector(circ_mat),
    linear_count = as.vector(linear_mat),
    stringsAsFactors = FALSE
  )
  plot_df$circRNA <- factor(plot_df$circRNA, levels = sig_ids)
  plot_df$group <- factor(plot_df$group, levels = unique(group_labels))

  p <- ggplot(plot_df, aes(x = group, y = ratio, color = group)) +
    geom_point(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.85) +
    facet_wrap(~ circRNA, scales = "free_y") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = paste0("DCC CircTest ratios: ", comparison),
      x = NULL,
      y = "circular / (circular + linear)"
    )
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_pdf, p, width = 10, height = max(5, ceiling(length(sig_ids) / 3) * 2.1), limitsize = FALSE)
}

circ_counts <- read_dcc_counts(snakemake@input[["circ_counts"]], all_samples)
linear_counts <- read_dcc_counts(snakemake@input[["linear_counts"]], all_samples)
coord_annot <- read_dcc_coordinates(snakemake@input[["circ_coordinates"]])

message("DCC CircRNACount rows: ", nrow(circ_counts))
message("DCC LinearCount rows: ", nrow(linear_counts))

if (nrow(circ_counts) == 0) {
  empty <- empty_results()
  filtered <- data.frame(circRNA = character(0), stringsAsFactors = FALSE)
  for (sample in selected_samples) {
    filtered[[sample]] <- integer(0)
  }
  write_table(empty, snakemake@output[["results"]])
  write_table(filtered, snakemake@output[["filtered_circ"]])
  write_table(filtered, snakemake@output[["filtered_linear"]])
  write_summary(0, 0, 0, 0, "no_dcc_circrnas")
  write_blank_pdf(snakemake@output[["ratio_plots"]], "No DCC circRNAs were available for CircTest.")
  quit(save = "no", status = 0)
}

tables <- make_circtest_tables(circ_counts, linear_counts, selected_samples)
Circ <- tables$Circ
Linear <- tables$Linear
if (tables$dropped > 0) {
  message("Dropped ", tables$dropped, " circRNAs missing matching LinearCount rows.")
}

if (nrow(Circ) == 0) {
  empty <- empty_results()
  write_table(empty, snakemake@output[["results"]])
  write_table(Circ, snakemake@output[["filtered_circ"]])
  write_table(Linear, snakemake@output[["filtered_linear"]])
  write_summary(nrow(circ_counts), 0, 0, 0, "no_matching_linear_counts")
  write_blank_pdf(snakemake@output[["ratio_plots"]], "No DCC circRNAs had matching LinearCount rows.")
  quit(save = "no", status = 0)
}

min_replicates <- min(table(group_labels))
message(
  "Filtering DCC CircTest input with filter.sample=", filter_sample,
  ", filter.count=", filter_count,
  ", percentage=", percentage,
  ", Nreplicates=", min_replicates
)

Circ_filtered <- tryCatch(
  {
    Circ.filter(
      circ = Circ,
      linear = Linear,
      Nreplicates = min_replicates,
      filter.sample = filter_sample,
      filter.count = filter_count,
      percentage = percentage,
      circle_description = 1
    )
  },
  error = function(e) {
    message("Circ.filter failed; using simple count/ratio filter. Error: ", conditionMessage(e))
    simple_filter(Circ, Linear)
  }
)
Linear_filtered <- Linear[rownames(Circ_filtered), , drop = FALSE]

write_table(Circ_filtered, snakemake@output[["filtered_circ"]])
write_table(Linear_filtered, snakemake@output[["filtered_linear"]])
message("CircRNAs after CircTest filtering: ", nrow(Circ_filtered))

if (nrow(Circ_filtered) == 0) {
  empty <- empty_results()
  write_table(empty, snakemake@output[["results"]])
  write_summary(nrow(circ_counts), 0, 0, 0, "no_circrnas_after_filtering")
  write_blank_pdf(snakemake@output[["ratio_plots"]], "No DCC circRNAs passed CircTest filtering.")
  quit(save = "no", status = 0)
}

group_numeric <- as.integer(factor(group_labels, levels = unique(group_labels)))
test <- Circ.test(
  Circ_filtered,
  Linear_filtered,
  group = group_numeric,
  circle_description = 1
)

pvalue <- rep(NA_real_, nrow(Circ_filtered))
padj <- rep(NA_real_, nrow(Circ_filtered))
if (!is.null(test[["p.val"]])) {
  pvalue <- as.numeric(test[["p.val"]])
}
if (!is.null(test[["p.adj"]])) {
  padj <- as.numeric(test[["p.adj"]])
}
if (length(pvalue) != nrow(Circ_filtered)) {
  warning("CircTest p.val length does not match filtered circRNA count; filling pvalue with NA.")
  pvalue <- rep(NA_real_, nrow(Circ_filtered))
}
if (length(padj) != nrow(Circ_filtered)) {
  warning("CircTest p.adj length does not match filtered circRNA count; filling padj with NA.")
  padj <- rep(NA_real_, nrow(Circ_filtered))
}

tested_ids <- Circ_filtered$circRNA
result_df <- data.frame(
  circRNA = tested_ids,
  pvalue = pvalue,
  padj = padj,
  regulation = ifelse(!is.na(padj) & padj < padj_cutoff, "Significant", "Not_significant"),
  stringsAsFactors = FALSE
)

coord_by_id <- coord_annot[match(result_df$circRNA, coord_annot$circRNA), , drop = FALSE]
coord_by_id$circRNA <- result_df$circRNA
result_df <- cbind(
  data.frame(
    chrom = sub(":.*$", "", result_df$circRNA),
    start = suppressWarnings(as.integer(sub("\\|.*$", "", sub("^[^:]+:", "", result_df$circRNA)))),
    end = suppressWarnings(as.integer(sub("^.*\\|", "", result_df$circRNA))),
    gene_id = coord_by_id$gene_id,
    strand = coord_by_id$strand,
    junction_type = coord_by_id$junction_type,
    stringsAsFactors = FALSE
  ),
  result_df
)
result_df <- result_df[, c(
  "circRNA", "chrom", "start", "end", "gene_id", "strand",
  "junction_type", "pvalue", "padj", "regulation"
)]
result_df <- result_df[order(result_df$padj, result_df$pvalue, na.last = TRUE), , drop = FALSE]

write_table(result_df, snakemake@output[["results"]])
significant_n <- sum(!is.na(result_df$padj) & result_df$padj < padj_cutoff)
write_summary(
  nrow(circ_counts),
  nrow(Circ_filtered),
  nrow(result_df),
  significant_n,
  "ok"
)
plot_ratios(Circ_filtered, Linear_filtered, result_df, snakemake@output[["ratio_plots"]])

message("DCC CircTest complete: ", nrow(result_df), " tested, ", significant_n, " significant.")
