suppressWarnings(options(stringsAsFactors = FALSE))

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      stop("Missing value for --", key)
    }
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

opts <- parse_args(args)

required <- c("method", "info", "bsj", "out")
missing <- required[!required %in% names(opts)]
if (length(missing) > 0) {
  stop("Missing required argument(s): ", paste(missing, collapse = ", "))
}

method <- opts$method
if (!method %in% c("de_ratio", "de_relative")) {
  stop("--method must be 'de_ratio' or 'de_relative'.")
}
if (method == "de_ratio" && (!"fsj" %in% names(opts))) {
  stop("--fsj is required for de_ratio.")
}
if (method == "de_relative" && (!"circ-gene" %in% names(opts))) {
  stop("--circ-gene is required for de_relative.")
}

read_table <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
}

to_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

read_count_matrix <- function(path) {
  df <- read_table(path)
  if (ncol(df) < 2) {
    stop("Count matrix must contain an ID column and at least one sample column: ", path)
  }
  ids <- as.character(df[[1]])
  mat <- as.matrix(data.frame(lapply(df[-1], to_numeric), check.names = FALSE))
  colnames(mat) <- colnames(df)[-1]
  rownames(mat) <- ids
  mat[is.na(mat)] <- 0
  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }
  mat
}

first_existing_column <- function(df, candidates) {
  for (candidate in candidates) {
    if (candidate %in% colnames(df)) return(candidate)
  }
  NULL
}

looks_like_circ <- function(x) {
  grepl("^[^:]+:[0-9]+\\|[0-9]+", as.character(x))
}

read_official_stats <- function(path, event_ids) {
  df <- tryCatch(read_table(path), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)

  p_col <- first_existing_column(
    df,
    c("PValue", "P.Value", "pvalue", "p_value", "Pvalue", "PVAL", "P.Val", "pval")
  )
  if (is.null(p_col)) return(NULL)

  id_col <- first_existing_column(df, c("circRNA_ID", "circRNA", "CircRNA", "ID"))
  if (is.null(id_col) || !any(looks_like_circ(df[[id_col]]))) {
    circ_cols <- colnames(df)[vapply(df, function(col) any(looks_like_circ(col)), logical(1))]
    if (length(circ_cols) == 0) return(NULL)
    id_col <- circ_cols[[1]]
  }

  stats <- data.frame(
    circRNA_ID = as.character(df[[id_col]]),
    PValue = to_numeric(df[[p_col]]),
    stringsAsFactors = FALSE
  )
  stats <- stats[!is.na(stats$PValue) & nzchar(stats$circRNA_ID), , drop = FALSE]
  if (nrow(stats) == 0) return(NULL)

  stats <- stats[!duplicated(stats$circRNA_ID), , drop = FALSE]
  stats$match_count <- sum(stats$circRNA_ID %in% event_ids)
  stats
}

find_official_stats <- function(out_path, event_ids) {
  candidates <- Sys.glob(paste0(out_path, "_*"))
  candidates <- candidates[!grepl("P-V\\.txt$", candidates)]
  candidates <- candidates[file.exists(candidates) & file.info(candidates)$size > 0]
  if (length(candidates) == 0) return(NULL)

  stats_list <- lapply(candidates, read_official_stats, event_ids = event_ids)
  stats_list <- Filter(Negate(is.null), stats_list)
  if (length(stats_list) == 0) return(NULL)

  match_counts <- vapply(stats_list, function(x) unique(x$match_count), numeric(1))
  best <- stats_list[[which.max(match_counts)]]
  best$match_count <- NULL
  if (sum(best$circRNA_ID %in% event_ids) == 0) return(NULL)
  best
}

safe_name <- function(x) {
  gsub("[^[:alnum:]_.-]+", "_", x)
}

native_pvalue <- function(success, failure, group, group_levels) {
  success <- as.numeric(success)
  failure <- as.numeric(failure)
  total <- success + failure
  keep <- is.finite(total) & total > 0 & is.finite(success) & is.finite(failure)
  if (sum(keep) < 2) return(NA_real_)
  if (length(unique(group[keep])) < 2) return(NA_real_)

  dat <- data.frame(
    success = round(success[keep]),
    failure = round(failure[keep]),
    group = factor(group[keep], levels = group_levels),
    stringsAsFactors = FALSE
  )
  if (length(unique(dat$group)) < 2) return(NA_real_)
  if (all(dat$success == 0) || all(dat$failure == 0)) return(1)

  p <- tryCatch(
    {
      fit0 <- suppressWarnings(glm(cbind(success, failure) ~ 1, family = quasibinomial(), data = dat))
      fit1 <- suppressWarnings(glm(cbind(success, failure) ~ group, family = quasibinomial(), data = dat))
      an <- suppressWarnings(anova(fit0, fit1, test = "F"))
      as.numeric(an[2, "Pr(>F)"])
    },
    error = function(e) NA_real_
  )
  if (is.finite(p)) return(max(min(p, 1), 0))

  mat <- rbind(
    c(sum(dat$success[dat$group == group_levels[[1]]]), sum(dat$failure[dat$group == group_levels[[1]]])),
    c(sum(dat$success[dat$group == group_levels[[2]]]), sum(dat$failure[dat$group == group_levels[[2]]]))
  )
  tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
}

build_metrics <- function(ids, values, successes, case_label, control_label, group) {
  case_samples <- group == case_label
  control_samples <- group == control_label
  mean_case <- rowMeans(values[, case_samples, drop = FALSE], na.rm = TRUE)
  mean_control <- rowMeans(values[, control_samples, drop = FALSE], na.rm = TRUE)
  positive_rate <- rowMeans(successes > 0, na.rm = TRUE)
  pseudocount <- 1e-6

  result <- data.frame(
    circRNA_ID = ids,
    stringsAsFactors = FALSE
  )
  result[[paste0("Mean_", safe_name(case_label))]] <- mean_case
  result[[paste0("Mean_", safe_name(control_label))]] <- mean_control
  result$FC <- (mean_case + pseudocount) / (mean_control + pseudocount)
  result$log2FC <- log2(result$FC)
  result$DeltaPSI <- mean_case - mean_control
  result$PositiveRate <- positive_rate
  result
}

info <- read_table(opts$info)
if (!all(c("Sample", "Class") %in% colnames(info))) {
  stop("CIRI3 DE infor.tsv must contain Sample and Class columns: ", opts$info)
}
info$Sample <- as.character(info$Sample)
info$Class <- as.character(info$Class)

if (all(c("Case", "Control") %in% unique(info$Class))) {
  case_label <- "Case"
  control_label <- "Control"
} else {
  class_levels <- unique(info$Class)
  if (length(class_levels) != 2) {
    stop("CIRI3 ratio/relative DEG requires exactly two classes in ", opts$info)
  }
  case_label <- class_levels[[1]]
  control_label <- class_levels[[2]]
}
group_levels <- c(control_label, case_label)
selected_samples <- info$Sample
group <- setNames(info$Class, info$Sample)

bsj <- read_count_matrix(opts$bsj)
missing_bsj <- setdiff(selected_samples, colnames(bsj))
if (length(missing_bsj) > 0) {
  stop("Selected samples missing from BSJ matrix: ", paste(missing_bsj, collapse = ", "))
}
bsj <- bsj[, selected_samples, drop = FALSE]
group <- group[selected_samples]

if (method == "de_ratio") {
  fsj <- read_count_matrix(opts$fsj)
  missing_fsj <- setdiff(selected_samples, colnames(fsj))
  if (length(missing_fsj) > 0) {
    stop("Selected samples missing from FSJ matrix: ", paste(missing_fsj, collapse = ", "))
  }
  common_rows <- intersect(rownames(bsj), rownames(fsj))
  if (length(common_rows) == 0) {
    stop("BSJ and FSJ matrices have no overlapping circRNAs.")
  }
  bsj <- bsj[common_rows, , drop = FALSE]
  fsj <- fsj[common_rows, selected_samples, drop = FALSE]
  successes <- 2 * bsj
  failures <- fsj
  totals <- successes + failures
  values <- successes / totals
  values[!is.finite(values)] <- 0
  result <- build_metrics(
    rownames(values),
    values,
    successes,
    case_label,
    control_label,
    as.character(group)
  )
  successes_for_test <- successes
  failures_for_test <- failures
} else {
  circ_gene <- read_table(opts[["circ-gene"]])
  if (!all(c("circRNA_ID", "gene_id") %in% colnames(circ_gene))) {
    stop("circ_Gene.txt must contain circRNA_ID and gene_id columns: ", opts[["circ-gene"]])
  }
  circ_gene$circRNA_ID <- as.character(circ_gene$circRNA_ID)
  circ_gene$gene_id <- as.character(circ_gene$gene_id)
  circ_gene <- circ_gene[
    nzchar(circ_gene$circRNA_ID) &
      nzchar(circ_gene$gene_id) &
      circ_gene$gene_id != "NA" &
      circ_gene$circRNA_ID %in% rownames(bsj),
    c("circRNA_ID", "gene_id"),
    drop = FALSE
  ]
  circ_gene <- circ_gene[!duplicated(circ_gene), , drop = FALSE]

  if (nrow(circ_gene) == 0) {
    result <- data.frame(
      circRNA_ID = character(0),
      Mean_Case = numeric(0),
      Mean_Control = numeric(0),
      FC = numeric(0),
      log2FC = numeric(0),
      DeltaPSI = numeric(0),
      PositiveRate = numeric(0),
      stringsAsFactors = FALSE
    )
    successes_for_test <- matrix(nrow = 0, ncol = length(selected_samples))
    failures_for_test <- matrix(nrow = 0, ncol = length(selected_samples))
    colnames(successes_for_test) <- selected_samples
    colnames(failures_for_test) <- selected_samples
  } else {
    pieces <- list()
    success_pieces <- list()
    failure_pieces <- list()
    gene_split <- split(circ_gene$circRNA_ID, circ_gene$gene_id)
    for (gene_id in names(gene_split)) {
      circs <- unique(gene_split[[gene_id]])
      gene_counts <- colSums(bsj[circs, , drop = FALSE])
      successes <- bsj[circs, , drop = FALSE]
      failures <- sweep(successes, 2, gene_counts, FUN = function(x, total) total - x)
      totals <- sweep(successes, 2, gene_counts, FUN = function(x, total) total)
      values <- successes / totals
      values[!is.finite(values)] <- 0
      rownames(values) <- paste0(rownames(successes), ";", gene_id)
      rownames(successes) <- rownames(values)
      rownames(failures) <- rownames(values)

      pieces[[gene_id]] <- build_metrics(
        rownames(values),
        values,
        successes,
        case_label,
        control_label,
        as.character(group)
      )
      success_pieces[[gene_id]] <- successes
      failure_pieces[[gene_id]] <- failures
    }
    result <- do.call(rbind, pieces)
    rownames(result) <- NULL
    successes_for_test <- do.call(rbind, success_pieces)
    failures_for_test <- do.call(rbind, failure_pieces)
  }
}

result$PValue <- NA_real_
official <- find_official_stats(opts$out, result$circRNA_ID)
if (!is.null(official)) {
  idx <- match(result$circRNA_ID, official$circRNA_ID)
  official_p <- official$PValue[idx]
  used <- sum(!is.na(official_p))
  message("Using CIRI3/rMATS PValue for ", used, " of ", nrow(result), " events.")
  result$PValue[!is.na(official_p)] <- official_p[!is.na(official_p)]
} else {
  message("No valid CIRI3/rMATS PValue output found; using native proportion GLM fallback.")
}

missing_p <- !is.finite(result$PValue)
if (any(missing_p)) {
  message("Computing native proportion GLM PValue for ", sum(missing_p), " event(s).")
  ids_to_test <- result$circRNA_ID[missing_p]
  result$PValue[missing_p] <- vapply(ids_to_test, function(id) {
    native_pvalue(successes_for_test[id, ], failures_for_test[id, ], as.character(group), group_levels)
  }, numeric(1))
}

result$FDR <- NA_real_
valid_p <- is.finite(result$PValue)
if (any(valid_p)) {
  result$FDR[valid_p] <- p.adjust(result$PValue[valid_p], method = "fdr")
}

if (nrow(result) > 0) {
  result <- result[order(is.na(result$PValue), result$PValue, result$circRNA_ID), , drop = FALSE]
}

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
write.table(result, file = opts$out, quote = FALSE, sep = "\t", row.names = FALSE, na = "")
message("Wrote ", method, " DEG result with ", nrow(result), " events: ", opts$out)
