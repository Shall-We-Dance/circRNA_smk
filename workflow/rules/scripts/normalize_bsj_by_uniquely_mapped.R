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
bsj_matrix_file <- input_path("bsj_matrix")
output_file     <- output_path("normalized")
star_dir        <- as.character(param("star_dir"))
per_million     <- as.numeric(param("per_million", 1e6))

if (is.null(star_dir) || !nzchar(star_dir)) {
  stop("snakemake.params.star_dir is required.")
}
if (!is.finite(per_million) || per_million <= 0) {
  stop("snakemake.params.per_million must be a positive number; got: ", per_million)
}
if (!dir.exists(star_dir)) {
  stop("STAR directory not found: ", star_dir)
}

message("BSJ matrix:          ", bsj_matrix_file)
message("STAR root directory: ", star_dir)
message("Output (CPM TSV):    ", output_file)
message("per_million scale:   ", per_million)

# -----------------------------
# Helpers
# -----------------------------
safe_read_tsv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read.delim(
    path,
    header           = TRUE,
    sep              = "\t",
    check.names      = FALSE,
    quote            = "",
    comment.char     = "",
    stringsAsFactors = FALSE
  )
}

write_tsv_base <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
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

# -----------------------------
# Main
# -----------------------------
bsj <- safe_read_tsv(bsj_matrix_file)

if (!"circRNA_ID" %in% colnames(bsj)) {
  stop("BSJ matrix must contain a 'circRNA_ID' column: ", bsj_matrix_file)
}
if (colnames(bsj)[1] != "circRNA_ID") {
  bsj <- bsj[, c("circRNA_ID", setdiff(colnames(bsj), "circRNA_ID")), drop = FALSE]
}

sample_cols <- setdiff(colnames(bsj), "circRNA_ID")
if (length(sample_cols) == 0) {
  stop("BSJ matrix has no sample columns alongside 'circRNA_ID': ", bsj_matrix_file)
}

message("Reading STAR Log.final for ", length(sample_cols), " sample(s)...")
unique_map_vec <- vapply(
  sample_cols,
  function(s) extract_uniquely_mapped_reads(find_star_log_file(s, star_dir)),
  FUN.VALUE = numeric(1)
)
names(unique_map_vec) <- sample_cols

if (any(is.na(unique_map_vec) | unique_map_vec <= 0)) {
  bad <- names(unique_map_vec)[is.na(unique_map_vec) | unique_map_vec <= 0]
  stop("Non-positive or missing uniquely mapped read counts for: ",
       paste(bad, collapse = ", "))
}

# Report per-sample unique-map counts for traceability in the log.
for (s in sample_cols) {
  message(sprintf("  %s: uniquely_mapped_reads = %s", s,
                  format(unique_map_vec[[s]], big.mark = ",", scientific = FALSE)))
}

# Coerce sample columns to numeric, normalize, write.
bsj_numeric <- bsj[, sample_cols, drop = FALSE]
for (j in seq_along(bsj_numeric)) bsj_numeric[[j]] <- as.numeric(bsj_numeric[[j]])
bsj_mat <- as.matrix(bsj_numeric)
rownames(bsj_mat) <- bsj$circRNA_ID

bsj_cpm <- sweep(bsj_mat, 2, unique_map_vec[colnames(bsj_mat)], "/") * per_million

out_df <- data.frame(
  circRNA_ID = rownames(bsj_cpm),
  bsj_cpm,
  check.names = FALSE
)

out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
write_tsv_base(out_df, output_file)

message("Wrote normalized BSJ matrix: ", output_file)
message("  ", nrow(out_df), " circRNAs x ", length(sample_cols), " samples")
