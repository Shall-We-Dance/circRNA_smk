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
gene_counts_file <- input_path("gene_counts")
out_cleaned      <- output_path("cleaned")
out_normalized   <- output_path("normalized")
star_dir         <- as.character(param("star_dir"))
per_million      <- as.numeric(param("per_million", 1e6))

if (is.null(star_dir) || !nzchar(star_dir)) {
  stop("snakemake.params.star_dir is required.")
}
if (!is.finite(per_million) || per_million <= 0) {
  stop("snakemake.params.per_million must be a positive number; got: ", per_million)
}
if (!dir.exists(star_dir)) {
  stop("STAR directory not found: ", star_dir)
}

message("featureCounts input:    ", gene_counts_file)
message("STAR root directory:    ", star_dir)
message("Cleaned counts output:  ", out_cleaned)
message("CPM-normalized output:  ", out_normalized)
message("per_million scale:      ", per_million)

# -----------------------------
# Helpers
# -----------------------------
safe_read_featurecounts <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  # comment.char = "#" skips the featureCounts command-line header line.
  read.delim(
    path, header = TRUE, sep = "\t",
    check.names = FALSE, quote = "", comment.char = "#",
    stringsAsFactors = FALSE
  )
}

write_tsv_base <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
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

# Strip the STAR coordinate-sorted BAM suffix (and any leading path) so column
# names become bare sample names. Matches the convention used by
# prepare_ciri3_de_bsj_inputs.py elsewhere in the pipeline.
clean_featurecounts_sample_name <- function(x) {
  x <- basename(x)
  x <- sub("\\.?Aligned\\.sortedByCoord\\.out\\.bam$", "", x)
  x
}

# -----------------------------
# Main
# -----------------------------
gene_df <- safe_read_featurecounts(gene_counts_file)

required_anno <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
missing_anno <- setdiff(required_anno, colnames(gene_df))
if (length(missing_anno) > 0) {
  stop("featureCounts file missing expected columns: ",
       paste(missing_anno, collapse = ", "))
}

sample_cols_raw <- setdiff(colnames(gene_df), required_anno)
if (length(sample_cols_raw) == 0) {
  stop("featureCounts file has no sample columns beyond the annotation columns: ",
       gene_counts_file)
}

sample_cols_clean <- vapply(sample_cols_raw, clean_featurecounts_sample_name,
                            FUN.VALUE = character(1))
if (anyDuplicated(sample_cols_clean)) {
  dupes <- unique(sample_cols_clean[duplicated(sample_cols_clean)])
  stop("Duplicate cleaned sample names after featureCounts cleanup: ",
       paste(dupes, collapse = ", "))
}
still_bam <- grep("\\.bam$", sample_cols_clean, value = TRUE)
if (length(still_bam) > 0) {
  warning("featureCounts column names still end in .bam after cleanup; ",
          "STAR output naming may differ from the assumed pattern. Examples: ",
          paste(head(still_bam, 3), collapse = ", "))
}
colnames(gene_df)[match(sample_cols_raw, colnames(gene_df))] <- sample_cols_clean

# Now sample_cols_clean is our canonical per-sample column name list.
sample_cols <- sample_cols_clean
message("Detected ", length(sample_cols), " sample column(s) after cleanup.")

# STAR uniquely mapped reads per sample
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

for (s in sample_cols) {
  message(sprintf("  %s: uniquely_mapped_reads = %s", s,
                  format(unique_map_vec[[s]], big.mark = ",", scientific = FALSE)))
}

# Strict numeric coercion on sample columns
for (nm in sample_cols) {
  gene_df[[nm]] <- strict_numeric(gene_df[[nm]],
                                  context = paste0("totalRNA column '", nm, "'"))
}

# Write cleaned counts (annotation cols + clean sample names, raw counts)
cleaned_out <- gene_df[, c(required_anno, sample_cols), drop = FALSE]
out_cleaned_dir <- dirname(out_cleaned)
if (!dir.exists(out_cleaned_dir)) {
  dir.create(out_cleaned_dir, recursive = TRUE, showWarnings = FALSE)
}
write_tsv_base(cleaned_out, out_cleaned)
message("Wrote cleaned totalRNA counts: ", out_cleaned)

# CPM-normalize the sample columns and write
raw_mat <- as.matrix(gene_df[, sample_cols, drop = FALSE])
rownames(raw_mat) <- as.character(gene_df$Geneid)

norm_mat <- sweep(raw_mat, 2, unique_map_vec[colnames(raw_mat)], "/") * per_million
normalized_out <- data.frame(
  gene_df[, required_anno, drop = FALSE],
  norm_mat,
  check.names = FALSE
)

out_norm_dir <- dirname(out_normalized)
if (!dir.exists(out_norm_dir)) {
  dir.create(out_norm_dir, recursive = TRUE, showWarnings = FALSE)
}
write_tsv_base(normalized_out, out_normalized)
message("Wrote normalized totalRNA matrix: ", out_normalized)
message("  ", nrow(normalized_out), " genes x ", length(sample_cols), " samples")
