#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
data_dir <- file.path(repo_root, "data")

hao_default <- file.path(data_dir, "Hao_sce.rds")
bulks_csv_default <- file.path(data_dir, "bulks_finotello_counts.csv")
bulks_rds_default <- file.path(data_dir, "bulks_finotello_counts.rds")

hao_path <- Sys.getenv("HAO_SCE_PATH", unset = hao_default)
bulks_path <- Sys.getenv(
  "BULKS_PATH",
  unset = if (file.exists(bulks_csv_default)) bulks_csv_default else bulks_rds_default
)

cran_repo <- getOption("repos")[["CRAN"]]
if (is.null(cran_repo) || !nzchar(cran_repo) || cran_repo %in% c("@CRAN", "@CRAN@")) {
  cran_repo <- "https://cloud.r-project.org"
}

ensure_cran_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = cran_repo)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' could not be installed.", pkg))
  }
}

ensure_bioc_pkg <- function(pkg) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repo)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required Bioconductor package '%s' could not be installed.", pkg))
  }
}

ensure_bisque_pkg <- function() {
  if (requireNamespace("BisqueRNA", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  # Preferred: Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repo)
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    try(BiocManager::install("BisqueRNA", ask = FALSE, update = FALSE), silent = TRUE)
  }
  if (requireNamespace("BisqueRNA", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  # Fallback: GitHub source install
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = cran_repo)
  }
  if (requireNamespace("remotes", quietly = TRUE)) {
    try(remotes::install_github("cozygene/bisque", upgrade = "never"), silent = TRUE)
  }
  if (!requireNamespace("BisqueRNA", quietly = TRUE)) {
    stop(
      paste(
        "Required package 'BisqueRNA' could not be installed.",
        "Please install it manually via BiocManager::install('BisqueRNA')",
        "or remotes::install_github('cozygene/bisque')."
      )
    )
  }
}

ensure_bioc_pkg("Biobase")
ensure_bisque_pkg()
ensure_cran_pkg("peakRAM")
ensure_cran_pkg("Matrix")

library(Biobase)
library(BisqueRNA)

# -----------------------------
# Data loading / preprocessing
# -----------------------------
if (!file.exists(hao_path)) {
  stop(sprintf("Missing Hao_sce input: %s", hao_path))
}
if (!file.exists(bulks_path)) {
  stop(sprintf("Missing bulk input: %s", bulks_path))
}

Hao_sce <- readRDS(hao_path)
genes <- Hao_sce@NAMES
cellnames <- Hao_sce@colData@rownames
counts <- Hao_sce@assays@data@listData[["counts"]]
if (grepl("\\.rds$", bulks_path, ignore.case = TRUE)) {
  bulks <- readRDS(bulks_path)
} else if (grepl("\\.csv$", bulks_path, ignore.case = TRUE)) {
  bulks <- utils::read.csv(bulks_path, row.names = 1, check.names = FALSE)
} else {
  stop(sprintf("Unsupported bulk input extension for: %s (use .csv or .rds)", bulks_path))
}
annotations <- as.character(Hao_sce@colData@listData[["cell_type"]])

sc_counts <- counts
bulk_mat <- bulks

if (is.null(rownames(sc_counts))) {
  rownames(sc_counts) <- genes
}
if (is.null(colnames(sc_counts))) {
  colnames(sc_counts) <- cellnames
}

common_genes <- intersect(rownames(bulk_mat), rownames(sc_counts))
bulk_mat <- bulk_mat[common_genes, , drop = FALSE]
sc_counts <- sc_counts[common_genes, , drop = FALSE]

# Downsample cells before dense conversion to keep memory bounded.
max_cells <- 35000L
if (!is.na(max_cells) && max_cells > 0 && ncol(sc_counts) > max_cells) {
  set.seed(42)
  ann <- as.character(annotations)
  idx_by_type <- split(seq_along(ann), ann)
  n_total <- length(ann)

  chosen <- integer(0)
  for (ct in names(idx_by_type)) {
    idx <- idx_by_type[[ct]]
    frac <- as.numeric(length(idx)) / as.numeric(n_total)
    n_pick <- max(1L, floor(as.numeric(max_cells) * frac))
    n_pick <- min(n_pick, length(idx))
    chosen <- c(chosen, sample(idx, as.integer(n_pick)))
  }

  if (length(chosen) > max_cells) {
    chosen <- sample(chosen, max_cells)
  } else if (length(chosen) < max_cells) {
    remaining <- setdiff(seq_along(ann), chosen)
    if (length(remaining) > 0) {
      chosen <- c(chosen, sample(remaining, min(max_cells - length(chosen), length(remaining))))
    }
  }

  chosen <- sort(unique(chosen))
  sc_counts <- sc_counts[, chosen, drop = FALSE]
  annotations <- annotations[chosen]
  cat(sprintf("Downsampled single-cell data to %d cells (target: %d).\n", ncol(sc_counts), max_cells))
}

# Keep all genes (hardcoded): no pre-filtering before dense conversion.
max_genes <- 0L
if (!is.na(max_genes) && max_genes > 0 && nrow(sc_counts) > max_genes) {
  mu <- Matrix::rowMeans(sc_counts)
  mu2 <- Matrix::rowMeans(sc_counts^2)
  vars <- pmax(mu2 - (mu^2), 0)
  keep <- order(vars, decreasing = TRUE)[seq_len(max_genes)]
  sc_counts <- sc_counts[keep, , drop = FALSE]
  bulk_mat <- bulk_mat[rownames(sc_counts), , drop = FALSE]
  cat(sprintf("Filtered to top %d variable genes before dense conversion.\n", nrow(sc_counts)))
}

# Bisque expects both cellType and subjectName in pData(sc.eset).
# Hao input has one donor here, so use pseudo-subjects to avoid degenerate fits.
set.seed(42)
pseudo_subject <- sample(c("S1_a", "S1_b"), size = ncol(sc_counts), replace = TRUE)

sc_pheno <- data.frame(
  row.names = colnames(sc_counts),
  cellType = annotations,
  SubjectName = pseudo_subject,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sc_meta <- data.frame(
  labelDescription = c("Cell type label", "Subject / donor id"),
  row.names = c("cellType", "SubjectName")
)
sc_pdata <- methods::new("AnnotatedDataFrame", data = sc_pheno, varMetadata = sc_meta)

# Keep matrices sparse until after gene intersection, then convert because
# ExpressionSet/Bisque workflows expect dense assay matrices.
sc_counts <- as.matrix(sc_counts)
bulk_mat <- as.matrix(bulk_mat)

sc_eset <- Biobase::ExpressionSet(assayData = sc_counts, phenoData = sc_pdata)

bulk_pdata <- Biobase::AnnotatedDataFrame(data.frame(row.names = colnames(bulk_mat)))
bulk_eset <- Biobase::ExpressionSet(assayData = bulk_mat, phenoData = bulk_pdata)

rm(Hao_sce, counts)

# -----------------------------
# Benchmarking Bisque (peakRAM)
# -----------------------------
library(peakRAM)

n_runs <- 3
runtimes <- numeric(n_runs)
memories <- numeric(n_runs)

cat(sprintf("Starting Bisque Benchmark (%d runs)...\n", n_runs))

for (i in 1:n_runs) {
  cat(sprintf("\nRun %d/%d...\n", i, n_runs))
  gc()

  start_time <- Sys.time()
  pr <- peakRAM::peakRAM({
    fit <- BisqueRNA::ReferenceBasedDecomposition(
      bulk.eset = bulk_eset,
      sc.eset = sc_eset,
      markers = NULL,
      use.overlap = FALSE,
      old.cpm = FALSE,
      verbose = FALSE
    )
  })
  end_time <- Sys.time()
  print(pr)

  duration_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))

  raw_peak_val <- NA_real_
  if ("Peak_RAM_Used_MiB" %in% colnames(pr)) {
    raw_peak_val <- pr$Peak_RAM_Used_MiB[1]
  } else {
    raw_peak_val <- max(as.numeric(pr[1, sapply(pr, is.numeric)]), na.rm = TRUE)
  }

  # Some systems report peakRAM in words despite MiB column naming.
  peak_mem_mib <- if (!is.na(raw_peak_val) && raw_peak_val > 1e6) {
    (raw_peak_val * 8) / (1024 * 1024)
  } else {
    raw_peak_val
  }

  runtimes[i] <- duration_secs
  memories[i] <- peak_mem_mib

  cat(sprintf("  Run %d: %.4fs, %.2f MiB (peakRAM)\n", i, duration_secs, peak_mem_mib))
}

cat("\n--- Bisque Benchmark Results Summary ---\n")
res_df <- data.frame(Run = 1:n_runs, Time_s = runtimes, Peak_Mem_MiB = memories)
print(res_df)

cat(sprintf("\nAverage Runtime: %.4f seconds (sd: %.4f)\n", mean(runtimes), sd(runtimes)))
cat(sprintf("Average Peak Memory: %.2f MiB (sd: %.2f)\n", mean(memories), sd(memories)))

csv_filename <- "bisque_benchmark_results.csv"
write.csv(res_df, csv_filename, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", csv_filename))
