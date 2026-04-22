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
  ensure_cran_pkg("BiocManager")
  options(repos = BiocManager::repositories())
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required Bioconductor package '%s' could not be installed.", pkg))
  }
}

ensure_bayesprism_pkg <- function() {
  if (requireNamespace("BayesPrism", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  # Prefer git clone + local install (avoids GitHub API dependency).
  if (nzchar(Sys.which("git"))) {
    tmp_dir <- tempfile("bayesprism_")
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    status <- try(system2("git", c("clone", "--depth", "1", "https://github.com/Danko-Lab/BayesPrism.git", tmp_dir)), silent = TRUE)
    if (!inherits(status, "try-error") && identical(status, 0L) && file.exists(file.path(tmp_dir, "BayesPrism", "DESCRIPTION"))) {
      try(system2("R", c("CMD", "INSTALL", file.path(tmp_dir, "BayesPrism"))), silent = TRUE)
    }
  }

  # Fallback to remotes in case git install path fails.
  if (!requireNamespace("BayesPrism", quietly = TRUE)) {
    ensure_cran_pkg("remotes")
    try(remotes::install_github("Danko-Lab/BayesPrism", subdir = "BayesPrism", upgrade = "never"), silent = TRUE)
  }

  if (!requireNamespace("BayesPrism", quietly = TRUE)) {
    stop(
      paste(
        "Required package 'BayesPrism' could not be installed.",
        "Please install it manually via git clone + R CMD INSTALL,",
        "or remotes::install_github('Danko-Lab/BayesPrism', subdir = 'BayesPrism')."
      )
    )
  }
}

ensure_cran_pkg("peakRAM")
ensure_cran_pkg("snowfall")
ensure_cran_pkg("gplots")
ensure_bioc_pkg("Biobase")
ensure_bioc_pkg("BiocParallel")
ensure_bioc_pkg("scran")
ensure_cran_pkg("NMF")
ensure_bayesprism_pkg()

library(BayesPrism)
library(peakRAM)

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
annotations <- as.character(Hao_sce@colData@listData[["cell_type"]])

if (grepl("\\.rds$", bulks_path, ignore.case = TRUE)) {
  bulks <- readRDS(bulks_path)
} else if (grepl("\\.csv$", bulks_path, ignore.case = TRUE)) {
  bulks <- utils::read.csv(bulks_path, row.names = 1, check.names = FALSE)
} else {
  stop(sprintf("Unsupported bulk input extension for: %s (use .csv or .rds)", bulks_path))
}

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

# Hardcoded benchmark settings:
# keep all genes, downsample cells to keep local Docker runs feasible.
max_cells <- 50000
if (ncol(sc_counts) > max_cells) {
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

# BayesPrism expects count matrices. Convert after subsetting.
sc_counts <- as.matrix(sc_counts)
bulk_mat <- as.matrix(bulk_mat)
sc_counts <- pmax(round(sc_counts), 0)
bulk_mat <- pmax(round(bulk_mat), 0)
storage.mode(sc_counts) <- "numeric"
storage.mode(bulk_mat) <- "numeric"

# BayesPrism expects:
# - reference: cell-by-gene
# - mixture: sample-by-gene
reference_mat <- t(sc_counts)
mixture_mat <- t(bulk_mat)

# labels must match number of rows in reference (cells).
cell_type_labels <- annotations
cell_state_labels <- annotations

rm(Hao_sce, counts)

# -----------------------------
# Benchmarking BayesPrism
# -----------------------------
n_runs <- 3
runtimes <- numeric(n_runs)
memories <- numeric(n_runs)

cat(sprintf("Starting BayesPrism Benchmark (%d runs)...\n", n_runs))

for (i in seq_len(n_runs)) {
  cat(sprintf("\nRun %d/%d...\n", i, n_runs))
  gc()

  start_time <- Sys.time()
  pr <- peakRAM::peakRAM({
    bp <- BayesPrism::new.prism(
      reference = reference_mat,
      mixture = mixture_mat,
      input.type = "count.matrix",
      cell.type.labels = cell_type_labels,
      cell.state.labels = cell_state_labels,
      key = NULL,
      outlier.cut = 0.01,
      outlier.fraction = 0.1
    )
    bp_res <- BayesPrism::run.prism(prism = bp, n.cores = 8)
    theta <- BayesPrism::get.fraction(
      bp_res,
      which.theta = "final",
      state.or.type = "type"
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

cat("\n--- BayesPrism Benchmark Results Summary ---\n")
res_df <- data.frame(Run = seq_len(n_runs), Time_s = runtimes, Peak_Mem_MiB = memories)
print(res_df)

cat(sprintf("\nAverage Runtime: %.4f seconds (sd: %.4f)\n", mean(runtimes), sd(runtimes)))
cat(sprintf("Average Peak Memory: %.2f MiB (sd: %.2f)\n", mean(memories), sd(memories)))

csv_filename <- "bayesprism_benchmark_results.csv"
write.csv(res_df, csv_filename, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", csv_filename))
