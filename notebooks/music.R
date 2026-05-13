args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
data_dir <- file.path(repo_root, "data")

hao_default <- file.path(data_dir, "Hao_sce.rds")
bulks_tpm_default <- file.path(data_dir, "bulks_finotello_tpm.csv")
bulks_counts_default <- file.path(data_dir, "bulks_finotello_counts.rds")

hao_path <- Sys.getenv("HAO_SCE_PATH", unset = hao_default)
bulks_path <- Sys.getenv(
  "BULKS_PATH",
  unset = if (file.exists(bulks_tpm_default)) bulks_tpm_default else bulks_counts_default
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
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    options(repos = BiocManager::repositories())
  }
  if (!requireNamespace(pkg, quietly = TRUE) && requireNamespace("BiocManager", quietly = TRUE)) {
    try(BiocManager::install(pkg, ask = FALSE, update = FALSE), silent = TRUE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required Bioconductor package '%s' could not be installed.", pkg))
  }
}

ensure_bioc_stack <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0L) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repo)
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    options(repos = BiocManager::repositories())
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    try(BiocManager::install(missing, ask = FALSE, update = FALSE), silent = TRUE)
  }
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0L) {
    stop(
      paste(
        sprintf("Failed to install Bioconductor stack: %s.", paste(still_missing, collapse = ", ")),
        "In conda environments, install binaries first:",
        "conda install -c conda-forge bioconductor-iranges bioconductor-s4vectors",
        "bioconductor-delayedarray bioconductor-genomicranges",
        "bioconductor-genomeinfodbdata bioconductor-genomeinfodb bioconductor-summarizedexperiment",
        "bioconductor-singlecellexperiment"
      )
    )
  }
}

ensure_toast_pkg <- function() {
  if (requireNamespace("TOAST", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repo)
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    options(repos = BiocManager::repositories())
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    try(BiocManager::install("TOAST", ask = FALSE, update = FALSE), silent = TRUE)
  }
  if (requireNamespace("TOAST", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = cran_repo)
  }
  if (requireNamespace("remotes", quietly = TRUE)) {
    try(remotes::install_github("xuranw/TOAST", upgrade = "never"), silent = TRUE)
  }

  if (!requireNamespace("TOAST", quietly = TRUE)) {
    stop(
      paste(
        "Required package 'TOAST' could not be installed.",
        "Please install it manually via BiocManager::install('TOAST')",
        "or remotes::install_github('xuranw/TOAST')."
      )
    )
  }
}

ensure_music_pkg <- function() {
  if (requireNamespace("MuSiC", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  # Install core Bioconductor stack first to avoid cascading source-build failures.
  ensure_bioc_stack(c(
    "IRanges", "S4Vectors", "S4Arrays", "DelayedArray",
    "GenomeInfoDbData", "GenomeInfoDb", "GenomicRanges", "SummarizedExperiment",
    "SingleCellExperiment"
  ))

  # Install key dependencies first to avoid GitHub build failures.
  ensure_bioc_pkg("Biobase")
  ensure_bioc_pkg("SingleCellExperiment")
  ensure_toast_pkg()

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repo)
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    try(BiocManager::install("MuSiC", ask = FALSE, update = FALSE), silent = TRUE)
  }
  if (requireNamespace("MuSiC", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = cran_repo)
  }
  if (requireNamespace("remotes", quietly = TRUE)) {
    try(remotes::install_github("xuranw/MuSiC", upgrade = "never"), silent = TRUE)
  }

  if (!requireNamespace("MuSiC", quietly = TRUE)) {
    stop(
      paste(
        "Required package 'MuSiC' could not be installed.",
        "Please install it manually via BiocManager::install('MuSiC')",
        "or remotes::install_github('xuranw/MuSiC')."
      )
    )
  }
}

ensure_music_pkg()
ensure_cran_pkg("peakRAM")
ensure_bioc_pkg("Biobase")
ensure_bioc_pkg("SingleCellExperiment")
library(MuSiC)
library(Biobase)
library(SingleCellExperiment)

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
annotations <- Hao_sce@colData@listData[["cell_type"]]

sc_pheno <- data.frame(
  check.names = FALSE, check.rows = FALSE,
  stringsAsFactors = FALSE,
  row.names = cellnames,
  cellType = annotations
)
sc_meta <- data.frame(
  labelDescription = c("cellType"),
  row.names = c("cellType")
)
sc_pdata <- methods::new("AnnotatedDataFrame",
  data = sc_pheno,
  varMetadata = sc_meta
)



sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = sc_pheno # <-- use data.frame, not AnnotatedDataFrame
)

bulk_mat <- as.matrix(bulks)


colData(sce)$cellType <- annotations

# MuSiC requires multiple subjects to calculate cross-subject variance.
# Since we only have one subject, we split it into two "pseudo-subjects"
# to satisfy the requirement and avoid the colMeans error.
set.seed(42)
colData(sce)$sampleID <- sample(c("S1_a", "S1_b"), size = ncol(sce), replace = TRUE)

if (is.null(rownames(counts(sce)))) {
  rownames(counts(sce)) <- genes
}


# -----------------------------
# Benchmarking SCDC (peakRAM)
# -----------------------------

library(peakRAM)

n_runs <- 5
runtimes <- numeric(n_runs)
memories <- numeric(n_runs)

cat(sprintf("Starting Benchmark (%d runs)...\n", n_runs))



for (i in 1:n_runs) {
  cat(sprintf("\nRun %d/%d...\n", i, n_runs))

  # 1) Warm-up / clean-up
  gc()

  # 2) Measure time + peak RAM around the target call
  start_time <- Sys.time()

  pr <- peakRAM::peakRAM({
    res <- music_prop(
      bulk.mtx  = bulk_mat,
      sc.sce    = sce,
      clusters  = "cellType",
      samples   = "sampleID",
      markers   = NULL, # uses common genes
      verbose   = TRUE
    )
  })
  print(pr)

  end_time <- Sys.time()

  # 3) Record results
  duration_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # peakRAM output handling:
  # On macOS/some systems, peakRAM reports values in 8-byte words (approx), not MiB.
  # We convert strictly based on the observed behavior:
  # Value * 8 bytes / 1024^2 = MiB

  raw_peak_val <- NA_real_
  if ("Peak_RAM_Used_MiB" %in% colnames(pr)) {
    raw_peak_val <- pr$Peak_RAM_Used_MiB[1]
  } else {
    # Fallback
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

cat("\n--- Benchmark Results Summary ---\n")
res_df <- data.frame(Run = 1:n_runs, Time_s = runtimes, Peak_Mem_MiB = memories)
print(res_df)

cat(sprintf("\nAverage Runtime: %.4f seconds (sd: %.4f)\n", mean(runtimes), sd(runtimes)))
cat(sprintf("Average Peak Memory: %.2f MiB (sd: %.2f)\n", mean(memories), sd(memories)))

# Save results to CSV
csv_filename <- "music_benchmark_results.csv"
write.csv(res_df, csv_filename, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", csv_filename))
