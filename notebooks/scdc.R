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
  ensure_cran_pkg("BiocManager")
  options(repos = BiocManager::repositories())
  if (!requireNamespace(pkg, quietly = TRUE)) {
    try(BiocManager::install(pkg, ask = FALSE, update = FALSE), silent = TRUE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required Bioconductor package '%s' could not be installed.", pkg))
  }
}

ensure_bioc_stack <- function(pkgs) {
  for (pkg in pkgs) {
    ensure_bioc_pkg(pkg)
  }
}

ensure_scdc_pkg <- function() {
  if (requireNamespace("SCDC", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  # Preferred: Bioconductor
  ensure_cran_pkg("BiocManager")
  options(repos = BiocManager::repositories())
  try(BiocManager::install("SCDC", ask = FALSE, update = FALSE), silent = TRUE)
  if (requireNamespace("SCDC", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  # Fallback: GitHub
  ensure_cran_pkg("remotes")
  try(remotes::install_github("meichendong/SCDC", upgrade = "never"), silent = TRUE)
  if (requireNamespace("SCDC", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  # Some SCDC installs require xbioc and deeper Bioconductor dependencies.
  ensure_bioc_stack(c(
    "Biobase",
    "BiocGenerics",
    "BiocVersion",
    "Biostrings",
    "AnnotationDbi",
    "KEGGREST"
  ))
  try(remotes::install_github("renozao/xbioc", upgrade = "never"), silent = TRUE)
  try(remotes::install_github("meichendong/SCDC", upgrade = "never"), silent = TRUE)

  if (!requireNamespace("SCDC", quietly = TRUE)) {
    stop(
      paste(
        "Required package 'SCDC' could not be installed.",
        "Please install it manually via remotes::install_github('meichendong/SCDC')."
      )
    )
  }
}

ensure_scdc_pkg()
ensure_bioc_pkg("Biobase")
library(SCDC)
library(Biobase)

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

max_cells <- 20000L
if (!is.na(max_cells) && max_cells > 0 && ncol(counts) > max_cells) {
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
  counts <- counts[, chosen, drop = FALSE]
  annotations <- annotations[chosen]
  cellnames <- colnames(counts)
  cat(sprintf("Downsampled single-cell data to %d cells (target: %d).\n", ncol(counts), max_cells))
}

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

eset <- (Biobase::ExpressionSet(
  assayData = as.matrix(counts),
  phenoData = sc_pdata
))

bulk_mat <- as.matrix(bulks)
common_genes <- intersect(rownames(bulk_mat), rownames(exprs(eset)))
bulk_mat2 <- bulk_mat[common_genes, , drop = FALSE]
eset <- eset[common_genes, ]
bulk_pdata <- AnnotatedDataFrame(data.frame(row.names = colnames(bulk_mat2)))
bulk_eset <- ExpressionSet(assayData = bulk_mat2, phenoData = bulk_pdata)
ct.sub <- sort(unique(pData(eset)$cellType))

rm(counts)
rm(Hao_sce)

# -----------------------------
# Benchmarking SCDC (peakRAM)
# -----------------------------

# Ensure peakRAM exists
if (!requireNamespace("peakRAM", quietly = TRUE)) {
  install.packages("peakRAM")
}
library(peakRAM)

n_runs <- 3
runtimes <- numeric(n_runs)
memories <- numeric(n_runs)

cat(sprintf("Starting Benchmark (%d runs)...\n", n_runs))

# Ensure sample_id exists for SCDC_prop_ONE
pData(eset)$sample_id <- "sample1"

for (i in 1:n_runs) {
  cat(sprintf("\nRun %d/%d...\n", i, n_runs))

  # 1) Warm-up / clean-up
  gc()

  # 2) Measure time + peak RAM around the target call
  start_time <- Sys.time()

  pr <- peakRAM::peakRAM({
    fit <- SCDC_prop_ONE(
      bulk.eset = bulk_eset,
      sc.eset = eset,
      ct.varname = "cellType",
      sample = "sample_id",
      ct.sub = ct.sub,
      weight.basis = TRUE
    )
  })
  print(pr)

  end_time <- Sys.time()

  # 3) Record results
  duration_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # peakRAM returns a data.frame with columns like:
  # Elapsed_Time_sec, Total_RAM_Used_MiB, Peak_RAM_Used_MiB
  peak_mem_mib <- NA_real_
  if ("Peak_RAM_Used_MiB" %in% colnames(pr)) {
    peak_mem_mib <- pr$Peak_RAM_Used_MiB[1]
  } else {
    # Fallback: take the max numeric value if column name differs
    peak_mem_mib <- max(as.numeric(pr[1, sapply(pr, is.numeric)]), na.rm = TRUE)
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
csv_filename <- "scdc_benchmark_results_20k.csv"
write.csv(res_df, csv_filename, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", csv_filename))
