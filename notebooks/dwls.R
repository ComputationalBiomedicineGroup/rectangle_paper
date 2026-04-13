# DWLS is not distributed as a standard R package on GitHub.
# The repository provides an R script with functions (Deconvolution_functions.R).
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()

dwls_repo_dir <- file.path(script_dir, "DWLS_repo")
dwls_script <- file.path(dwls_repo_dir, "Deconvolution_functions.R")

if (!file.exists(dwls_script)) {
  if (!nzchar(Sys.which("git"))) {
    stop("DWLS script not found and git is unavailable. Please install git.")
  }
  if (dir.exists(dwls_repo_dir)) {
    unlink(dwls_repo_dir, recursive = TRUE)
  }
  status <- system2("git", c("clone", "--depth", "1", "https://github.com/dtsoucas/DWLS.git", dwls_repo_dir))
  if (!identical(status, 0L) || !file.exists(dwls_script)) {
    stop("Failed to fetch DWLS sources. Check internet access and repository availability.")
  }
}

dwls_env <- new.env(parent = globalenv())
source(dwls_script, local = dwls_env)

if (!exists("solveDampenedWLS", envir = dwls_env, mode = "function")) {
  stop("solveDampenedWLS() not found in Deconvolution_functions.R.")
}
dwls_fn <- get("solveDampenedWLS", envir = dwls_env, mode = "function")
trim_fn <- get("trimData", envir = dwls_env, mode = "function")
build_sig_fn <- NULL
if (exists("buildSignatureMatrixMAST", envir = dwls_env, mode = "function")) {
  build_sig_fn <- get("buildSignatureMatrixMAST", envir = dwls_env, mode = "function")
} else if (exists("buildSignatureMatrixUsingSeurat", envir = dwls_env, mode = "function")) {
  build_sig_fn <- get("buildSignatureMatrixUsingSeurat", envir = dwls_env, mode = "function")
} else {
  stop("No DWLS signature-building function found (buildSignatureMatrixMAST/buildSignatureMatrixUsingSeurat).")
}

# -----------------------------
# Data loading / preprocessing
# -----------------------------
Hao_sce <- readRDS("~/Uni/Rectangle/omnideconv_runtimes/Hao_sce.rds")
genes <- Hao_sce@NAMES
cellnames <- Hao_sce@colData@rownames
counts <- Hao_sce@assays@data@listData[["counts"]]
bulks <- readRDS("~/Uni/Rectangle/omnideconv_runtimes/finotello_tpm.rds")
annotations <- Hao_sce@colData@listData[["cell_type"]]
# DWLS internals construct symbols from cell-type names and use eval(parse()).
# Make labels syntactically safe (e.g., "B cells" -> "B.cells") to avoid parse errors.
annotations <- make.names(as.character(annotations))

# Keep sparse counts as long as possible; convert to dense only after subsetting.
sc_counts <- counts
bulk_mat <- as.matrix(bulks)

if (is.null(rownames(sc_counts))) {
  rownames(sc_counts) <- genes
}
if (is.null(colnames(sc_counts))) {
  colnames(sc_counts) <- cellnames
}

# Downsample to control memory usage in signature-building/deconvolution.
max_cells <- 20000
if (ncol(sc_counts) > max_cells) {
  set.seed(42)
  ann <- as.character(annotations)
  idx_by_type <- split(seq_along(ann), ann)
  n_total <- length(ann)

  chosen <- integer(0)
  for (ct in names(idx_by_type)) {
    idx <- idx_by_type[[ct]]
    n_pick <- max(1L, floor(max_cells * length(idx) / n_total))
    n_pick <- min(n_pick, length(idx))
    chosen <- c(chosen, sample(idx, n_pick))
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
  cat(sprintf("Downsampled Hao single-cell data to %d cells (target: %d).\n", ncol(sc_counts), max_cells))
}

common_genes <- intersect(rownames(bulk_mat), rownames(sc_counts))
bulk_mat <- bulk_mat[common_genes, , drop = FALSE]
sc_counts <- sc_counts[common_genes, , drop = FALSE]

# Keep all genes. Convert after cell subsetting.
sc_counts <- as.matrix(sc_counts)

# -----------------------------
# Benchmarking DWLS (peakRAM)
# -----------------------------
if (!requireNamespace("peakRAM", quietly = TRUE)) {
  install.packages("peakRAM")
}
library(peakRAM)

n_runs <- 1
runtimes <- numeric(n_runs)
memories <- numeric(n_runs)

cat(sprintf("Starting DWLS Benchmark (%d runs)...\n", n_runs))

for (i in 1:n_runs) {
  cat(sprintf("\nRun %d/%d...\n", i, n_runs))

  gc()
  start_time <- Sys.time()

  run_tmp <- file.path(script_dir, sprintf("dwls_run_%d_tmp", i))
  if (dir.exists(run_tmp)) {
    unlink(run_tmp, recursive = TRUE)
  }
  dir.create(run_tmp, recursive = TRUE, showWarnings = FALSE)

  pr <- peakRAM::peakRAM({
    signature <- build_sig_fn(
      scdata = sc_counts,
      id = annotations,
      path = run_tmp,
      diff.cutoff = 0.5,
      pval.cutoff = 0.01
    )
    signature <- as.matrix(signature)

    res <- lapply(seq_len(ncol(bulk_mat)), function(j) {
      bulk_vec <- bulk_mat[, j]
      names(bulk_vec) <- rownames(bulk_mat)
      td <- trim_fn(signature, bulk_vec)
      as.numeric(dwls_fn(td$sig, td$bulk))
    })
    names(res) <- colnames(bulk_mat)
  })
  print(pr)

  end_time <- Sys.time()
  duration_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))

  peak_mem_mib <- NA_real_
  raw_peak_val <- NA_real_
  if ("Peak_RAM_Used_MiB" %in% colnames(pr)) {
    raw_peak_val <- pr$Peak_RAM_Used_MiB[1]
  } else {
    raw_peak_val <- max(as.numeric(pr[1, sapply(pr, is.numeric)]), na.rm = TRUE)
  }

  # Some systems report peakRAM in 8-byte words despite MiB column naming.
  # If the value is implausibly large, convert words -> MiB.
  if (!is.na(raw_peak_val) && raw_peak_val > 1e6) {
    peak_mem_mib <- (raw_peak_val * 8) / (1024 * 1024)
  } else {
    peak_mem_mib <- raw_peak_val
  }

  runtimes[i] <- duration_secs
  memories[i] <- peak_mem_mib

  cat(sprintf("  Run %d: %.4fs, %.2f MiB (peakRAM)\n", i, duration_secs, peak_mem_mib))
}

cat("\n--- DWLS Benchmark Results Summary ---\n")
res_df <- data.frame(Run = 1:n_runs, Time_s = runtimes, Peak_Mem_MiB = memories)
print(res_df)

cat(sprintf("\nAverage Runtime: %.4f seconds (sd: %.4f)\n", mean(runtimes), sd(runtimes)))
cat(sprintf("Average Peak Memory: %.2f MiB (sd: %.2f)\n", mean(memories), sd(memories)))

csv_filename <- "dwls_benchmark_results.csv"
write.csv(res_df, csv_filename, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", csv_filename))
