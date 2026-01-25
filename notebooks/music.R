library(MuSiC)
library(Biobase)
library(SingleCellExperiment)

# -----------------------------
# Data loading / preprocessing
# -----------------------------
Hao_sce <- readRDS("~/Uni/Rectangle/omnideconv_runtimes/Hao_sce.rds")
genes <- Hao_sce@NAMES
cellnames <- Hao_sce@colData@rownames
counts <- Hao_sce@assays@data@listData[["counts"]]
bulks <- readRDS("~/Uni/Rectangle/omnideconv_runtimes/finotello_tpm.rds")
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

# Ensure peakRAM exists
if (!requireNamespace("peakRAM", quietly = TRUE)) {
  install.packages("peakRAM")
}
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

  # Conversion: Words (or weird unit) -> Bytes -> MiB
  peak_mem_mib <- (raw_peak_val * 8) / (1024 * 1024)

  runtimes[i] <- duration_secs
  memories[i] <- peak_mem_mib

  cat(sprintf("  Run %d: %.4fs, %.2f MiB (converted from raw %.0f)\n", i, duration_secs, peak_mem_mib, raw_peak_val))
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
