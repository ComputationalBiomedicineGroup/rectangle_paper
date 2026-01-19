library(SCDC)
library(Biobase)

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
csv_filename <- "scdc_benchmark_results.csv"
write.csv(res_df, csv_filename, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", csv_filename))
