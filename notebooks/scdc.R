library(SCDC)


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
  labelDescription = c(
    "cellType"
  ),
  row.names = c(
    "cellType"
  )
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


# ------------------------------------------------------------------------------
# Benchmarking SCDC
# ------------------------------------------------------------------------------

n_runs <- 3
runtimes <- numeric(n_runs)
memories <- numeric(n_runs)

cat(sprintf("Starting Benchmark (%d runs)...\n", n_runs))

# Ensure sample_id exists for SCDC_basis_ONE and SCDC_prop_ONE
pData(eset)$sample_id <- "sample1"

for (i in 1:n_runs) {
  cat(sprintf("\nRun %d/%d...\n", i, n_runs))

  # 1. Warm-up / Clean-up
  gc()

  # 2. Start Measurement
  start_time <- Sys.time()
  gc(reset = TRUE) # Reset max used stats


  # 3. Execution
  # SCDC_prop_ONE will calculate the basis internally since weight.basis = TRUE
  fit <- SCDC_prop_ONE(
    bulk.eset = bulk_eset,
    sc.eset = eset,
    ct.varname = "cellType",
    sample = "sample_id",
    ct.sub = ct.sub,
    weight.basis = TRUE
  )

  # 4. End Measurement
  end_time <- Sys.time()
  gc_stats <- gc()

  # 5. Record Results
  duration_secs <- as.numeric(end_time - start_time, units = "secs")

  # Extract Peak Memory (Max Used)
  # gc() returns columns: used, gc trigger, max used.
  # Columns 2, 4, 6 are in Mb (usually).
  # We use the 6th column (Max Used Mb) if available.
  peak_mem_mb <- NA
  if (ncol(gc_stats) >= 6) {
    peak_mem_mb <- sum(gc_stats[, 6])
  } else {
    # Fallback: estimate from used if max not available (unlikely in modern R)
    peak_mem_mb <- sum(gc_stats[, 2])
  }

  runtimes[i] <- duration_secs
  memories[i] <- peak_mem_mb

  cat(sprintf("  Run %d: %.4fs, %.2f MB\n", i, duration_secs, peak_mem_mb))
}

cat("\n--- Benchmark Results Summary ---\n")
res_df <- data.frame(Run = 1:n_runs, Time_s = runtimes, Mem_MB = memories)
print(res_df)

cat(sprintf("\nAverage Runtime: %.4f seconds (sd: %.4f)\n", mean(runtimes), sd(runtimes)))
cat(sprintf("Average Peak Memory: %.2f MB (sd: %.2f)\n", mean(memories), sd(memories)))

# Save results to CSV
csv_filename <- "scdc_benchmark_results.csv"
write.csv(res_df, csv_filename, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", csv_filename))

