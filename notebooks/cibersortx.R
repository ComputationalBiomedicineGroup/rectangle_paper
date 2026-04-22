#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
data_dir <- file.path(repo_root, "data")
fig3_dir <- file.path(data_dir, "fig3")
dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)

hao_default <- file.path(data_dir, "Hao_sce.rds")
bulks_csv_default <- file.path(data_dir, "bulks_finotello_tpm.csv")
bulks_rds_default <- file.path(data_dir, "bulks_finotello_tpm.rds")

hao_path <- Sys.getenv("HAO_SCE_PATH", unset = hao_default)
bulks_path <- Sys.getenv(
  "BULKS_PATH",
  unset = if (file.exists(bulks_csv_default)) bulks_csv_default else bulks_rds_default
)

cran_repo <- getOption("repos")[["CRAN"]]
if (is.null(cran_repo) || !nzchar(cran_repo) || cran_repo %in% c("@CRAN", "@CRAN@")) {
  cran_repo <- "https://cloud.r-project.org"
}

ensure_package <- function(pkg, archive_version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = cran_repo)
  }
  if (!requireNamespace(pkg, quietly = TRUE) && !is.null(archive_version)) {
    archive_url <- sprintf("%s/src/contrib/Archive/%s/%s_%s.tar.gz", cran_repo, pkg, pkg, archive_version)
    install.packages(archive_url, repos = NULL, type = "source")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' could not be loaded.", pkg))
  }
}

ensure_package("Matrix", archive_version = "1.6-5")
ensure_package("processx")
library(Matrix)
library(processx)

to_matrix_like <- function(x) {
  if (inherits(x, "Matrix") || is.matrix(x)) {
    return(x)
  }
  as.matrix(x)
}

stratified_subset_indices <- function(labels, max_cells, seed = 42L) {
  labels <- as.character(labels)
  n_total <- length(labels)
  if (n_total <= max_cells) {
    return(seq_len(n_total))
  }
  set.seed(as.integer(seed))
  idx_by_type <- split(seq_along(labels), labels)
  chosen <- integer(0)
  for (ct in names(idx_by_type)) {
    idx <- idx_by_type[[ct]]
    n_pick <- max(1L, floor(as.numeric(max_cells) * length(idx) / n_total))
    n_pick <- min(n_pick, length(idx))
    chosen <- c(chosen, sample(idx, n_pick))
  }
  if (length(chosen) > max_cells) {
    chosen <- sample(chosen, max_cells)
  } else if (length(chosen) < max_cells) {
    remaining <- setdiff(seq_len(n_total), chosen)
    if (length(remaining) > 0) {
      chosen <- c(chosen, sample(remaining, min(max_cells - length(chosen), length(remaining))))
    }
  }
  sort(unique(chosen))
}

write_tsv_matrix <- function(mat, file_path, chunk_rows = 250L) {
  mat <- to_matrix_like(mat)
  gene_ids <- rownames(mat)
  sample_ids <- colnames(mat)
  if (is.null(gene_ids) || is.null(sample_ids)) {
    stop("Matrix must have rownames and colnames before writing TSV.")
  }
  con <- file(file_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(c("GeneSymbol", sample_ids), collapse = "\t"), con = con)

  n <- nrow(mat)
  starts <- seq.int(1L, n, by = as.integer(chunk_rows))
  for (start in starts) {
    end <- min(start + as.integer(chunk_rows) - 1L, n)
    block <- mat[start:end, , drop = FALSE]
    if (inherits(block, "sparseMatrix")) {
      block <- as.matrix(block)
    }
    for (i in seq_len(nrow(block))) {
      vals <- sprintf("%.10g", block[i, ])
      writeLines(paste(c(gene_ids[start + i - 1L], vals), collapse = "\t"), con = con)
    }
  }
}

parse_mem_usage_to_mib <- function(mem_usage) {
  raw <- trimws(strsplit(mem_usage, "/", fixed = TRUE)[[1]][1])
  if (!nzchar(raw)) {
    return(0)
  }
  num <- gsub("[^0-9.]", "", raw)
  unit <- tolower(gsub("[0-9.[:space:]]", "", raw))
  if (!nzchar(num)) {
    return(0)
  }
  value <- as.numeric(num)
  factors <- c(
    b = 1 / (1024^2),
    kb = 1000 / (1024^2),
    mb = 1000^2 / (1024^2),
    gb = 1000^3 / (1024^2),
    kib = 1 / 1024,
    mib = 1,
    gib = 1024
  )
  value * ifelse(unit %in% names(factors), factors[[unit]], 1 / (1024^2))
}

run_docker_with_memory_sampling <- function(docker_args, container_name, stats_interval = 0.5) {
  proc <- process$new("docker", docker_args, stdout = "|", stderr = "|", cleanup = TRUE)
  t0 <- proc.time()[["elapsed"]]
  peak_mib <- 0

  repeat {
    if (!proc$is_alive()) {
      break
    }
    stats_out <- tryCatch(
      system2("docker", c("stats", "--no-stream", "--format", "{{.MemUsage}}", container_name), stdout = TRUE, stderr = FALSE),
      warning = function(...) character(0),
      error = function(...) character(0)
    )
    if (length(stats_out) > 0 && nzchar(stats_out[1])) {
      peak_mib <- max(peak_mib, parse_mem_usage_to_mib(stats_out[1]))
    }
    Sys.sleep(stats_interval)
  }

  status <- proc$get_exit_status()
  if (is.null(status)) {
    status <- proc$wait(timeout = Inf)
  }
  elapsed_s <- proc.time()[["elapsed"]] - t0
  stdout_txt <- paste(proc$read_all_output_lines(), collapse = "\n")
  stderr_txt <- paste(proc$read_all_error_lines(), collapse = "\n")

  list(status = status, elapsed_s = elapsed_s, peak_mib = peak_mib, stdout = stdout_txt, stderr = stderr_txt)
}

run_checked_step <- function(step_name, docker_args, container_name, stats_interval, expected_files = character(0)) {
  cat(sprintf("  Step: %s\n", step_name))
  res <- run_docker_with_memory_sampling(
    docker_args = docker_args,
    container_name = container_name,
    stats_interval = stats_interval
  )

  if (!identical(as.integer(res$status), 0L)) {
    if (nzchar(res$stdout)) {
      cat(res$stdout, "\n")
    }
    if (nzchar(res$stderr)) {
      cat(res$stderr, "\n")
    }
    stop(sprintf("CIBERSORTx step '%s' failed with exit code %s.", step_name, as.character(res$status)))
  }

  missing <- expected_files[!file.exists(expected_files)]
  if (length(missing) > 0 && length(missing) == length(expected_files)) {
    if (nzchar(res$stdout)) {
      cat(res$stdout, "\n")
    }
    if (nzchar(res$stderr)) {
      cat(res$stderr, "\n")
    }
    stop(
      sprintf(
        "CIBERSORTx step '%s' finished without required output files:\n%s",
        step_name,
        paste(missing, collapse = "\n")
      )
    )
  }

  res
}

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
if (grepl("\\.csv$", bulks_path, ignore.case = TRUE)) {
  bulks <- utils::read.csv(bulks_path, row.names = 1, check.names = FALSE)
} else if (grepl("\\.rds$", bulks_path, ignore.case = TRUE)) {
  bulks <- readRDS(bulks_path)
} else {
  stop(sprintf("Unsupported BULKS_PATH extension for: %s", bulks_path))
}

sc_counts <- counts
bulk_mat <- to_matrix_like(bulks)

if (is.null(rownames(sc_counts))) {
  rownames(sc_counts) <- genes
}
if (is.null(colnames(sc_counts))) {
  colnames(sc_counts) <- cellnames
}

# Normalize to CPM using full library size before gene subsetting.
cat("Normalizing single-cell counts to CPM...\n")
cell_lib_sizes <- Matrix::colSums(sc_counts)
cell_lib_sizes[cell_lib_sizes == 0] <- 1  # guard against empty cells
sc_counts <- Matrix::t(Matrix::t(sc_counts) / cell_lib_sizes) * 1e6

max_cells <- 80000
subset_idx <- stratified_subset_indices(annotations, max_cells = max_cells, seed = 42L)
sc_counts <- sc_counts[, subset_idx, drop = FALSE]
annotations <- annotations[subset_idx]
cat(sprintf("Using %d single cells for CIBERSORTx signature generation (target=%d).\n", ncol(sc_counts), max_cells))

common_genes <- intersect(rownames(bulk_mat), rownames(sc_counts))
if (length(common_genes) == 0) {
  stop("No overlapping genes between single-cell data and bulk data.")
}
bulk_mat <- bulk_mat[common_genes, , drop = FALSE]
sc_counts <- sc_counts[common_genes, , drop = FALSE]

signature_path <- file.path(fig3_dir, sprintf("cibersortx_signature_max%d.txt", ncol(sc_counts)))
mixture_path <- file.path(fig3_dir, "cibersortx_mixture_full.txt")
reuse_inputs <- tolower(Sys.getenv("CIBERSORTX_REUSE_INPUTS", unset = "true")) %in% c("1", "true", "yes", "y")
inputs_exist <- all(file.exists(c(signature_path, mixture_path)))

if (reuse_inputs && inputs_exist) {
  cat("Reusing existing CIBERSORTx input files.\n")
} else {
  cat("Writing CIBERSORTx signature (refsample) matrix...\n")
  # CIBERSORTx requires cell type labels as column names; no separate phenoclasses file.
  sc_for_cibx <- sc_counts
  colnames(sc_for_cibx) <- annotations
  write_tsv_matrix(sc_for_cibx, signature_path)
  cat("Writing CIBERSORTx mixture matrix...\n")
  write_tsv_matrix(bulk_mat, mixture_path)
}

# CIBERSORTx binary mount checks are strict; stage files at /src/data root.
# Use .txt extension — CIBERSORTx expects .txt, not .tsv.
staged_refsample <- file.path(data_dir, sprintf("cibx_refsample_max%d.txt", ncol(sc_counts)))
staged_mixture <- file.path(data_dir, "cibx_mixture_full.txt")
ok_ref <- file.copy(signature_path, staged_refsample, overwrite = TRUE)
ok_mix <- file.copy(mixture_path, staged_mixture, overwrite = TRUE)
if (!all(c(ok_ref, ok_mix))) {
  stop("Failed to stage CIBERSORTx input files into data/ for Docker execution.")
}

docker_bin <- Sys.which("docker")
if (!nzchar(docker_bin)) {
  stop("Docker is required, but 'docker' was not found in PATH.")
}
docker_ok <- suppressWarnings(system2(docker_bin, c("info"), stdout = FALSE, stderr = FALSE) == 0L)
if (!docker_ok) {
  stop("Docker command is available but not usable (docker info failed). Start Docker and retry.")
}

username <- ""
api_key <- ""
token <- Sys.getenv("CIBERSORTX_TOKEN", unset = "")
if (!nzchar(token) && nzchar(api_key)) {
  if (grepl(":", api_key, fixed = TRUE)) {
    parts <- strsplit(api_key, ":", fixed = TRUE)[[1]]
    if (!nzchar(username) && length(parts) >= 1) {
      username <- parts[1]
    }
    if (length(parts) >= 2) {
      token <- paste(parts[-1], collapse = ":")
    }
  } else {
    token <- api_key
  }
}
if (!nzchar(username)) {
  stop("Set CIBERSORTX_USERNAME before running this benchmark.")
}
if (!nzchar(token)) {
  stop("Set CIBERSORTX_TOKEN (or CIBERSORTX_API_KEY) before running this benchmark.")
}

cibersortx_image <- Sys.getenv("CIBERSORTX_DOCKER_IMAGE", unset = "cibersortx/fractions")
default_platform <- if (grepl("arm|aarch64", tolower(Sys.info()[["machine"]]))) "linux/amd64" else ""
docker_platform <- Sys.getenv("CIBERSORTX_DOCKER_PLATFORM", unset = default_platform)
n_runs <- 3
n_runs <- as.integer(Sys.getenv("CIBERSORTX_N_RUNS", unset = as.character(n_runs)))
stats_interval <- as.numeric(Sys.getenv("CIBERSORTX_DOCKER_STATS_INTERVAL", unset = "0.5"))

runtime_records <- list()
memory_records <- list()
n_cells <- ncol(sc_counts)
cat(sprintf("Benchmarking CIBERSORTx with %d-cell signature subset.\n", n_cells))

for (run in seq_len(n_runs)) {
  run_outdir <- file.path(fig3_dir, sprintf("cibersortx_run_%d", run))
  dir.create(run_outdir, recursive = TRUE, showWarnings = FALSE)
  container_name_base <- sprintf("cibx-bench-%d-%d", run, as.integer(Sys.time()))

  rel_data <- function(path) {
    gsub(
      "^/",
      "",
      sub(
        paste0("^", normalizePath(data_dir, mustWork = FALSE), "/?"),
        "",
        normalizePath(path, mustWork = FALSE)
      )
    )
  }

  signature_in <- basename(staged_refsample)   # binary resolves relative to /src/data
  mixture_in <- basename(staged_mixture)
  outdir_in <- "/src/outdir"
  run_sigmatrix_base <- sprintf("cibersortx_sigmatrix_run_%d.txt", run)
  run_sigmatrix_path <- file.path(run_outdir, run_sigmatrix_base)
  run_sigmatrix_in_frac <- run_sigmatrix_base   # will be copied to data_dir before frac step
  if (file.exists(run_sigmatrix_path)) {
    file.remove(run_sigmatrix_path)
  }
  fractions_out <- file.path(run_outdir, "CIBERSORTx_Results.txt")
  fractions_out_alt <- file.path(run_outdir, "CIBERSORTx_Adjusted.txt")
  if (file.exists(fractions_out)) {
    file.remove(fractions_out)
  }

  docker_prefix <- c("run", "--rm")
  if (nzchar(docker_platform)) {
    docker_prefix <- c(docker_prefix, "--platform", docker_platform)
  }
  docker_prefix <- c(
    docker_prefix,
    "-v", sprintf("%s:/src/data", normalizePath(data_dir, mustWork = FALSE)),
    "-v", sprintf("%s:/src/outdir", normalizePath(run_outdir, mustWork = FALSE))
  )

  cat(sprintf("Run %d/%d...\n", run, n_runs))
  gc()

  # 1) Build signature matrix from single-cell reference.
  build_args <- c(
    docker_prefix,
    "--name", paste0(container_name_base, "-build"),
    cibersortx_image,
    "--username", username,
    "--token", token,
    "--single_cell", "TRUE",
    "--refsample", signature_in,
    "--G.min", "300",
    "--G.max", "500",
    "--q.value", "0.01",
    "--filter", "FALSE",
    "--k.max", "999",
    "--remake", "FALSE",
    "--replicates", "5",
    "--sampling", "0.5",
    "--fraction", "0.75",
    "--outdir", outdir_in
  )
  build_res <- run_checked_step(
    step_name = "build_signature",
    docker_args = build_args,
    container_name = paste0(container_name_base, "-build"),
    stats_interval = stats_interval
  )
  cat(sprintf("    build_signature: %.4fs, %.2f MiB\n", build_res$elapsed_s, build_res$peak_mib))

  # Locate the sigmatrix written by CIBERSORTx (named like *inferred_refsample.bm.K*.txt)
  generated_sig <- list.files(run_outdir, pattern = "bm\\.K[0-9]+\\.txt$", full.names = TRUE)
  if (length(generated_sig) == 0) {
    cat("Warning: could not locate sigmatrix in", run_outdir, "- skipping rename.\n")
  } else {
    file.copy(generated_sig[1], run_sigmatrix_path, overwrite = TRUE)
  }

  # 2) Run fractions — pass the sigmatrix by its original CIBERSORTx output name via /src/outdir.
  actual_sigmatrix_name <- if (length(generated_sig) > 0) basename(generated_sig[1]) else run_sigmatrix_base
  # Stage into data_dir so fractions can read it via /src/data
  staged_sigmatrix <- file.path(data_dir, actual_sigmatrix_name)
  if (length(generated_sig) > 0) file.copy(generated_sig[1], staged_sigmatrix, overwrite = TRUE)

  frac_args <- c(
    docker_prefix,
    "--name", paste0(container_name_base, "-frac"),
    cibersortx_image,
    "--username", username,
    "--token", token,
    "--mixture", mixture_in,
    "--sigmatrix", actual_sigmatrix_name,
    "--perm", "100",
    "--rmbatchBmode", "FALSE",
    "--rmbatchSmode", "FALSE",
    "--QN", "FALSE",
    "--absolute", "FALSE",
    "--abs_method", "sig.score",
    "--outdir", outdir_in
  )
  frac_res <- run_checked_step(
    step_name = "run_fractions",
    docker_args = frac_args,
    container_name = paste0(container_name_base, "-frac"),
    stats_interval = stats_interval,
    expected_files = c(fractions_out, fractions_out_alt)
  )

  # Outputs land directly in run_outdir via /src/outdir mount.

  elapsed_total <- build_res$elapsed_s + frac_res$elapsed_s
  peak_mib <- max(build_res$peak_mib, frac_res$peak_mib)
  peak_gib <- peak_mib / 1024
  runtime_records[[length(runtime_records) + 1L]] <- data.frame(
    n_cells = n_cells,
    run = run,
    time = elapsed_total
  )
  memory_records[[length(memory_records) + 1L]] <- data.frame(
    n_cells = n_cells,
    run = run,
    memory_mib = peak_mib,
    memory_gib = peak_gib
  )
  cat(sprintf("  %.4fs, %.2f MiB (%.3f GiB)\n", elapsed_total, peak_mib, peak_gib))
}

runtime_df <- do.call(rbind, runtime_records)
memory_df <- do.call(rbind, memory_records)

runtime_file <- file.path(fig3_dir, "cibersortx_runtime_80k.csv")
memory_file <- file.path(fig3_dir, "cibersortx_memory_80k.csv")
write.csv(runtime_df, runtime_file, row.names = FALSE)
write.csv(memory_df, memory_file, row.names = FALSE)

runtime_summary <- c(mean = mean(runtime_df$time), sd = sd(runtime_df$time))
memory_summary_mib <- c(mean = mean(memory_df$memory_mib), sd = sd(memory_df$memory_mib))
memory_summary_gib <- c(mean = mean(memory_df$memory_gib), sd = sd(memory_df$memory_gib))

cat("\n--- CIBERSORTx Benchmark Summary ---\n")
cat(sprintf("Runtime mean=%.4fs sd=%.4f\n", runtime_summary["mean"], runtime_summary["sd"]))
cat(sprintf("Memory mean=%.2f MiB sd=%.2f\n", memory_summary_mib["mean"], memory_summary_mib["sd"]))
cat(sprintf("Memory mean=%.3f GiB sd=%.3f\n", memory_summary_gib["mean"], memory_summary_gib["sd"]))
cat(sprintf("Saved runtime: %s\n", runtime_file))
cat(sprintf("Saved memory: %s\n", memory_file))
