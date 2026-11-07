# SECRET (github.com/YyLu5/SECRET) on the Figure 2 bulks with the references Rectangle uses (HaoSub, Tabula Muris).
#  - signature: marker genes per cell type with Seurat FindAllMarkers (bimod, min.pct 0.25, logFC 0.25) on the
#    log-normalised reference; gene weight 1 / number of cell types the gene is a marker of (as in the SECRET vignette)
#  - reference: mean expression per cell type from SECRET::scRefer, restricted to the marker genes
#  - deconvolution: SECRET on marker genes shared by bulk and reference, both rescaled to CPM inside SECRET
#    (yNorm = bNorm = 'cpm'), with the extra unknown component; with the argument no_unknown, without it
#    (fractions sum to 1).
# With the argument spillover, it runs the purified spillover bulks of figure_3.ipynb (Hao pseudo-bulks, Hoek,
# Linsley) with the HaoSub reference instead, reusing its cached marker genes.
# Inputs from scripts/export_secret_fig2_inputs.py.
# Usage: Rscript scripts/run_secret_fig2.R [spillover] [no_unknown] [dataset ...]
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Biobase)
  library(SECRET)
})
set.seed(1)
with_unknown <- !("no_unknown" %in% commandArgs(trailingOnly = TRUE))
suffix <- if (with_unknown) "" else "_no_unknown"
spillover <- "spillover" %in% commandArgs(trailingOnly = TRUE)
only <- setdiff(commandArgs(trailingOnly = TRUE), c("no_unknown", "spillover"))  # optional subset of bulk datasets
root <- "output/secret_fig2"
inputs <- file.path(root, "inputs")

datasets <- list(
  hao_sub = c("finotello", "hoek", "morandini", "dietrich"),
  tabula_muris = c("chen", "petitprez")
)
if (spillover) datasets <- list(hao_sub = c("spillover_hao", "spillover_hoek", "spillover_linsley"))

load_reference <- function(name) {
  folder <- file.path(inputs, "references", name)
  counts <- as(readMM(file.path(folder, "counts.mtx")), "CsparseMatrix")
  cells <- read.csv(file.path(folder, "cells.csv"), check.names = FALSE)
  rownames(counts) <- read.csv(file.path(folder, "genes.csv"))$gene
  colnames(counts) <- as.character(cells$cell)
  cell_type <- setNames(cells$cell_type, colnames(counts))
  list(counts = counts, cell_type = cell_type)
}

find_markers <- function(name, reference) {
  path <- file.path(root, paste0("markers_", name, ".csv"))
  if (file.exists(path)) return(read.csv(path, check.names = FALSE))
  seurat <- CreateSeuratObject(reference$counts)
  seurat <- NormalizeData(seurat, verbose = FALSE)
  Idents(seurat) <- factor(reference$cell_type[colnames(seurat)])
  markers <- FindAllMarkers(seurat, test.use = "bimod", min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  write.csv(markers, path, row.names = FALSE)
  markers
}

# SECRET::scRefer averages per sample and needs at least two samples (it calls apply() on the per-sample
# columns of a cell type). The references have no donor column, so the cells are passed twice as two
# identical samples; the result is then the mean count per cell type, which is checked below.
build_reference <- function(reference, genes) {
  counts <- as.matrix(reference$counts[genes, , drop = FALSE])
  meta <- data.frame(cell_type = rep(reference$cell_type, 2), sample = rep(c("s1", "s2"), each = ncol(counts)))
  doubled <- cbind(counts, counts)
  colnames(doubled) <- rownames(meta) <- c(paste0("s1_", colnames(counts)), paste0("s2_", colnames(counts)))
  eset <- ExpressionSet(assayData = doubled, phenoData = AnnotatedDataFrame(meta))
  ref <- scRefer(eset, ct_var = "cell_type", sample_var = "sample")
  means <- sapply(colnames(ref), function(ct) rowMeans(counts[rownames(ref), reference$cell_type == ct, drop = FALSE]))
  stopifnot(isTRUE(all.equal(ref, means, check.attributes = FALSE)))
  ref
}

for (ref_name in names(datasets)) {
  if (length(only) && !any(datasets[[ref_name]] %in% only)) next
  reference <- load_reference(ref_name)
  markers <- find_markers(ref_name, reference)
  cat(ref_name, ":", nrow(markers), "marker rows,", length(unique(markers$gene)), "genes\n")
  frequency <- table(markers$gene)
  for (dataset in datasets[[ref_name]]) {
    if (length(only) && !(dataset %in% only)) next
    bulk <- as.matrix(read.csv(file.path(inputs, "bulks", paste0(dataset, ".csv")), row.names = 1, check.names = FALSE))
    genes <- Reduce(intersect, list(names(frequency), rownames(bulk), rownames(reference$counts)))
    genes <- genes[rowSums(reference$counts[genes, , drop = FALSE]) > 0]  # scRefer drops unexpressed genes
    ref <- build_reference(reference, genes)
    genes <- rownames(ref)
    weights <- as.numeric(1 / frequency[genes])
    cat(" ", dataset, ":", length(genes), "genes,", ncol(bulk), "samples\n")
    # SECRET fits every sample on its own, so the samples are split over cores
    estimates <- do.call(rbind, parallel::mclapply(colnames(bulk), function(sample) {
      SECRET(bulk[genes, sample, drop = FALSE], ref[genes, ], withUnknown = with_unknown, w = weights,
             yNorm = "cpm", bNorm = "cpm")[[1]]
    }, mc.cores = parallel::detectCores()))
    stopifnot(identical(rownames(estimates), colnames(bulk)))
    colnames(estimates)[colnames(estimates) == "unknown"] <- "Unknown"
    write.csv(estimates, file.path(root, paste0("secret_estimations_", dataset, suffix, ".csv")))
  }
}
