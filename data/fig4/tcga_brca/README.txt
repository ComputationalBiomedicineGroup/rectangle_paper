TCGA-BRCA deconvolved with wu (rectangle revision), run 21-22 Oct 2026
Pipeline: deconvBench, branch refactor, workflow main; configs pipeline/configs/tcga_brca.config and
tcga_brca_dwls.config (DWLS with one core); methods: autogenes, bayesprism, bisque, cibersortx, dwls,
music, scaden, scdc, rectangle (rectanglepy 1.7.0), secret

reference/wu-normal-tcga-full/   wu at the normal annotation level, all cells (66,559 cells, 26 donors,
                                 18 cell types, 23,542 genes): without the cell types Cycling T-cells,
                                 Cycling_Myeloid, Cancer Cycling, Cycling PVL, Myoepithelial, Luminal
                                 Progenitors, Mature Luminal, Cancer Her2 SC, Cancer LumB SC, Cancer Basal SC,
                                 Cancer LumA SC, and without the genes matching ^RP[SL]|^RP\d|^RPSA;
                                 matrix_norm_counts.rds is the CPM of the remaining genes.
prepare_tcga_brca.R              the script that built the reference and the bulk dataset
results/<method>_wu-normal-tcga-full_<sc norm>_tcga-brca_<bulk norm>_ct0_rep0/
                                 deconvolution.rds: estimated fractions, 1,091 samples x 18 cell types
                                 (rectangle and secret also an unknown fraction); signature and runtimes.
                                 No metrics: the bulks have no ground truth. The trained Scaden network
                                 (scaden_*/model) is not included.
Bulk: 1,091 primary tumours (TCGA barcodes ending -01), 20,472 genes, TPM and counts as provided
(counts are not integer, TPM columns sum to 0.57-0.98 million); not included, in
/nfs/data/omnideconv_benchmarking_clean/data/TCGA-BRCA/tcga-brca/.
Cancer epithelial cells are not in the reference: rectangle and secret report them as unknown, the
other methods distribute them over the 18 cell types.
