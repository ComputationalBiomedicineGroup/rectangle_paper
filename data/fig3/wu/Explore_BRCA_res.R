rm(list = ls())

## ------------------------------
## Packages
## ------------------------------
library(ggplot2)

## ------------------------------
## Paths
## ------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
}

## ------------------------------
## Settings
## ------------------------------
subtype_levels <- c("D", "F", "IE/F", "IE")

cell_types_to_plot <- c(
  "T-cells",
  make.names("T cells CD4+"),
  make.names("T cells CD8+"),
  "CAFs",
  make.names("CAFs myCAF-like"),
  make.names("CAFs MSC iCAF-like")
)

cell_type_labels <- c(
  "T-cells" = "T-cells",
  "T.cells.CD4." = "T.cells.CD4.",
  "T.cells.CD8." = "T.cells.CD8.",
  "CAFs" = "CAFs",
  "CAFs.myCAF.like" = "CAFs.myCAF.like",
  "CAFs.MSC.iCAF.like" = "CAFs.MSC.iCAF.like"
)

all_comparisons <- data.frame(
  group1 = c("D", "D", "D", "F", "F", "IE/F"),
  group2 = c("F", "IE/F", "IE", "IE/F", "IE", "IE"),
  stringsAsFactors = FALSE
)

## ------------------------------
## Data loading and preprocessing
## ------------------------------
data <- read.csv("wu_estimations.csv", row.names = 1, check.names = TRUE)

cellmap <- read.csv2("Cell_type_map.csv", fileEncoding = "UTF-8-BOM",
                     stringsAsFactors = FALSE, check.names = FALSE)

for (celltype in unique(cellmap$celltype_major)) {
  sel_types <- make.names(cellmap$celltype_minor[cellmap$celltype_major == celltype])
  sel_types <- intersect(sel_types, colnames(data))
  if (length(sel_types) > 0) {
    data[[celltype]] <- rowSums(data[, sel_types, drop = FALSE], na.rm = TRUE)
  }
}

subtypes <- read.csv2("BRCA_immune_subtypes.csv", fileEncoding = "UTF-8-BOM",
                      row.names = 1, stringsAsFactors = FALSE,
                      check.names = FALSE)

data$subtype <- subtypes$Class[match(substr(rownames(data), 1, 12),
                                     rownames(subtypes))]
data <- data[!is.na(data$subtype), , drop = FALSE]
data$subtype <- factor(data$subtype, levels = subtype_levels)

missing_cell_types <- setdiff(cell_types_to_plot, colnames(data))
if (length(missing_cell_types) > 0) {
  stop("Missing requested cell type columns: ",
       paste(missing_cell_types, collapse = ", "))
}

plot_data <- do.call(rbind, lapply(cell_types_to_plot, function(ct) {
  data.frame(
    sample_id = rownames(data),
    subtype = data$subtype,
    cell_type = factor(ct, levels = cell_types_to_plot),
    proportion = data[[ct]],
    stringsAsFactors = FALSE
  )
}))

plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_types_to_plot)

## ------------------------------
## Statistical testing
## ------------------------------
format_p <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  }
  if (p < 2e-16) {
    return("p < 2e-16")
  }
  paste0("p = ", format.pval(p, digits = 2, eps = 2e-16))
}

signif_label <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.0001, "****",
                ifelse(p < 0.001, "***",
                       ifelse(p < 0.01, "**",
                              ifelse(p < 0.05, "*", "ns")))))
}

kw_rows <- list()
pairwise_rows <- list()

for (ct in cell_types_to_plot) {
  ct_data <- plot_data[plot_data$cell_type == ct, , drop = FALSE]

  kw_test <- kruskal.test(proportion ~ subtype, data = ct_data)
  kw_rows[[ct]] <- data.frame(
    cell_type = ct,
    test = "Kruskal-Wallis",
    group1 = NA_character_,
    group2 = NA_character_,
    comparison = "overall",
    p_raw = kw_test$p.value,
    p_adjusted = NA_real_,
    stringsAsFactors = FALSE
  )

  raw_p <- numeric(nrow(all_comparisons))
  for (i in seq_len(nrow(all_comparisons))) {
    vals1 <- ct_data$proportion[ct_data$subtype == all_comparisons$group1[i]]
    vals2 <- ct_data$proportion[ct_data$subtype == all_comparisons$group2[i]]
    raw_p[i] <- wilcox.test(vals1, vals2, exact = FALSE)$p.value
  }

  pairwise_rows[[ct]] <- data.frame(
    cell_type = ct,
    test = "Wilcoxon rank-sum",
    group1 = all_comparisons$group1,
    group2 = all_comparisons$group2,
    comparison = paste(all_comparisons$group1, "vs", all_comparisons$group2),
    p_raw = raw_p,
    p_adjusted = p.adjust(raw_p, method = "BH"),
    stringsAsFactors = FALSE
  )
}

stats_table <- rbind(do.call(rbind, kw_rows), do.call(rbind, pairwise_rows))
stats_table$significance <- signif_label(ifelse(is.na(stats_table$p_adjusted),
                                                stats_table$p_raw,
                                                stats_table$p_adjusted))
write.csv(stats_table, "statistical_tests_results.csv", row.names = FALSE)

## ------------------------------
## Plot annotations
## ------------------------------
pairwise_stats <- do.call(rbind, pairwise_rows)
pairwise_stats$significance <- signif_label(pairwise_stats$p_adjusted)
pairwise_stats <- pairwise_stats[pairwise_stats$significance != "ns", , drop = FALSE]
pairwise_stats$xmin <- match(pairwise_stats$group1, subtype_levels)
pairwise_stats$xmax <- match(pairwise_stats$group2, subtype_levels)

kw_stats <- do.call(rbind, kw_rows)
kw_stats$label <- vapply(kw_stats$p_raw, format_p, character(1))

annotation_data <- data.frame()
kw_annotation_data <- data.frame()

for (ct in cell_types_to_plot) {
  ct_values <- plot_data$proportion[plot_data$cell_type == ct]
  y_min <- min(ct_values, na.rm = TRUE)
  y_max <- max(ct_values, na.rm = TRUE)
  y_range <- y_max - y_min
  step <- max(y_range * 0.09, y_max * 0.04, 0.01)

  ct_pairwise <- pairwise_stats[pairwise_stats$cell_type == ct, , drop = FALSE]
  if (nrow(ct_pairwise) > 0) {
    ct_pairwise$y.position <- y_max + step * seq_len(nrow(ct_pairwise)) + step
    annotation_data <- rbind(annotation_data, ct_pairwise)
  }

  ct_kw <- kw_stats[kw_stats$cell_type == ct, , drop = FALSE]
  ct_kw$x <- 0.82
  ct_kw$y <- y_max + step
  kw_annotation_data <- rbind(kw_annotation_data, ct_kw)
}

## ------------------------------
## Plot
## ------------------------------
subtype_colors <- c(
  "D" = "#F28E8B",
  "F" = "#9BC45C",
  "IE/F" = "#54C7CD",
  "IE" = "#C98AED"
)

p <- ggplot(plot_data, aes(x = subtype, y = proportion, fill = subtype)) +
  geom_boxplot(width = 0.72, alpha = 0.82, outlier.size = 1.4,
               outlier.alpha = 0.7, linewidth = 0.35) +
  geom_segment(data = annotation_data,
               aes(x = xmin, xend = xmax, y = y.position, yend = y.position),
               inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = annotation_data,
               aes(x = xmin, xend = xmin, y = y.position,
                   yend = y.position - 0.012),
               inherit.aes = FALSE, linewidth = 0.35) +
  geom_segment(data = annotation_data,
               aes(x = xmax, xend = xmax, y = y.position,
                   yend = y.position - 0.012),
               inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = annotation_data,
            aes(x = (xmin + xmax) / 2, y = y.position,
                label = significance),
            inherit.aes = FALSE, vjust = -0.25, size = 3.0) +
  geom_text(data = kw_annotation_data,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3.2) +
  facet_wrap(
    ~ cell_type,
    ncol = 3,
    scales = "free_y",
    labeller = as_labeller(cell_type_labels)
  ) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.22))) +
  labs(x = "Subtype", y = "Proportion") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 11),
    axis.title = element_text(size = 10),
    axis.text = element_text(color = "grey25"),
    panel.spacing = unit(1.2, "lines")
  )

ggsave("boxplots_by_celltype.pdf", p, width = 12, height = 7.2,
       units = "in", device = cairo_pdf)
ggsave("boxplots_by_celltype.png", p, width = 12, height = 7.2,
       units = "in", dpi = 300)

cat("\nSaved plots:\n")
cat("  boxplots_by_celltype.pdf\n")
cat("  boxplots_by_celltype.png\n")
cat("\nSaved statistics:\n")
cat("  statistical_tests_results.csv\n")
