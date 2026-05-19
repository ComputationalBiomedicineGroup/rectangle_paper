rm(list=ls())

## ------------------------------
## Packages
## ------------------------------
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(ggpubr)  # For stat_compare_means and easy significance annotation

p_to_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

## ------------------------------
## Data loading and preprocessing
## ------------------------------
# Load data
data <- data.frame(read_csv("wu_estimations.csv"),
                   row.names = 1)

# Add aggregated data
cellmap <- data.frame(read_delim("Cell_type_map.csv"))
for (celltype in cellmap$celltype_major) {
  sel_types <- make.names(cellmap$celltype_minor[cellmap$celltype_major == celltype])
  data[, celltype] <- apply(data[, sel_types], 1, sum)
}

# Add tumor immune subtypes
subtypes <- data.frame(read_delim("BRCA_immune_subtypes.csv"),
                       row.names = 1)
data$subtype <- subtypes$Class[match(substr(rownames(data), 1, 12), 
                                     rownames(subtypes))]

## ------------------------------
## Plots with statistical testing
## ------------------------------
# Remove rows with NA subtype
data <- data %>% filter(!is.na(subtype))

# Set subtype order
data$subtype <- factor(data$subtype, levels = c("D", "F", "IE/F", "IE"))
subtype_labels <- c(
  "D" = "Desert",
  "F" = "Fibrotic",
  "IE/F" = "Fibrotic /\nInflamed",
  "IE" = "Inflamed"
)

# Get cell type columns to plot, ordered as requested
cell_types <- c(
  "T-cells",
  "T.cells.CD4.",
  "T.cells.CD8.",
  "CAFs",
  "CAFs.myCAF.like",
  "CAFs.MSC.iCAF.like"
)

cell_type_labels <- c(
  "T-cells" = "T cells",
  "T.cells.CD4." = "CD4+ T cells",
  "T.cells.CD8." = "CD8+ T cells",
  "CAFs" = "CAFs",
  "CAFs.myCAF.like" = "myCAFs",
  "CAFs.MSC.iCAF.like" = "iCAFs"
)

# Reshape data to long format
plot_data <- data %>%
  pivot_longer(cols = all_of(cell_types),
               names_to = "cell_type",
               values_to = "proportion")

# Define pairwise comparisons of interest
# Option 1: All pairwise comparisons
all_comparisons <- list(
  c("D", "F"), c("D", "IE/F"), c("D", "IE"),
  c("F", "IE/F"), c("F", "IE"), c("IE/F", "IE")
)

all_comparisons_df <- data.frame(
  group1 = vapply(all_comparisons, `[`, character(1), 1),
  group2 = vapply(all_comparisons, `[`, character(1), 2)
)

# Option 2: Only comparisons vs reference (e.g., D)
# comparisons_vs_D <- list(c("D", "F"), c("D", "IE/F"), c("D", "IE"))

# Pre-compute pairwise statistics so plot annotations and exported results use
# the same BH-adjusted p-values. BH correction is applied globally across all
# 6 comparisons x 6 cell types.
pairwise_stats <- data.frame()

for (ct in cell_types) {
  
  ct_data <- plot_data %>% filter(cell_type == ct)
  
  p_raw <- mapply(function(group1, group2) {
    wilcox.test(
      ct_data$proportion[ct_data$subtype == group1],
      ct_data$proportion[ct_data$subtype == group2]
    )$p.value
  }, all_comparisons_df$group1, all_comparisons_df$group2)
  
  y_range <- range(ct_data$proportion, na.rm = TRUE)
  y_step <- max(diff(y_range) * 0.08, max(y_range) * 0.05, 0.01)
  
  pairwise_stats <- rbind(
    pairwise_stats,
    data.frame(
      cell_type = ct,
      group1 = all_comparisons_df$group1,
      group2 = all_comparisons_df$group2,
      comparison = paste(all_comparisons_df$group1, "vs", all_comparisons_df$group2),
      p = p_raw,
      y.position = max(y_range) + seq_along(p_raw) * y_step
    )
  )
}

pairwise_stats <- pairwise_stats %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  mutate(p.adj.signif = p_to_stars(p.adj))

# Create list of plots with statistics
plot_list <- list()

for (ct in cell_types) {
  ct_pairwise_stats <- pairwise_stats %>%
    filter(cell_type == ct, p.adj < 0.05)
  
  p <- ggplot(plot_data %>% filter(cell_type == ct), 
              aes(x = subtype, y = proportion, fill = subtype)) +
    geom_boxplot(alpha = 0.7) +
    labs(x = "Subtype", y = "Cell fractions (%)", title = cell_type_labels[[ct]]) +
    scale_x_discrete(labels = subtype_labels) +
    scale_y_continuous(labels = function(x) x * 100) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 9, face = "plain"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9)
    )
  
  if (nrow(ct_pairwise_stats) > 0) {
    p <- p +
      stat_pvalue_manual(ct_pairwise_stats,
                         label = "p.adj.signif",
                         xmin = "group1",
                         xmax = "group2",
                         y.position = "y.position",
                         tip.length = 0.01,
                         size = 3.5)
  }
  
  plot_list[[ct]] <- p
}

# Save to PDF
pdf("boxplots_by_celltype.pdf", width = 9, height = 6)
n_per_page <- 6
n_pages <- ceiling(length(plot_list) / n_per_page)

for (i in seq_len(n_pages)) {
  start_idx <- (i - 1) * n_per_page + 1
  end_idx <- min(i * n_per_page, length(plot_list))
  
  grid.arrange(grobs = plot_list[start_idx:end_idx], 
               ncol = 3, nrow = 2)
}
dev.off()

## ------------------------------
## Export statistics table
## ------------------------------
# Create a summary table with all statistics
stats_table <- data.frame()

for (ct in cell_types) {
  
  ct_data <- plot_data %>% filter(cell_type == ct)
  
  # Kruskal-Wallis test
  kw_test <- kruskal.test(proportion ~ subtype, data = ct_data)
  
  pw_pvals <- pairwise_stats %>%
    filter(cell_type == ct) %>%
    transmute(cell_type, comparison, p_adjusted = p.adj)
  
  # Add KW result
  kw_row <- data.frame(
    cell_type = ct,
    comparison = "Kruskal-Wallis (overall)",
    p_adjusted = kw_test$p.value
  )
  
  stats_table <- rbind(stats_table, kw_row, pw_pvals)
}

# Add significance stars
stats_table <- stats_table %>%
  mutate(significance = p_to_stars(p_adjusted))

# Save statistics
write.csv(stats_table, "statistical_tests_results.csv", row.names = FALSE)

# Print summary for CAF subtypes
cat("\n=== CAF subtype statistics ===\n")
stats_table %>% 
  filter(grepl("CAF", cell_type)) %>%
  print()

## ------------------------------
## Verify statistics
## ------------------------------

# Pick a specific cell type to check
ct_check <- "CAFs.myCAF.like"  # Change to your cell type of interest

ct_data <- plot_data %>% filter(cell_type == ct_check)

# 1. Kruskal-Wallis test (overall)
kw_test <- kruskal.test(proportion ~ subtype, data = ct_data)
cat("\n=== Kruskal-Wallis test for", ct_check, "===\n")
print(kw_test)

# 2. Pairwise Wilcoxon tests WITHOUT correction
cat("\n=== Pairwise Wilcoxon (no correction) ===\n")
pw_uncorrected <- pairwise.wilcox.test(ct_data$proportion, 
                                       ct_data$subtype, 
                                       p.adjust.method = "none")
print(pw_uncorrected)

# 3. Pairwise Wilcoxon tests WITH global BH correction
cat("\n=== Pairwise Wilcoxon (global BH correction across all cell types) ===\n")
pw_corrected <- pairwise_stats %>%
  filter(cell_type == ct_check) %>%
  select(group1, group2, p, p.adj, p.adj.signif)
print(pw_corrected)

# 4. Manual check of a specific comparison
cat("\n=== Manual check: D vs F ===\n")
d_vals <- ct_data$proportion[ct_data$subtype == "D"]
f_vals <- ct_data$proportion[ct_data$subtype == "F"]
manual_test <- wilcox.test(d_vals, f_vals)
print(manual_test)

# 5. Check ggpubr raw p-values
cat("\n=== ggpubr raw p-values ===\n")
ggpubr_stats <- compare_means(proportion ~ subtype, 
                              data = ct_data,
                              method = "wilcox.test",
                              p.adjust.method = "none")
print(as.data.frame(ggpubr_stats))
