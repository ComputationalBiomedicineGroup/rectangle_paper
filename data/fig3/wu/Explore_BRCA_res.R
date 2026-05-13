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
  "T-cells" = "T-cells",
  "T.cells.CD4." = "T cells CD4",
  "T.cells.CD8." = "T cells CD8",
  "CAFs" = "CAFs",
  "CAFs.myCAF.like" = "CAFs myCAF like",
  "CAFs.MSC.iCAF.like" = "CAFs MSC iCAF like"
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

# Option 2: Only comparisons vs reference (e.g., D)
# comparisons_vs_D <- list(c("D", "F"), c("D", "IE/F"), c("D", "IE"))

# Create list of plots with statistics
plot_list <- list()

for (ct in cell_types) {
  
  p <- ggplot(plot_data %>% filter(cell_type == ct), 
              aes(x = subtype, y = proportion, fill = subtype)) +
    geom_boxplot(alpha = 0.7) +
    # Add pairwise Wilcoxon tests (uncomment one option below)
    # Option A: Show all pairwise comparisons
    # stat_compare_means(comparisons = all_comparisons, 
    #                    method = "wilcox.test",
    #                    p.adjust.method = "BH",
    #                    label = "p.signif",
    #                    hide.ns = TRUE) +
    # Option B: Show only significant comparisons with exact p-values
    stat_compare_means(comparisons = all_comparisons,
                       method = "wilcox.test",
                       p.adjust.method = "BH",
                       label = "p.signif",
                       hide.ns = TRUE,
                       size = 3) +
    labs(x = "Subtype", y = "Proportion", title = cell_type_labels[[ct]]) +
    theme_minimal() +
    theme(legend.position = "none")
  
  plot_list[[ct]] <- p
}

# Save to PDF
pdf("boxplots_by_celltype.pdf", width = 12, height = 8)
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
  
  # Pairwise Wilcoxon tests with BH correction
  pw_test <- pairwise.wilcox.test(ct_data$proportion, 
                                  ct_data$subtype, 
                                  p.adjust.method = "BH")
  
  # Extract pairwise p-values
  pw_pvals <- as.data.frame(as.table(pw_test$p.value))
  pw_pvals <- pw_pvals %>% 
    filter(!is.na(Freq)) %>%
    mutate(comparison = paste(Var1, "vs", Var2),
           cell_type = ct,
           p_adjusted = Freq) %>%
    select(cell_type, comparison, p_adjusted)
  
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
  mutate(significance = case_when(
    p_adjusted < 0.001 ~ "***",
    p_adjusted < 0.01 ~ "**",
    p_adjusted < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Save statistics
write.csv(stats_table, "statistical_tests_results.csv", row.names = FALSE)

# Print summary for CAF subtypes
cat("\n=== CAF subtype statistics ===\n")
stats_table %>% 
  filter(grepl("CAF", cell_type)) %>%
  print()

## ------------------------------
## Verify ggpubr statistics
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

# 3. Pairwise Wilcoxon tests WITH BH correction
cat("\n=== Pairwise Wilcoxon (BH correction) ===\n")
pw_corrected <- pairwise.wilcox.test(ct_data$proportion, 
                                     ct_data$subtype, 
                                     p.adjust.method = "BH")
print(pw_corrected)

# 4. Manual check of a specific comparison
cat("\n=== Manual check: D vs F ===\n")
d_vals <- ct_data$proportion[ct_data$subtype == "D"]
f_vals <- ct_data$proportion[ct_data$subtype == "F"]
manual_test <- wilcox.test(d_vals, f_vals)
print(manual_test)

# 5. Check what ggpubr computes internally
cat("\n=== ggpubr internal computation ===\n")
ggpubr_stats <- compare_means(proportion ~ subtype, 
                              data = ct_data,
                              method = "wilcox.test",
                              p.adjust.method = "BH")
print(as.data.frame(ggpubr_stats))

