# "braino_rectangle_estimations.csv"
# "braino_rectangle_err.csv"

rm(list=ls())

## ------------------------------
## Packages
## ------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

## ------------------------------
## Data loading and preprocessing
## ------------------------------

# Load data
data <- data.frame(read_csv("braino_rectangle_estimations.csv"),
                   row.names = 1)

# Create sample info table
split_info <- strsplit(rownames(data), "-")
metadata <- data.frame(replicate = paste0("R", sapply(split_info, `[`, 1)))
metadata$cell_line <- sapply(split_info, `[`, 2)
metadata$protocol <- sapply(split_info, `[`, 3)
metadata$day <- as.numeric(gsub("day", "", sapply(split_info, `[`, 4)))
rownames(metadata) <- rownames(data)

# Remove "pluripotent" samples
idrm <- which(metadata$protocol == "pluripotent")
metadata <- metadata[-idrm,]
data <- data[-idrm,]

## ------------------------------
## Define categories
## ------------------------------

cell_types <- colnames(data)
cell_lines <- unique(metadata$cell_line)
protocols <- unique(metadata$protocol)

cat("Cell types:", length(cell_types), "\n")
cat("Cell lines:", length(cell_lines), "\n")
cat("Protocols:", length(protocols), "\n")

## ------------------------------
## Define palettes
## ------------------------------

# Cell type palette (16 colors)
cell_type_palette <- c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
  "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
  "#cab2d6", "#6a3d9a", "#ffff99", "#b15928",
  "#8dd3c7", "#bebada", "#fb8072", "#80b1d3"
)
names(cell_type_palette) <- cell_types

# Cell line palette
cell_line_palette <- c(
  "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
  "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a",
  "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99"
)[1:length(cell_lines)]
names(cell_line_palette) <- cell_lines

## ------------------------------
## Plots: replicates vs. time
## ------------------------------

# Helper function to extract legend from a ggplot
get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) tmp$grobs[[leg]] else NULL
}

# By cell type
pdf("rectangle_by_celltype.pdf", width = 14, height = 9)

for (sel_cell_type in cell_types) {
  
  y_values <- data[, sel_cell_type]
  y_min <- min(y_values, na.rm = TRUE)
  y_max <- max(y_values, na.rm = TRUE)
  
  plot_list <- list()
  
  for (sel_protocol in protocols) {
    
    plot_data <- metadata %>%
      mutate(y_value = data[rownames(metadata), sel_cell_type]) %>%
      filter(protocol == sel_protocol)
    
    p <- ggplot(plot_data, aes(x = day, y = y_value, color = cell_line, 
                               group = interaction(cell_line, replicate))) +
      geom_line() +
      geom_point() +
      scale_color_manual(values = cell_line_palette) +
      labs(x = "Day", y = "Proportion", title = sel_protocol) +
      ylim(y_min, y_max) +
      theme_minimal() +
      theme(legend.position = "none")
    
    plot_list[[sel_protocol]] <- p
  }
  
  # Create a dummy plot to extract the legend
  legend_plot <- ggplot(metadata, aes(x = day, y = 1, color = cell_line)) +
    geom_point() +
    scale_color_manual(values = cell_line_palette, name = "Cell Line") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_text(face = "bold"))
  shared_legend <- get_legend(legend_plot)
  
  grid.arrange(
    shared_legend,
    arrangeGrob(grobs = plot_list, ncol = 3, nrow = 2),
    ncol = 1,
    heights = c(1, 9),
    top = sel_cell_type
  )
}

dev.off()

## ------------------------------
## Day 13 top-10 gene heatmaps
## ------------------------------

top10_error <- read_csv("braino_day013_top10_mean_error_by_protocol_cellline.csv",
                        show_col_types = FALSE)

top10_genes <- unique(top10_error$gene)

top10_protocol_order <- top10_error %>%
  group_by(protocol) %>%
  summarise(protocol_hits = n_distinct(gene), .groups = "drop") %>%
  arrange(desc(protocol_hits), protocol) %>%
  pull(protocol)

top10_cell_line_order <- top10_error %>%
  group_by(protocol, cell_line) %>%
  summarise(cell_line_hits = n_distinct(gene), .groups = "drop") %>%
  arrange(factor(protocol, levels = top10_protocol_order),
          desc(cell_line_hits), cell_line)

top10_column_levels <- top10_cell_line_order %>%
  mutate(protocol_cell_line = paste(protocol, cell_line, sep = "__")) %>%
  pull(protocol_cell_line)

top10_gene_order <- top10_error %>%
  count(gene, name = "gene_hits") %>%
  arrange(desc(gene_hits), gene) %>%
  pull(gene)

top10_grid <- expand_grid(
  gene = top10_genes,
  protocol_cell_line = top10_column_levels
) %>%
  separate(protocol_cell_line, into = c("protocol", "cell_line"),
           sep = "__", remove = FALSE) %>%
  left_join(
    top10_error %>%
      mutate(protocol_cell_line = paste(protocol, cell_line, sep = "__"),
             present = 1L) %>%
      select(protocol_cell_line, gene, present, mean_error),
    by = c("protocol_cell_line", "gene")
  ) %>%
  mutate(
    present = replace_na(present, 0L),
    protocol = factor(protocol, levels = top10_protocol_order),
    protocol_cell_line = factor(protocol_cell_line, levels = top10_column_levels),
    cell_line = factor(cell_line, levels = unique(top10_cell_line_order$cell_line)),
    gene = factor(gene, levels = rev(top10_gene_order))
  )

heatmap_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_top10_presence <- ggplot(top10_grid,
                           aes(x = protocol_cell_line, y = gene,
                               fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(. ~ protocol, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub(".*__", "", x), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(
    values = c("0" = "#f0f0f0", "1" = "#2b8cbe"),
    breaks = c("0", "1"),
    labels = c("Absent", "Present"),
    name = "Top-10"
  ) +
  labs(title = "Day 13 top-10 gene presence") +
  heatmap_theme

ggsave("braino_day013_top10_gene_presence_heatmap.pdf",
       p_top10_presence, width = 11, height = 8)

p_top10_mean_error <- ggplot(top10_grid,
                             aes(x = protocol_cell_line, y = gene,
                                 fill = mean_error)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(. ~ protocol, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub(".*__", "", x), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient(
    low = "#f7fbff", high = "#08306b",
    na.value = "#d9d9d9",
    name = "Mean error"
  ) +
  labs(title = "Day 13 top-10 gene mean error") +
  heatmap_theme

ggsave("braino_day013_top10_gene_mean_error_heatmap.pdf",
       p_top10_mean_error, width = 11, height = 8)

err_day13 <- read_csv("braino_rectangle_err.csv",
                      col_select = c(1, all_of(top10_genes)),
                      show_col_types = FALSE)

err_ids <- err_day13[[1]]
err_split_info <- strsplit(err_ids, "-")
err_metadata <- data.frame(
  replicate = paste0("R", sapply(err_split_info, `[`, 1)),
  cell_line = sapply(err_split_info, `[`, 2),
  protocol = sapply(err_split_info, `[`, 3),
  day = as.numeric(gsub("day", "", sapply(err_split_info, `[`, 4)))
)

err_all_mean <- cbind(err_metadata, err_day13[, top10_genes]) %>%
  as.data.frame() %>%
  filter(protocol != "pluripotent", day == 13) %>%
  pivot_longer(cols = all_of(top10_genes),
               names_to = "gene",
               values_to = "mean_error") %>%
  group_by(protocol, cell_line, gene) %>%
  summarise(mean_error = mean(mean_error, na.rm = TRUE), .groups = "drop") %>%
  mutate(protocol_cell_line = paste(protocol, cell_line, sep = "__"))

all_error_grid <- top10_grid %>%
  select(protocol, cell_line, protocol_cell_line, gene) %>%
  left_join(err_all_mean,
            by = c("protocol", "cell_line", "protocol_cell_line", "gene")) %>%
  mutate(
    protocol = factor(protocol, levels = top10_protocol_order),
    protocol_cell_line = factor(protocol_cell_line, levels = top10_column_levels),
    cell_line = factor(cell_line, levels = unique(top10_cell_line_order$cell_line)),
    gene = factor(gene, levels = rev(top10_gene_order))
  )

p_all_mean_error <- ggplot(all_error_grid,
                           aes(x = protocol_cell_line, y = gene,
                               fill = mean_error)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(. ~ protocol, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub(".*__", "", x), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0,
    name = "Mean error"
  ) +
  labs(title = "Day 13 mean error for top-10 gene union") +
  heatmap_theme

ggsave("braino_day013_top10_gene_all_mean_error_heatmap.pdf",
       p_all_mean_error, width = 11, height = 8)

# By protocol
pdf("rectangle_by_protocol.pdf", width = 16, height = 13)

for (sel_protocol in protocols) {
  
  plot_list <- list()
  
  for (sel_cell_type in cell_types) {
    
    y_values <- data[, sel_cell_type]
    y_min <- min(y_values, na.rm = TRUE)
    y_max <- max(y_values, na.rm = TRUE)
    
    plot_data <- metadata %>%
      mutate(y_value = data[rownames(metadata), sel_cell_type]) %>%
      filter(protocol == sel_protocol)
    
    p <- ggplot(plot_data, aes(x = day, y = y_value, color = cell_line, 
                               group = interaction(cell_line, replicate))) +
      geom_line() +
      geom_point() +
      scale_color_manual(values = cell_line_palette) +
      labs(x = "Day", y = "Proportion", title = sel_cell_type) +
      ylim(y_min, y_max) +
      theme_minimal() +
      theme(legend.position = "none")
    
    plot_list[[sel_cell_type]] <- p
  }
  
  # Create a dummy plot to extract the legend
  legend_plot <- ggplot(metadata, aes(x = day, y = 1, color = cell_line)) +
    geom_point() +
    scale_color_manual(values = cell_line_palette, name = "Cell Line") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_text(face = "bold"))
  shared_legend <- get_legend(legend_plot)
  
  grid.arrange(
    shared_legend,
    arrangeGrob(grobs = plot_list, ncol = 4, nrow = 4),
    ncol = 1,
    heights = c(1, 12),
    top = sel_protocol
  )
}

dev.off()

## ------------------------------
## Plots: mean vs. time
## ------------------------------

# By cell type
pdf("rectangle_by_celltype_mean.pdf", width = 14, height = 9)

for (sel_cell_type in cell_types) {
  
  summary_data <- metadata %>%
    mutate(y_value = data[rownames(metadata), sel_cell_type]) %>%
    group_by(protocol, cell_line, day) %>%
    summarise(
      mean_value = mean(y_value, na.rm = TRUE),
      min_value = min(y_value, na.rm = TRUE),
      max_value = max(y_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  y_min <- min(summary_data$min_value, na.rm = TRUE)
  y_max <- max(summary_data$max_value, na.rm = TRUE)
  
  plot_list <- list()
  
  for (sel_protocol in protocols) {
    
    plot_data <- summary_data %>%
      filter(protocol == sel_protocol)
    
    p <- ggplot(plot_data, aes(x = day, y = mean_value, color = cell_line, 
                               group = cell_line)) +
      geom_line() +
      geom_point() +
      geom_errorbar(aes(ymin = min_value, ymax = max_value), width = 0.5) +
      scale_color_manual(values = cell_line_palette) +
      labs(x = "Day", y = "Proportion", title = sel_protocol) +
      ylim(y_min, y_max) +
      theme_minimal() +
      theme(legend.position = "none")
    
    plot_list[[sel_protocol]] <- p
  }
  
  # Create a dummy plot to extract the legend
  legend_plot <- ggplot(metadata, aes(x = day, y = 1, color = cell_line)) +
    geom_point() +
    scale_color_manual(values = cell_line_palette, name = "Cell Line") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_text(face = "bold"))
  shared_legend <- get_legend(legend_plot)
  
  grid.arrange(
    shared_legend,
    arrangeGrob(grobs = plot_list, ncol = 3, nrow = 2),
    ncol = 1,
    heights = c(1, 9),
    top = sel_cell_type
  )
}

dev.off()

# By protocol
pdf("rectangle_by_protocol_mean.pdf", width = 16, height = 13)

for (sel_protocol in protocols) {
  
  plot_list <- list()
  
  for (sel_cell_type in cell_types) {
    
    summary_data <- metadata %>%
      mutate(y_value = data[rownames(metadata), sel_cell_type]) %>%
      group_by(protocol, cell_line, day) %>%
      summarise(
        mean_value = mean(y_value, na.rm = TRUE),
        min_value = min(y_value, na.rm = TRUE),
        max_value = max(y_value, na.rm = TRUE),
        .groups = "drop"
      )
    
    y_min <- min(summary_data$min_value, na.rm = TRUE)
    y_max <- max(summary_data$max_value, na.rm = TRUE)
    
    plot_data <- summary_data %>%
      filter(protocol == sel_protocol)
    
    p <- ggplot(plot_data, aes(x = day, y = mean_value, color = cell_line, 
                               group = cell_line)) +
      geom_line() +
      geom_point() +
      geom_errorbar(aes(ymin = min_value, ymax = max_value), width = 0.5) +
      scale_color_manual(values = cell_line_palette) +
      labs(x = "Day", y = "Proportion", title = sel_cell_type) +
      ylim(y_min, y_max) +
      theme_minimal() +
      theme(legend.position = "none")
    
    plot_list[[sel_cell_type]] <- p
  }
  
  # Create a dummy plot to extract the legend
  legend_plot <- ggplot(metadata, aes(x = day, y = 1, color = cell_line)) +
    geom_point() +
    scale_color_manual(values = cell_line_palette, name = "Cell Line") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_text(face = "bold"))
  shared_legend <- get_legend(legend_plot)
  
  grid.arrange(
    shared_legend,
    arrangeGrob(grobs = plot_list, ncol = 4, nrow = 4),
    ncol = 1,
    heights = c(1, 12),
    top = sel_protocol
  )
}

dev.off()

## ------------------------------
## Plots: stacked barplots
## ------------------------------

# Combine data
plot_data <- cbind(metadata, data[rownames(metadata), ]) %>%
  as.data.frame() %>%
  pivot_longer(cols = all_of(cell_types),
               names_to = "cell_type",
               values_to = "proportion")

# Set cell_type factor levels to ensure consistent ordering
plot_data <- plot_data %>%
  mutate(cell_type = factor(cell_type, levels = cell_types))

# Create sample label and order
plot_data <- plot_data %>%
  arrange(day, replicate) %>%
  mutate(sample_label = paste0("D", day, "_", replicate),
         sample_label = factor(sample_label, levels = unique(sample_label)))

# Add spacing between day groups
plot_data <- plot_data %>%
  mutate(day_group = as.numeric(factor(day)),
         x_pos = as.numeric(sample_label) + (day_group - 1) * 0.5)

# Calculate day label positions
day_labels <- plot_data %>%
  group_by(protocol, cell_line, day) %>%
  summarise(x_center = mean(x_pos), .groups = "drop")

ggplot(plot_data, aes(x = x_pos, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(data = day_labels, 
            aes(x = x_center, y = -0.05, label = paste0("D", day), fill = NULL),
            size = 2.5, vjust = 1) +
  facet_grid(cell_line ~ protocol, scales = "free_x", space = "free_x") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(limits = c(-0.1, 1), oob = scales::oob_keep) +
  scale_fill_manual(values = cell_type_palette) +
  labs(x = NULL, y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    plot.margin = margin(5, 5, 10, 5)
  ) +
  guides(fill = guide_legend(nrow = 2))

ggsave("rectangle_stacked_barplots.pdf", width = 16, height = 14)

# Stacked barplot focused on dorsal protocol for H9 and 176
sel_protocol <- "dorsal"
sel_cell_lines <- c("H9", "176")

plot_data_sel <- cbind(metadata, data[rownames(metadata), ]) %>%
  as.data.frame() %>%
  filter(protocol == sel_protocol, cell_line %in% sel_cell_lines) %>%
  pivot_longer(cols = all_of(cell_types),
               names_to = "cell_type",
               values_to = "proportion") %>%
  mutate(cell_type = factor(cell_type, levels = cell_types)) %>%
  group_by(cell_line) %>%
  arrange(day, replicate, .by_group = TRUE) %>%
  mutate(sample_label = paste0("D", day, "_", replicate),
         sample_order = as.numeric(factor(sample_label, levels = unique(sample_label))),
         day_group = as.numeric(factor(day)),
         x_pos = sample_order + (day_group - 1) * 0.5) %>%
  ungroup()

day_labels_sel <- plot_data_sel %>%
  group_by(cell_line, day) %>%
  summarise(x_center = mean(x_pos), .groups = "drop")

p_sel <- ggplot(plot_data_sel, aes(x = x_pos, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(data = day_labels_sel,
            aes(x = x_center, y = -0.05, label = paste0("D", day), fill = NULL),
            size = 3, vjust = 1) +
  facet_grid(cell_line ~ ., scales = "free_x", space = "free_x") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(limits = c(-0.1, 1), oob = scales::oob_keep) +
  scale_fill_manual(values = cell_type_palette) +
  labs(x = NULL, y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    plot.margin = margin(5, 5, 10, 5)
  ) +
  guides(fill = guide_legend(nrow = 2))

ggsave("rectangle_stacked_barplots_dorsal_H9_176.pdf", p_sel, width = 12, height = 8)

# Focus on selected cell lines

# Dorsal protocol only, H9 and 176 cell lines
pdf("rectangle_dorsal_H9_176.pdf", width = 16, height = 13)

sel_protocol <- "dorsal"
sel_cell_lines <- c("H9", "176")

plot_list <- list()

for (sel_cell_type in cell_types) {
  
  plot_data <- metadata %>%
    mutate(y_value = data[rownames(metadata), sel_cell_type]) %>%
    filter(protocol == sel_protocol, cell_line %in% sel_cell_lines)
  
  p <- ggplot(plot_data, aes(x = day, y = y_value, color = cell_line, 
                             group = interaction(cell_line, replicate))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = cell_line_palette[sel_cell_lines]) +
    labs(x = "Day", y = "Proportion", title = sel_cell_type) +
    ylim(0, 1) +
    theme_minimal() +
    theme(legend.position = "none")
  
  plot_list[[sel_cell_type]] <- p
}

# Create a dummy plot to extract the legend
legend_plot <- ggplot(metadata %>% filter(cell_line %in% sel_cell_lines), 
                      aes(x = day, y = 1, color = cell_line)) +
  geom_point() +
  scale_color_manual(values = cell_line_palette[sel_cell_lines], name = "Cell Line") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"))
shared_legend <- get_legend(legend_plot)

grid.arrange(
  shared_legend,
  arrangeGrob(grobs = plot_list, ncol = 4, nrow = 4),
  ncol = 1,
  heights = c(1, 12),
  top = sel_protocol
)

dev.off()

# Dorsal protocol only, H9 and 176 cell lines - selected cell types

sel_protocol <- "dorsal"
sel_cell_lines <- c("H9", "176")
sel_cell_types_requested <- c(
  "dorsal forebrain",
  "dorsal forebrain prog.",
  "hindbrain prog.",
  "Unknown"
)
sel_cell_type_labels <- c(
  "Dorsal forebrain",
  "Dorsal forebrain progenitors",
  "Hindbrain progenitors",
  "Unknown"
)

normalize_cell_type <- function(x) {
  tolower(trimws(gsub("\\.", " ", trimws(x))))
}

matched_cell_types <- sapply(
  normalize_cell_type(sel_cell_types_requested),
  function(sel_cell_type) {
    match_idx <- match(sel_cell_type, normalize_cell_type(cell_types))
    if (is.na(match_idx)) NA_character_ else cell_types[match_idx]
  },
  USE.NAMES = FALSE
)

if (all(is.na(matched_cell_types))) {
  stop("None of the requested selected cell types were found in the data")
}

missing_cell_types <- sel_cell_types_requested[is.na(matched_cell_types)]
if (length(missing_cell_types) > 0) {
  message("Requested cell types not found: ",
          paste(missing_cell_types, collapse = ", "))
}

sel_cell_type_labels <- sel_cell_type_labels[!is.na(matched_cell_types)]
matched_cell_types <- matched_cell_types[!is.na(matched_cell_types)]
names(sel_cell_type_labels) <- matched_cell_types

plot_list <- list()

for (sel_cell_type in matched_cell_types) {

  plot_data <- metadata %>%
    mutate(y_value = data[rownames(metadata), sel_cell_type]) %>%
    filter(protocol == sel_protocol, cell_line %in% sel_cell_lines)

  p <- ggplot(plot_data, aes(x = day, y = y_value, color = cell_line,
                             group = interaction(cell_line, replicate))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = cell_line_palette[sel_cell_lines]) +
    scale_y_continuous(limits = c(0, 0.8),
                       labels = function(x) x * 100) +
    labs(x = "Day", y = "Unknown content (%)",
         title = sel_cell_type_labels[[sel_cell_type]]) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.margin = margin(5.5, 18, 5.5, 5.5)
    )

  plot_list[[sel_cell_type]] <- p
}

legend_plot <- ggplot(metadata %>% filter(cell_line %in% sel_cell_lines),
                      aes(x = day, y = 1, color = cell_line)) +
  geom_point() +
  scale_color_manual(values = cell_line_palette[sel_cell_lines], name = "Cell Line") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"))
shared_legend <- get_legend(legend_plot)

selected_celltypes_grob <- arrangeGrob(
  shared_legend,
  arrangeGrob(grobs = plot_list,
              ncol = 4,
              nrow = ceiling(length(plot_list) / 4)),
  ncol = 1,
  heights = c(1, 10)
)

ggsave("rectangle_dorsal_H9_176_selected_celltypes.pdf",
       selected_celltypes_grob, width = 10, height = 3)

## ------------------------------
## Manuscript min-max values
## ------------------------------

format_pct <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f"), 100 * x)
}

format_pct_range <- function(x, digits = 1) {
  range_values <- range(x, na.rm = TRUE)
  paste0(format_pct(range_values[1], digits), "-", format_pct(range_values[2], digits), "%")
}

print_range <- function(label, values) {
  cat(sprintf("%-75s %s\n", label, format_pct_range(values)))
}

find_cell_type <- function(requested_cell_type) {
  match_idx <- match(normalize_cell_type(requested_cell_type),
                     normalize_cell_type(cell_types))
  if (is.na(match_idx)) {
    stop(sprintf("Column '%s' not found in braino_rectangle_estimations.csv",
                 requested_cell_type))
  }
  cell_types[match_idx]
}

dorsal_forebrain_col <- find_cell_type("dorsal forebrain")
dorsal_forebrain_prog_col <- find_cell_type("dorsal forebrain prog.")
hindbrain_col <- find_cell_type("hindbrain")
hindbrain_prog_col <- find_cell_type("hindbrain prog.")
unknown_col <- find_cell_type("Unknown")

manuscript_data <- cbind(metadata, data[rownames(metadata), ]) %>%
  as.data.frame() %>%
  filter(protocol == "dorsal", cell_line %in% c("H9", "176")) %>%
  mutate(
    dorsal_forebrain_total = .data[[dorsal_forebrain_col]] +
      .data[[dorsal_forebrain_prog_col]],
    hindbrain_total = .data[[hindbrain_col]] + .data[[hindbrain_prog_col]]
  )

first_day <- min(manuscript_data$day, na.rm = TRUE)
final_day <- max(manuscript_data$day, na.rm = TRUE)

cat("\nManuscript MIN-MAX values (dorsal protocol; replicate-level Rectangle fractions)\n")
cat(sprintf("First time point: day %s\n", first_day))
cat(sprintf("Final time point: day %s\n", final_day))

print_range(
  sprintf("176 Unknown, day %s", first_day),
  manuscript_data %>%
    filter(cell_line == "176", day == first_day) %>%
    pull(all_of(unknown_col))
)
print_range(
  sprintf("H9 Unknown, day %s", first_day),
  manuscript_data %>%
    filter(cell_line == "H9", day == first_day) %>%
    pull(all_of(unknown_col))
)
print_range(
  sprintf("H9 dorsal forebrain + dorsal forebrain prog., day %s", final_day),
  manuscript_data %>%
    filter(cell_line == "H9", day == final_day) %>%
    pull(dorsal_forebrain_total)
)
print_range(
  "H9 dorsal forebrain + dorsal forebrain prog., day 80",
  manuscript_data %>%
    filter(cell_line == "H9", day == 80) %>%
    pull(dorsal_forebrain_total)
)
print_range(
  sprintf("176 hindbrain + hindbrain prog., day %s", final_day),
  manuscript_data %>%
    filter(cell_line == "176", day == final_day) %>%
    pull(hindbrain_total)
)
cat("\n")

# Plot with bundles

# By protocol with shaded ribbon for SD
pdf("rectangle_by_protocol_mean_ribbon.pdf", width = 16, height = 13)

for (sel_protocol in protocols) {
  
  plot_list <- list()
  
  for (sel_cell_type in cell_types) {
    
    summary_data <- metadata %>%
      mutate(y_value = data[rownames(metadata), sel_cell_type]) %>%
      group_by(protocol, cell_line, day) %>%
      summarise(
        mean_value = mean(y_value, na.rm = TRUE),
        sd_value = sd(y_value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(sd_value = replace_na(sd_value, 0))
    
    y_min <- min(summary_data$mean_value - summary_data$sd_value, na.rm = TRUE)
    y_max <- max(summary_data$mean_value + summary_data$sd_value, na.rm = TRUE)
    
    plot_data <- summary_data %>%
      filter(protocol == sel_protocol)
    
    p <- ggplot(plot_data, aes(x = day, y = mean_value, color = cell_line, 
                               fill = cell_line, group = cell_line)) +
      geom_ribbon(aes(ymin = mean_value - sd_value, 
                      ymax = mean_value + sd_value), 
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      scale_color_manual(values = cell_line_palette) +
      scale_fill_manual(values = cell_line_palette) +
      labs(x = "Day", y = "Proportion", title = sel_cell_type) +
      ylim(y_min, y_max) +
      theme_minimal() +
      theme(legend.position = "none")
    
    plot_list[[sel_cell_type]] <- p
  }
  
  # Create a dummy plot to extract the legend
  legend_plot <- ggplot(metadata, aes(x = day, y = 1, color = cell_line, fill = cell_line)) +
    geom_point() +
    scale_color_manual(values = cell_line_palette, name = "Cell Line") +
    scale_fill_manual(values = cell_line_palette, name = "Cell Line") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_text(face = "bold"))
  shared_legend <- get_legend(legend_plot)
  
  grid.arrange(
    shared_legend,
    arrangeGrob(grobs = plot_list, ncol = 4, nrow = 4),
    ncol = 1,
    heights = c(1, 12),
    top = sel_protocol
  )
}

dev.off()

# Unknown only: all protocols side by side with shaded ribbon for SD
if (!"Unknown" %in% cell_types) {
  stop("Column 'Unknown' not found in braino_rectangle_estimations.csv")
}

unknown_summary <- metadata %>%
  mutate(y_value = data[rownames(metadata), "Unknown"]) %>%
  group_by(protocol, cell_line, day) %>%
  summarise(
    mean_value = mean(y_value, na.rm = TRUE),
    sd_value = sd(y_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd_value = replace_na(sd_value, 0),
    ribbon_min = pmax(mean_value - sd_value, 0),
    ribbon_max = mean_value + sd_value
  )

unknown_y_max <- max(unknown_summary$ribbon_max, na.rm = TRUE)
unknown_y_limit <- min(1, unknown_y_max * 1.08)

unknown_plot_list <- list()

for (sel_protocol in protocols) {

  plot_data <- unknown_summary %>%
    filter(protocol == sel_protocol)

  protocol_label <- tools::toTitleCase(sel_protocol)

  p <- ggplot(plot_data, aes(x = day, y = mean_value, color = cell_line,
                             fill = cell_line, group = cell_line)) +
    geom_ribbon(aes(ymin = ribbon_min, ymax = ribbon_max),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.7, stroke = 0) +
    scale_color_manual(values = cell_line_palette) +
    scale_fill_manual(values = cell_line_palette) +
    scale_y_continuous(labels = function(x) x * 100,
                       expand = expansion(mult = c(0.02, 0.06))) +
    scale_x_continuous(breaks = sort(unique(plot_data$day))) +
    coord_cartesian(ylim = c(0, unknown_y_limit)) +
    labs(x = "Day", y = paste0(protocol_label, " unknown content (%)")) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      plot.margin = margin(7, 8, 7, 8)
    )

  unknown_plot_list[[sel_protocol]] <- p
}

legend_plot <- ggplot(unknown_summary,
                      aes(x = day, y = mean_value, color = cell_line,
                          fill = cell_line, group = cell_line)) +
  geom_ribbon(aes(ymin = ribbon_min, ymax = ribbon_max),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.7, stroke = 0) +
  scale_color_manual(values = cell_line_palette, name = "Cell line") +
  scale_fill_manual(values = cell_line_palette, name = "Cell line") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1))
shared_legend <- get_legend(legend_plot)

pdf("rectangle_unknown_by_protocol_mean_ribbon.pdf", width = 10, height = 3)

grid.arrange(
  arrangeGrob(grobs = unknown_plot_list, ncol = length(unknown_plot_list)),
  shared_legend,
  ncol = 1,
  heights = c(8, 1.15)
)

dev.off()
