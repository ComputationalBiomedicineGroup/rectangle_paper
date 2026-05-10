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
pdf("rectangle_dorsal_H9_176_selected_celltypes.pdf", width = 15, height = 9)

sel_protocol <- "dorsal"
sel_cell_lines <- c("H9", "176")
sel_cell_types_requested <- c(
  "Unknown",
  "dorsal.forebrain",
  "dorsal.forebrain.prog.",
  "hindbrain.prog.",
  "muscle",
  "eye"
)

normalize_cell_type <- function(x) {
  tolower(gsub("\\.", " ", trimws(x)))
}

matched_cell_types <- cell_types[
  normalize_cell_type(cell_types) %in% normalize_cell_type(sel_cell_types_requested)
]

if (length(matched_cell_types) == 0) {
  stop("None of the requested selected cell types were found in the data")
}

missing_cell_types <- setdiff(
  normalize_cell_type(sel_cell_types_requested),
  normalize_cell_type(matched_cell_types)
)
if (length(missing_cell_types) > 0) {
  message("Requested cell types not found (normalized): ",
          paste(missing_cell_types, collapse = ", "))
}

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
    labs(x = "Day", y = "Proportion", title = sel_cell_type) +
    ylim(0, 1) +
    theme_minimal() +
    theme(
      legend.position = "none",
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

grid.arrange(
  shared_legend,
  arrangeGrob(grobs = plot_list,
              ncol = 3,
              nrow = ceiling(length(plot_list) / 3)),
  ncol = 1,
  heights = c(1, 10),
  top = paste0(sel_protocol, " - selected cell types")
)

dev.off()

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

# Unknown only: by protocol with shaded ribbon for SD
if (!"Unknown" %in% cell_types) {
  stop("Column 'Unknown' not found in braino_rectangle_estimations.csv")
}

pdf("rectangle_unknown_by_protocol_mean_ribbon.pdf", width = 12, height = 8)

for (sel_protocol in protocols) {

  summary_data <- metadata %>%
    mutate(y_value = data[rownames(metadata), "Unknown"]) %>%
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
    labs(x = "Day", y = "Proportion", title = "Unknown") +
    ylim(y_min, y_max) +
    theme_minimal() +
    theme(legend.position = "none")

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
    p,
    ncol = 1,
    heights = c(1, 8),
    top = sel_protocol
  )
}

dev.off()

