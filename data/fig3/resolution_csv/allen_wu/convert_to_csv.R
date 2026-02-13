# Find all deconvolution.rds files in replicate folders
rds_files <- list.files(
    path = ".",
    pattern = "deconvolution\\.rds$",
    recursive = TRUE,
    full.names = TRUE
)

# Filter to only include replicate folders
rds_files <- rds_files[grepl("replicate_", rds_files)]

for (rds_path in rds_files) {
    message(paste("Processing:", rds_path))

    # Read the RDS file
    data <- readRDS(rds_path)

    # Get the directory of the current RDS file
    current_dir <- dirname(rds_path)

    # Check if deconvolution table exists and write to CSV
    if (!is.null(data$deconvolution)) {
        write.csv(data$deconvolution,
            file = file.path(current_dir, "deconvolution.csv"),
            row.names = TRUE
        )
        message(paste("  - Saved deconvolution.csv"))
    }

    # Check if true_cell_fractions table exists and write to CSV
    if (!is.null(data$true_cell_fractions)) {
        write.csv(data$true_cell_fractions,
            file = file.path(current_dir, "true_cell_fractions.csv"),
            row.names = TRUE
        )
        message(paste("  - Saved true_cell_fractions.csv"))
    }
}

message("Conversion complete.")
