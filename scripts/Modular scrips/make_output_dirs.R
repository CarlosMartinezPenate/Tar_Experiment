# make_output_dirs.R
make_output_dirs <- function(root_dir, control_group, comparison_group, filter_type) {
  comp_name <- paste0(comparison_group, "_vs_", control_group, "_", filter_type)
  base_dir <- file.path(root_dir, comp_name)
  
  # Define subdirectories
  dirs_to_create <- c(
    base_dir,
    file.path(base_dir, "results"),
    file.path(base_dir, "plots"),
    file.path(base_dir, "heatmaps"),
    file.path(base_dir, "volcanoes"),
    file.path(base_dir, "deseq2_only")
  )
  
  # Create each directory if it doesn't exist
  for (dir_path in dirs_to_create) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      message("📁 Created: ", dir_path)
    }
  }
  
  return(base_dir)  # useful if you want to refer to it later
}