generate_all_combined_volcanoes <- function(results_dir,
                                            physeq,
                                            group_var,
                                            tax_level,
                                            out_subdir = "combined_volcano_plots",
                                            label_top_n = 20) {
  files <- list.files(results_dir, pattern = "^results_asvs_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop("❌ No result files found.")
  
  out_dir <- file.path(results_dir, out_subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (file in files) {
    fname <- basename(file)
    parts <- stringr::str_match(fname, "^results_asvs_(.*?)_vs_(.*?)_(.*)\\.csv$")
    if (any(is.na(parts))) next
    
    contrast_group <- parts[2]
    control_group <- parts[3]
    filter_type <- parts[4]
    
    cat(sprintf("📊 Plotting: %s vs %s (%s)\n", contrast_group, control_group, filter_type))
    
    flagged_df <- read.csv(file)
    
    volcano_ready <- prepare_volcano_data(
      flagged_df = flagged_df,
      physeq = physeq,
      group_var = group_var,
      control_group = control_group,
      contrast_group = contrast_group,
      tax_level = tax_level
    )
    
    # Save table + volcano
    write.csv(volcano_ready,
              file.path(out_dir, sprintf("volcano_data_%s_vs_%s_%s.csv", contrast_group, control_group, filter_type)),
              row.names = FALSE)
    
    plot_combined_volcano(
      volcano_df = volcano_ready,
      title = sprintf("Combined Volcano: %s vs %s (%s)", contrast_group, control_group, filter_type),
      out_path = file.path(out_dir, sprintf("combined_volcano_%s_vs_%s_%s.png", contrast_group, control_group, filter_type)),
      label_top_n = label_top_n
    )
  }
}