run_super_comparisons <- function(physeq,
                                  control_groups,
                                  comparison_groups,
                                  filter_types = c("sterivex", "11_um"),
                                  thresholds,
                                  methods = c("DESeq2", "ALDEx2", "LEfSe"),
                                  tax_level = "Genus",
                                  transform_type = "log",
                                  output_root = "./super_comparisons") {
  group_var <- "Station_treatment"
  
  for (control in control_groups) {
    for (comparison in comparison_groups) {
      for (filt in filter_types) {
        
        comp_id <- paste0(comparison, "_vs_", control, "_", filt)
        comp_dir <- file.path(output_root, comp_id)
        
        make_output_dirs(
          comp_dir,
          control_group = control,
          comparison_group = comparison,
          filter_type = filt
        )
        
        message("🚀 Running comparison: ", comp_id)
        
        phy_sub <- safe_subset_phyloseq(
          physeq,
          group_var = group_var,
          control_group = control,
          contrast_group = comparison,
          filter_type = filt,
          min_total_reads = 120
        )
        
        if (is.null(phy_sub) || !inherits(phy_sub, "phyloseq")) {
          message(sprintf("⛔ Skipping comparison %s: invalid or NULL phyloseq object", comp_id))
          next
        }
        
        cat("📦 Subset class:\n")
        print(class(phy_sub))
        cat("📊 Taxa × Samples:\n")
        print(dim(as(otu_table(phy_sub), "matrix")))
        
        da_df <- run_differential_abundance_all(
          physeq_obj = phy_sub,
          group_var = group_var,
          contrast_group = comparison,
          control_group = control,
          filter_type = filt,
          method_list = methods,
          thresholds = thresholds,
          output_dir = comp_dir,
          plot = FALSE
        )
        
        flagged <- flag_significance(da_df, thresholds)
        flagged_path <- file.path(comp_dir, sprintf("results_asvs_%s_vs_%s_%s.csv", comparison, control, filt))
        write.csv(flagged, flagged_path, row.names = FALSE)
        message("✅ Flagged results saved to: ", flagged_path)
        
        # Volcano Plot
        volcano_df <- prepare_volcano_data(
          flagged_df = flagged,
          physeq = phy_sub,
          group_var = group_var,
          control_group = control,
          contrast_group = comparison,
          tax_level = tax_level
        )
        cat("📋 Columns in flagged_df:\n")
        print(names(flagged))
        volcano_path <- file.path(comp_dir, sprintf("volcano_%s_vs_%s_%s.png", comparison, control, filt))
        message("📊 Saving volcano to: ", volcano_path)
        
        plot_combined_volcano(
          volcano_df,
          out_path = volcano_path,
          title = paste("Combined Volcano:", comparison, "vs", control, "(", filt, ")")
        )
        
        # 🌡️ Heatmap (robust-only)
        heatmap_path <- file.path(comp_dir, sprintf("heatmap_%s_vs_%s_%s.png", comparison, control, filt))
        message("🌡️ Saving heatmap to: ", heatmap_path)
        message("📊 Generating heatmap with ", nrow(flagged), " entries in flagged_df")
        tryCatch({
          result <- generate_heatmap_from_flagged(
            flagged_df = flagged,
            physeq = phy_sub,
            control_group = control,
            contrast_group = comparison,
            group_var = group_var,
            filter_type = filt,
            transform_type = transform_type,
            robust_only = TRUE,
            save_path = heatmap_path
          )
          if (is.null(result)) {
            message("⚠️ Heatmap generation returned NULL. Possibly due to zero significant taxa.")
          } else {
            message("✅ Heatmap successfully generated.")
          }
        }, error = function(e) {
          warning(sprintf("❌ Error generating heatmap for %s: %s", comp_id, e$message))
        })
      }
    }
  }
  
  # 📈 Run once after all comparisons are done
  generate_deseq2_only_volcanoes_pretty(
    results_dir = output_root,
    physeq = physeq,
    group_var = group_var,
    tax_level = tax_level
  )
  
  message("✅ All super comparisons completed.")
}