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
        
        message("\n🔧 Setting up output directories for: ", comp_id)
        make_output_dirs(
          comp_dir,
          control_group = control,
          comparison_group = comparison,
          filter_type = filt
        )
        
        message("🚀 Starting comparison: ", comp_id)
        
        message("🔍 Subsetting phyloseq...")
        phy_sub <- safe_subset_phyloseq(
          physeq,
          group_var = group_var,
          control_group = control,
          contrast_group = comparison,
          filter_type = filt,
          min_total_reads = 120
        )
        
        if (is.null(phy_sub) || !inherits(phy_sub, "phyloseq")) {
          message(sprintf("⛔ Skipping %s: NULL or invalid phyloseq object returned.", comp_id))
          next
        }
        
        message("📦 Phyloseq subset created.")
        cat("📦 Subset class:\n")
        print(class(phy_sub))
        cat("📊 Taxa × Samples:\n")
        print(dim(as(otu_table(phy_sub), "matrix")))
        
        message("📈 Running differential abundance methods: ", paste(methods, collapse = ", "))
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
        message("✅ Differential abundance completed.")
        
        message("🚩 Flagging significant taxa...")
        flagged <- flag_significance(da_df, thresholds)
        flagged_path <- file.path(comp_dir, sprintf("results_asvs_%s_vs_%s_%s.csv", comparison, control, filt))
        write.csv(flagged, flagged_path, row.names = FALSE)
        message("📝 Flagged taxa written to: ", flagged_path)
        
        message("🧬 Preparing volcano data...")
        volcano_df <- prepare_volcano_data(
          flagged_df = flagged,
          physeq = phy_sub,
          group_var = group_var,
          control_group = control,
          contrast_group = comparison,
          tax_level = tax_level
        )
        
        volcano_path <- file.path(comp_dir, sprintf("volcano_%s_vs_%s_%s.png", comparison, control, filt))
        message("🧨 Generating volcano plot...")
        plot_combined_volcano(
          volcano_df,
          out_path = volcano_path,
          title = paste("Combined Volcano:", comparison, "vs", control, "(", filt, ")")
        )
        message("✅ Volcano plot saved: ", volcano_path)
        
        if (exists("generate_deseq2_only_volcanoes_pretty")) {
          message("🎨 Generating DESeq2-only pretty volcano plots...")
          generate_deseq2_only_volcanoes_pretty(
            results_dir = output_root,
            physeq = physeq,
            group_var = group_var,
            tax_level = tax_level
          )
        }
        
        message("🌡️ Starting heatmap generation...")
        heatmap_path <- file.path(comp_dir, sprintf("heatmap_%s_vs_%s_%s.png", comparison, control, filt))
        
        tryCatch({
          generate_heatmap_from_flagged(
            flagged_df = flagged,
            physeq = phy_sub,
            control_group = control,
            contrast_group = comparison,
            group_var = group_var,
            filter_type = filt,
            transform_type = transform_type,
            save_path = heatmap_path
          )
          message("✅ Heatmap saved: ", heatmap_path)
        }, error = function(e) {
          warning(sprintf("❌ Heatmap generation failed for %s: %s", comp_id, e$message))
        })
        
        message("🎯 Comparison finished: ", comp_id)
      }
    }
  }
  
  message("\n🏁 All super comparisons complete.")
}