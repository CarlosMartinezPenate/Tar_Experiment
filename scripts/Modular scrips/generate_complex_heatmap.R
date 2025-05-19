generate_complex_heatmap <- function(flagged_df,
                                     physeq,
                                     control_group,
                                     contrast_group,
                                     group_var,
                                     filter_type,
                                     transform_type = "ra",
                                     tax_level = "Genus",
                                     robust_only = TRUE,
                                     top_n_taxa = 50,
                                     output_dir = "complex_heatmaps") {
  cat("🧼 Starting ComplexHeatmap heatmap generation...\n")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  comp_id <- paste0(contrast_group, "_vs_", control_group, "_", filter_type)
  save_path <- file.path(output_dir, sprintf("complex_heatmap_%s.png", comp_id))
  cat("🧾 Saving to directory:", output_dir, "\n")
  
  cat("📊 Flagged DF dimensions: ", dim(flagged_df), "\n")
  
  if (robust_only && "n_methods" %in% colnames(flagged_df)) {
    flagged_df <- flagged_df[flagged_df$n_methods >= 2 & flagged_df$is_significant == TRUE, ]
    cat("🔎 Filtered robust taxa (n_methods >=  2 ): ", nrow(flagged_df), " remaining\n")
  }
  
  sig_taxa <- unique(flagged_df$TaxaID)
  physeq_sig <- prune_taxa(sig_taxa, physeq)
  cat("📦 Filtered phyloseq → taxa: ", ntaxa(physeq_sig), ", samples: ", nsamples(physeq_sig), "\n")
  
  mat <- as(otu_table(physeq_sig), "matrix")
  if (!taxa_are_rows(physeq_sig)) mat <- t(mat)
  cat("✅ Extracted OTU matrix: ", dim(mat)[1], " taxa × ", dim(mat)[2], " samples\n")
  
  # Transform
  if (transform_type == "zscore") {
    mat <- t(scale(t(mat)))
  } else if (transform_type == "log2") {
    mat <- log2(mat + 1e-6)
  } else {
    mat <- sweep(mat, 2, colSums(mat), FUN = "/")
  }
  
  mat[is.na(mat)] <- 0
  cat("🔄 Transformed to", transform_type, "\n")
  
  # Filter top N most variable taxa
  if (!is.null(top_n_taxa) && nrow(mat) > top_n_taxa) {
    vars <- apply(mat, 1, var, na.rm = TRUE)
    top_taxa <- names(sort(vars, decreasing = TRUE))[1:top_n_taxa]
    mat <- mat[top_taxa, , drop = FALSE]
    cat("🔪 Selected top ", top_n_taxa, " variable taxa\n")
  }
  
  # Taxa labels
  tax_df <- as.data.frame(tax_table(physeq_sig))
  taxa_labels <- paste0(rownames(mat), " | ", tax_df[rownames(mat), tax_level])
  rownames(mat) <- taxa_labels
  cat("🧬 Constructed taxa labels with level: ", tax_level, "\n")
  
  # Metadata
  meta <- data.frame(sample_data(physeq_sig))
  meta$SampleID <- rownames(meta)
  cat("📋 Sample metadata dimensions: ", dim(meta)[1], "\n")
  
  group_colors <- structure(RColorBrewer::brewer.pal(8, "Dark2"), names = unique(meta[[group_var]]))
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = data.frame(Group = meta[[group_var]]),
    col = list(Group = group_colors),
    annotation_name_side = "left"
  )
  
  # Save plot
  tryCatch({
    if (file.exists(save_path)) unlink(save_path)
    png(save_path, width = 1000, height = 500 + nrow(mat) * 12, res = 150)
    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = transform_type,
      top_annotation = ha,
      cluster_columns = TRUE,
      cluster_rows = TRUE,
      show_column_names = TRUE,
      show_row_names = TRUE,
      column_names_gp = grid::gpar(fontsize = 10),
      row_names_gp = grid::gpar(fontsize = 8),
      heatmap_legend_param = list(title = transform_type),
      column_title = sprintf("ComplexHeatmap: %s vs %s (%s)", contrast_group, control_group, filter_type)
    )
    draw(ht)
    dev.off()
    cat("✅ ComplexHeatmap saved to:", save_path, "\n")
  }, error = function(e) {
    warning("❌ Failed to save ComplexHeatmap:", e$message)
  })
  
  invisible(NULL)
}