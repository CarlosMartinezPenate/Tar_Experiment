generate_heatmap_from_flagged <- function(flagged_df,
                                          physeq,
                                          control_group,
                                          contrast_group,
                                          group_var,
                                          filter_type,
                                          transform_type = "log",
                                          save_path = NULL) {
  message("Starting heatmap_from_flagged () 🧼 Subsetting phyloseq for heatmap...")
  
  # Debug print: dimensions and head
  cat("📋 flagged_df dimensions:", dim(flagged_df), "\n")
  print(head(flagged_df, 3))
  
  if (nrow(flagged_df) == 0) {
    warning("❌ No significant taxa in flagged_df — skipping heatmap.")
    return(invisible(NULL))
  }
  
  # Get vector of significant TaxaIDs
  sig_taxa <- flagged_df$TaxaID
  cat("🔍 Number of significant taxa:", length(sig_taxa), "\n")
  
  physeq_subset <- phyloseq::prune_taxa(sig_taxa, physeq)
  
  # Transform counts if requested
  if (transform_type == "log") {
    message("🔧 Applying log1p transformation to counts...")
    physeq_subset <- phyloseq::transform_sample_counts(physeq_subset, function(x) log1p(x))
  } else if (transform_type == "clr") {
    stop("CLR transformation not implemented yet in this function.")
  }
  
  # Create matrix
  mat <- as(phyloseq::otu_table(physeq_subset), "matrix")
  
  # Debug: show taxa names and sample names
  cat("🧬 Taxa names (first 5):", paste(rownames(mat)[1:5], collapse = ", "), "\n")
  cat("🧪 Sample names (first 5):", paste(colnames(mat)[1:5], collapse = ", "), "\n")
  
  # Reorder samples based on metadata
  sample_metadata <- phyloseq::sample_data(physeq_subset)[[group_var]]
  cat("📑 Sample group variable levels:", paste(levels(factor(sample_metadata)), collapse = ", "), "\n")
  
  sample_order <- order(sample_metadata)
  mat <- mat[, sample_order, drop = FALSE]
  
  # Diagnostics
  message("🧪 ", nrow(mat), " taxa | ", ncol(mat), " samples in heatmap phyloseq")
  message("📏 Heatmap matrix dimensions: ", paste(dim(mat), collapse = " × "))
  message("🔍 Any NA? ", anyNA(mat), " | All NA? ", all(is.na(mat)))
  
  if (nrow(mat) == 0 || ncol(mat) == 0 || all(is.na(mat))) {
    warning("❌ Heatmap matrix has zero dimensions or all NA — skipping.")
    return(invisible(NULL))
  }
  
  # Plot heatmap safely
  message("🌡️ Generating heatmap...")
  
  heatmap_expr <- function() {
    pheatmap::pheatmap(mat,
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       fontsize_row = 6,
                       fontsize_col = 8,
                       main = paste("Heatmap:", contrast_group, "vs", control_group, "(", filter_type, ")"))
  }
  
  tryCatch({
    if (!is.null(save_path)) {
      grDevices::png(filename = save_path, width = 1000, height = 800)
      heatmap_expr()
      grDevices::dev.off()
      message("✅ Heatmap saved to: ", save_path)
    } else {
      heatmap_expr()
    }
  }, error = function(e) {
    message("❌ Heatmap generation failed: ", conditionMessage(e))
  })
  
  invisible(NULL)
}