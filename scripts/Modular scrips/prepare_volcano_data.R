prepare_volcano_data <- function(flagged_df,
                                 physeq,
                                 group_var,
                                 control_group,
                                 contrast_group,
                                 tax_level = "Genus") {
  cat("🔍 Starting prepare_volcano_data()\n")
  
  # Validate inputs
  if (!"TaxaID" %in% colnames(flagged_df)) {
    stop("❌ flagged_df must contain 'TaxaID' column")
  }
  
  # Handle missing n_methods
  if (!"n_methods" %in% names(flagged_df)) {
    warning("⚠️ 'n_methods' not found in flagged_df. Defaulting to NA.")
    flagged_df$n_methods <- NA
  }
  
  # Extract taxonomy table
  tax_tab <- tax_table(physeq)
  
  if (inherits(tax_tab, "matrix") || inherits(tax_tab, "array")) {
    tax_tab <- as.data.frame(tax_tab)
  } else {
    tax_tab <- as.data.frame(as.matrix(tax_tab))
  }
  
  cat("📦 Converting tax_table to matrix and then data.frame...\n")
  cat("📏 Tax table dimensions:", dim(tax_tab), "\n")
  
  tax_tab$TaxaID <- rownames(tax_tab)
  
  # Reduce tax table to desired level
  if (!(tax_level %in% colnames(tax_tab))) {
    stop(paste("❌ tax_level", tax_level, "not found in taxonomy table"))
  }
  
  tax_sub <- tax_tab[, c("TaxaID", tax_level)]
  colnames(tax_sub)[2] <- "TaxLabel"
  
  cat("✅ Tax data preview:\n")
  print(head(tax_sub, 3))
  
  # Merge taxonomy info
  volcano_df <- merge(flagged_df, tax_sub, by = "TaxaID", all.x = TRUE)
  
  cat("🔗 Merged taxonomy. Dimensions:\n")
  print(dim(volcano_df))
  
  # Set defaults for plotting aesthetics
  if (!"score" %in% names(volcano_df)) {
    volcano_df$score <- volcano_df$log2FoldChange
    cat("⚠️ No 'score' column found. Using log2FoldChange instead.\n")
  }
  
  # Set significance factor
  volcano_df$significance_label <- ifelse(volcano_df$is_significant, "Significant", "Not Significant")
  
  # Replace NA n_methods if needed
  if (any(is.na(volcano_df$n_methods))) {
    volcano_df$n_methods[is.na(volcano_df$n_methods)] <- 0
  }
  
  # Final preview
  cat("✅ Prepared volcano_df:\n")
  print(head(volcano_df, 3))
  
  cat("✅ Finished prepare_volcano_data()\n")
  
  return(volcano_df)
}