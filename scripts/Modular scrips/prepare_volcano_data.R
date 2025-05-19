<<<<<<< HEAD

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

=======
>>>>>>> parent of 05ae7b2 (changes to make it work)
prepare_volcano_data <- function(flagged_df, physeq, group_var, control_group, contrast_group, tax_level = "Genus") {
  relab_df <- otu_table(physeq)
  if (!taxa_are_rows(physeq)) relab_df <- t(relab_df)
  relab_df <- sweep(relab_df, 2, colSums(relab_df), FUN = "/")
  relab_df <- as.data.frame(relab_df)

  
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
  
<<<<<<< HEAD

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

n_methods_df <- flagged_df %>%
=======
  n_methods_df <- flagged_df %>%
>>>>>>> parent of 05ae7b2 (changes to make it work)
    filter(is_significant) %>%
    group_by(TaxaID) %>%
    summarise(n_methods = n_distinct(method), .groups = "drop")
  
  tax_df <- as.data.frame(tax_table(physeq)) %>%
    rownames_to_column("TaxaID") %>%
    mutate(tax_group = .data[[tax_level]]) %>%
    select(TaxaID, tax_group)
  
  joined_df <- flagged_df %>%
    left_join(relab_df_combined, by = "TaxaID") %>%
    left_join(n_methods_df, by = "TaxaID") %>%
    left_join(tax_df, by = "TaxaID")
  
  joined_df %>%
    mutate(
      adjusted_p_value = as.numeric(as.character(adjusted_p_value)),
      score = as.numeric(as.character(score)),
      n_methods = tidyr::replace_na(n_methods, 0)
    ) %>%
    group_by(TaxaID) %>%
    arrange(adjusted_p_value, desc(abs(score)), .by_group = TRUE) %>%
    filter(row_number() == 1) %>%
    ungroup()
}