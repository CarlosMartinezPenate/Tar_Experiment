prepare_volcano_data <- function(flagged_df, physeq, group_var, control_group, contrast_group, tax_level = "Genus") {
  message("🔍 Starting prepare_volcano_data()")
  
  # ── Relative Abundance Calculation ──
  relab_df <- otu_table(physeq)
  if (!taxa_are_rows(physeq)) relab_df <- t(relab_df)
  relab_df <- sweep(relab_df, 2, colSums(relab_df), FUN = "/")
  relab_df <- as.data.frame(relab_df)
  
  meta <- data.frame(sample_data(physeq))
  control_samples <- rownames(meta)[meta[[group_var]] == control_group]
  contrast_samples <- rownames(meta)[meta[[group_var]] == contrast_group]
  
  avg_relab_control <- rowMeans(relab_df[, control_samples, drop = FALSE])
  avg_relab_contrast <- rowMeans(relab_df[, contrast_samples, drop = FALSE])
  
  relab_df_combined <- tibble(
    TaxaID = rownames(relab_df),
    relab_control = avg_relab_control,
    relab_contrast = avg_relab_contrast
  ) %>%
    mutate(
      dominant_group = ifelse(relab_contrast > relab_control, "contrast", "control"),
      avg_relab = ifelse(dominant_group == "contrast", relab_contrast, relab_control)
    )
  
  # ── Method Count ──
  n_methods_df <- flagged_df %>%
    filter(is_significant) %>%
    group_by(TaxaID) %>%
    summarise(n_methods = n_distinct(method), .groups = "drop")
  
  # ── Tax Table Mapping ──
  message("📦 Converting tax_table to matrix and then tibble...")
  tax_df <- as(tax_table(physeq), "matrix") %>%
    as.data.frame() %>%
    rownames_to_column("TaxaID") %>%
    tibble()
  
  message("📏 Tax table dimensions: ", paste(dim(tax_df), collapse = " × "))
  cat("✅ Tax data preview:\n")
  print(head(tax_df[, c("TaxaID", tax_level)], 3))
  
  tax_df <- tax_df %>%
    mutate(tax_group = .data[[tax_level]]) %>%
    dplyr::select(TaxaID, tax_group)
  
  # ── Join All ──
  joined_df <- flagged_df %>%
    left_join(relab_df_combined, by = "TaxaID") %>%
    left_join(n_methods_df, by = "TaxaID") %>%
    left_join(tax_df, by = "TaxaID")
  
  # ── Final Clean and Rank ──
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