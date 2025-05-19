# ────────────────────────────────────────────────────────────────
# 🌡️ Modular Heatmap Generator — ggplot-based
# File: generate_heatmap_from_flagged_ggplot.R
# ────────────────────────────────────────────────────────────────
generate_heatmap_from_flagged <- function(flagged_df,
                                          physeq,
                                          control_group,
                                          contrast_group,
                                          group_var,
                                          filter_type,
                                          transform_type = "log",
                                          robust_only = TRUE,
                                          top_n_taxa = 50,
                                          save_path = "heatmap.png") {
  message("🧼 Starting generate_heatmap_from_flagged()")
  message("🧾 Intended save_path: ", save_path)
  message("📋 Input flagged_df dimensions: ", paste(dim(flagged_df), collapse = " x "))
  
  # 🔍 Filter robust taxa
  if (robust_only && "is_significant" %in% names(flagged_df)) {
    flagged_df <- dplyr::filter(flagged_df, is_significant)
    if (nrow(flagged_df) == 0) {
      warning("⚠️ No significant taxa remain after filtering. Skipping heatmap.")
      return(NULL)
    }
  }
  
  sig_taxa <- unique(flagged_df$TaxaID)
  physeq_sig <- prune_taxa(sig_taxa, physeq)
  
  if (ntaxa(physeq_sig) == 0 || nsamples(physeq_sig) == 0) {
    warning("⚠️ No taxa or samples in subset phyloseq object. Skipping heatmap.")
    return(NULL)
  }
  
  # 🧮 Extract and orient matrix
  mat <- as(otu_table(physeq_sig), "matrix")
  if (!taxa_are_rows(physeq_sig)) {
    mat <- t(mat)
  }
  
  # 🚿 Relative abundance
  mat <- sweep(mat, 2, colSums(mat), FUN = "/")
  mat[is.na(mat)] <- 0
  
  # 🔄 Transform
  if (transform_type == "log") {
    mat <- log1p(mat)
  } else if (transform_type == "zscore") {
    mat <- t(scale(t(mat)))
  } else if (transform_type != "ra") {
    warning("⚠️ Unknown transform_type. Defaulting to relative abundance.")
  }
  
  # 🔪 Limit to top N variable taxa
  if (!is.null(top_n_taxa) && top_n_taxa > 0 && nrow(mat) > top_n_taxa) {
    vars <- apply(mat, 1, var, na.rm = TRUE)
    top_taxa <- names(sort(vars, decreasing = TRUE))[1:top_n_taxa]
    mat <- mat[top_taxa, , drop = FALSE]
  }
  
  # 📦 Convert to tidy format
  df_long <- as.data.frame(mat) %>%
    tibble::rownames_to_column("TaxaID") %>%
    tidyr::pivot_longer(-TaxaID, names_to = "SampleID", values_to = "Abundance")
  
  # 📋 Sample metadata
  meta <- data.frame(sample_data(physeq_sig)) %>%
    tibble::rownames_to_column("SampleID")
  
  if (!group_var %in% names(meta)) {
    stop("❌ group_var not found in sample_data()")
  }
  
  df_long <- dplyr::left_join(df_long, meta[, c("SampleID", group_var)], by = "SampleID")
  colnames(df_long)[colnames(df_long) == group_var] <- "Group"
  
  # 🧭 Ordering
  taxa_order <- df_long %>%
    dplyr::group_by(TaxaID) %>%
    dplyr::summarise(var = var(Abundance, na.rm = TRUE)) %>%
    dplyr::arrange(desc(var)) %>%
    dplyr::pull(TaxaID)
  
  df_long$TaxaID <- factor(df_long$TaxaID, levels = taxa_order)
  df_long$SampleID <- factor(df_long$SampleID)
  
  # 🎨 Plot
  p <- ggplot(df_long, aes(x = SampleID, y = TaxaID, fill = Abundance)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient2(
      low = "navy", mid = "white", high = "firebrick",
      midpoint = 0, name = transform_type
    ) +
    theme_minimal(base_size = 10) +
    labs(
      title = sprintf("Heatmap: %s vs %s (%s)", contrast_group, control_group, filter_type),
      x = "Sample",
      y = "Taxa"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank()
    )
  
  # 💾 Save
  tryCatch({
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300)
    message("✅ Heatmap saved to: ", save_path)
  }, error = function(e) {
    warning("❌ Failed to save heatmap: ", e$message)
  })
  
  invisible(p)
}