# Ensure the pals package is loaded
# library(pals) # Or use pals:: explicitly

plot_combined_volcano <- function(
    volcano_df,
    title = "Combined Volcano Plot",
    out_path = NULL,
    label_top_n = 20
) {
  message("🔍 Starting plot_combined_volcano()")
  
  # ── Validate input ──
  required_cols <- c("score", "adjusted_p_value", "tax_group", "TaxaID", "is_significant", "n_methods", "avg_relab")
  missing_cols <- setdiff(required_cols, names(volcano_df))
  if (length(missing_cols) > 0) {
    message("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
    return(invisible(NULL))
  }
  
  # ── Diagnostics before filtering ──
  message("🔢 Rows before filtering: ", nrow(volcano_df))
  message("🧪 Non-NA score: ", sum(!is.na(volcano_df$score)))
  message("🧪 Non-NA adjusted_p_value: ", sum(!is.na(volcano_df$adjusted_p_value)))
  message("🧪 Range of score: ", paste(range(volcano_df$score, na.rm = TRUE), collapse = " to "))
  message("🧪 Range of adjusted_p_value: ", paste(range(volcano_df$adjusted_p_value, na.rm = TRUE), collapse = " to "))
  message("🧪 Non-NA avg_relab: ", sum(!is.na(volcano_df$avg_relab)))
  # Add diagnostic for tax_group before filtering
  message("🧪 Unique tax_group values before filtering: ", paste(na.omit(unique(volcano_df$tax_group)), collapse = ", "))
  
  
  # ── Filter valid rows ──
  volcano_df <- volcano_df %>%
    filter(!is.na(avg_relab), avg_relab > 0, is.finite(score), is.finite(adjusted_p_value))
  
  message("✅ Rows after filtering: ", nrow(volcano_df))
  if (nrow(volcano_df) == 0) {
    message("⚠️ No valid data after filtering. Skipping plot.")
    return(invisible(NULL))
  }
  
  # ── Calculate y-axis values ──
  volcano_df <- volcano_df %>%
    mutate(neg_log_p = -log10(adjusted_p_value))
  
  # ── Check ranges to avoid viewport errors ──
  if (diff(range(volcano_df$score, na.rm = TRUE)) == 0) {
    message("⚠️ Skipping plot: constant 'score' values.")
    return(invisible(NULL))
  }
  if (diff(range(volcano_df$neg_log_p, na.rm = TRUE)) == 0) {
    message("⚠️ Skipping plot: constant -log10(p) values.")
    return(invisible(NULL))
  }
  
  # ── Label top taxa ──
  label_df <- volcano_df %>%
    filter(is_significant) %>%
    arrange(desc(n_methods), desc(avg_relab), desc(abs(score))) %>%
    slice_head(n = label_top_n)
  
  message("🏷️ Labeling top ", nrow(label_df), " taxa")
  
  # ── Handle color palette ──
  unique_tax_groups <- na.omit(unique(volcano_df$tax_group))
  
  # --- Added Diagnostics for Color Palette ---
  message("🎨 Unique Tax Groups found *after filtering*: ", paste(unique_tax_groups, collapse = ", "))
  message("🎨 Number of unique Tax Groups *after filtering*: ", length(unique_tax_groups))
  # --- End Added Diagnostics ---
  
  custom_colors <- rep("grey80", length(unique_tax_groups))
  names(custom_colors) <- unique_tax_groups
  if (length(unique_tax_groups) > 0) {
    raw_colors <- pals::polychrome(n = max(30, length(unique_tax_groups)))
    if (length(raw_colors) >= length(unique_tax_groups)) {
      custom_colors <- setNames(raw_colors[1:length(unique_tax_groups)], unique_tax_groups)
    } else {
      # If somehow not enough colors are generated, warn and use what's available
      warning(sprintf("pals::polychrome generated %d colors, less than the %d needed unique groups.", length(raw_colors), length(unique_tax_groups)))
      # Use generated colors for as many groups as possible
      custom_colors <- setNames(raw_colors, unique_tax_groups[1:length(raw_colors)])
      # Remaining groups will get na.value
    }
  }
  
  # --- Added Diagnostics for Generated Colors ---
  message("🎨 Number of colors generated for palette: ", length(custom_colors))
  message("🎨 First 5 generated colors: ", paste(head(custom_colors, 5), collapse = ", "))
  message("🎨 Names of first 5 generated colors: ", paste(head(names(custom_colors), 5), collapse = ", "))
  # --- End Added Diagnostics ---
  
  
  # ── Plot ──
  p <- ggplot(volcano_df, aes(x = score, y = neg_log_p)) +
    geom_point(aes(color = tax_group, shape = factor(n_methods), size = avg_relab), alpha = 0.75) +
    geom_text_repel(data = label_df, aes(label = TaxaID), size = 2.5, max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 8), drop = FALSE) +
    scale_color_manual(values = custom_colors, na.value = "grey80", name = "Tax Group") +
    scale_size_continuous(
      range = c(1.5, 6),
      limits = if (length(label_df$avg_relab) > 1) range(label_df$avg_relab, na.rm = TRUE) else c(1.5, 6),
      trans = "sqrt"
    ) +
    labs(
      x = "Effect Size / Log2FC / LDA",
      y = "-log10 Adjusted P-value",
      shape = "# Methods",
      size = "Avg Relab.",
      title = title
    ) +
    theme_minimal() + theme(
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.position = "right"
  )
  
  # ── Print plot for interactive use ──
  if (interactive()) {
    print(p)
  }
  
  # ── Save plot if requested ──
  if (!is.null(out_path)) {
    tryCatch({
      ggsave(out_path, plot = p, width = 12, height = 6, dpi = 300)
      message("✅ Plot saved to: ", out_path)
    }, error = function(e) {
      warning("❌ Failed to save plot: ", e$message)
    })
  }
  
  message("✅ Finished plot_combined_volcano()")
  return(invisible(p))
}