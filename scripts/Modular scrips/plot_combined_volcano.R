
# Ensure the pals package is loaded
library(pals) 

# Or use pals:: explicitly

plot_combined_volcano <- function(
    volcano_df,
    title = "Combined Volcano Plot",
    out_path = NULL,
    label_top_n = 20
) {
  # Filter NA/zero relab
  volcano_df <- volcano_df %>%
    filter(!is.na(avg_relab), avg_relab > 0)
  
  # Check essential columns exist and are finite
  if (!all(c("score", "adjusted_p_value") %in% names(volcano_df))) {
    message("⚠️ Skipping plot: required columns missing.")
    return(invisible(NULL))
  }
  
  volcano_df <- volcano_df %>%
    filter(is.finite(score), is.finite(adjusted_p_value))
  
  if (nrow(volcano_df) == 0) {
    message("⚠️ Skipping plot: no valid rows to plot.")
    return(invisible(NULL))
  }
  
  label_df <- volcano_df %>%
    filter(is_significant) %>%
    arrange(desc(n_methods), desc(avg_relab), desc(abs(score))) %>%
    slice_head(n = label_top_n)
  
  size_training_data <- volcano_df %>%
    filter(is_significant) %>%
    pull(avg_relab)
  
  p <- ggplot(volcano_df, aes(x = score, y = -log10(adjusted_p_value))) +
    geom_point(aes(color = tax_group, shape = factor(n_methods), size = avg_relab), alpha = 0.75) +
    geom_text_repel(data = label_df, aes(label = TaxaID), size = 2.5, max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 8), drop = FALSE) +
    scale_color_brewer(palette = "Set3", na.value = "grey80") +
    scale_size_continuous(
      range = c(1.5, 6),
      limits = range(size_training_data, na.rm = TRUE),
      trans = "sqrt"
    ) +
    labs(
      x = "Effect Size / Log2FC / LDA",
      y = "-log10 Adjusted P-value",
      shape = "# Methods",
      size = "Avg Relab.",
      color = "Tax Group",
      title = title
    ) +

    theme_minimal() + theme(
      legend.key.size = unit(0.4, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.position = "right"
    )

    theme_minimal()

    theme_minimal()

    theme_minimal()

  
  print(p)
  
  if (!is.null(out_path)) {
    ggsave(out_path, p, width = 8, height = 6, dpi = 300)
  }
}