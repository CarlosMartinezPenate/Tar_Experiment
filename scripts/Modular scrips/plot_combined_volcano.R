plot_combined_volcano <- function(volcano_df,
                                  output_path,
                                  title = "Volcano Plot",
                                  label_top_n = 10,
                                  point_size = 1.5,
                                  text_size = 3) {
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # Ensure required columns are present
  required_cols <- c("TaxaID", "score", "adjusted_p_value", "is_significant")
  missing_cols <- setdiff(required_cols, colnames(volcano_df))
  if (length(missing_cols) > 0) {
    stop("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Subset top N labels by absolute score (e.g., |log2FC|)
  top_labels <- volcano_df %>%
    filter(is_significant) %>%
    arrange(desc(abs(score))) %>%
    slice_head(n = label_top_n)
  
  # Create base plot
  p <- ggplot(volcano_df, aes(x = score, y = -log10(adjusted_p_value))) +
    geom_point(aes(color = is_significant), size = point_size, alpha = 0.7) +
    scale_color_manual(values = c("gray60", "#D55E00")) +
    theme_minimal(base_size = 12) +
    labs(
      x = "Log2 Fold Change (score)",
      y = expression(-log[10](adj.~p~value)),
      title = title,
      color = "Significant"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
  
  # Add top taxa labels
  if (nrow(top_labels) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = TaxaID),
      size = text_size,
      max.overlaps = Inf,
      box.padding = 0.4,
      segment.color = "grey50"
    )
  }
  
  # Save plot
  ggsave(filename = output_path, plot = p, width = 8, height = 6, dpi = 300)
  message("✅ Saved combined volcano to: ", output_path)
}