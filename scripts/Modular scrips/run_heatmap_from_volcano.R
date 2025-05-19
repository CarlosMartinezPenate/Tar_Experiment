run_heatmap_from_volcano <- function(
    physeq_obj,
    volcano_csv_path,
    transformation = "zscore",
    output_dir,
    html_dir,
    tax_level = "Genus",
    filter_type = NULL  # ✅ Added argument
) {
  # Load necessary libraries
  library(phyloseq)
  library(ggplot2)
  library(heatmaply)
  library(tibble)
  library(dplyr)
  
  # Read volcano data
  volcano_df <- read.csv(volcano_csv_path)
  
  # Extract ASVs to highlight
  sig_asvs <- volcano_df %>%
    filter(significant == TRUE) %>%
    pull(TaxaID) %>%
    unique()
  
  if (length(sig_asvs) == 0) {
    message("⚠️ No significant ASVs found in: ", basename(volcano_csv_path))
    return(NULL)
  }
  
  # Subset phyloseq to significant ASVs
  physeq_sub <- prune_taxa(sig_asvs, physeq_obj)
  
  # Get transformed abundance matrix
  otu_mat <- switch(
    transformation,
    zscore = scale(otu_table(physeq_sub)),
    clr = {
      library(compositions)
      clr(t(otu_table(physeq_sub)) + 1e-6)  # add pseudo count
    },
    ra = prop.table(otu_table(physeq_sub), 2),
    log2 = log2(otu_table(physeq_sub) + 1),
    arcsine = asin(sqrt(prop.table(otu_table(physeq_sub), 2))),
    vst = {
      library(DESeq2)
      dds <- phyloseq_to_deseq2(physeq_sub, ~ 1)
      vst_counts <- varianceStabilizingTransformation(dds)
      assay(vst_counts)
    },
    stop("❌ Unknown transformation method: ", transformation)
  )
  
  # Convert to matrix
  heatmap_mat <- as.matrix(otu_mat)
  
  # Prepare plot title and filenames
  suffix <- if (!is.null(filter_type)) filter_type else "all"
  plot_title <- paste("Heatmap:", transformation, "—", suffix)
  fname_base <- paste0("heatmap_", transformation, "_", suffix)
  
  # Define save paths
  png_path <- file.path(output_dir, paste0(fname_base, ".png"))
  html_path <- file.path(html_dir, paste0(fname_base, ".html"))
  
  # Generate heatmap (static)
  heatmap_plot <- pheatmap::pheatmap(
    heatmap_mat,
    main = plot_title,
    fontsize = 10,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean"
  )
  
  # Save static plot
  ggsave(filename = png_path, plot = heatmap_plot$gtable, width = 10, height = 8)
  
  # Generate interactive plot
  heatmaply_plot <- heatmaply::heatmaply(
    heatmap_mat,
    main = plot_title,
    file = html_path
  )
  
  message("✅ Heatmap saved: ", png_path)
}