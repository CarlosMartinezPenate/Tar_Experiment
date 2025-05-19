generate_complex_heatmap_from_flagged <- function(flagged_df,
                                                  physeq,
                                                  control_group,
                                                  contrast_group,
                                                  group_var,
                                                  filter_type,
                                                  transform_type = "zscore",
                                                  tax_level = "Genus",
                                                  robust_only = TRUE,
                                                  output_dir = "./complex_heatmaps",
                                                  min_methods = 2,
                                                  top_n_taxa = 50) {
  # Load dependencies
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(phyloseq)
  })
  
  message("🧼 Starting ComplexHeatmap heatmap generation...")
  message("🧾 Saving to directory: ", output_dir)
  cat("📊 Flagged DF dimensions: ", dim(flagged_df), "\n")
  
  # Filter for robust taxa
  if (robust_only) {
    if (!("n_methods" %in% colnames(flagged_df))) stop("Missing 'n_methods' in flagged_df")
    flagged_df <- flagged_df[flagged_df$is_significant & flagged_df$n_methods >= min_methods, , drop = FALSE]
    cat("🔎 Filtered robust taxa (n_methods >= ", min_methods, "): ", nrow(flagged_df), " remaining\n")
    if (nrow(flagged_df) == 0) {
      message("⚠️ No robust significant taxa. Skipping.")
      return(NULL)
    }
  }
  
  sig_taxa <- unique(flagged_df$TaxaID)
  physeq_sig <- prune_taxa(sig_taxa, physeq)
  cat("📦 Filtered phyloseq → taxa: ", ntaxa(physeq_sig), ", samples: ", nsamples(physeq_sig), "\n")
  
  if (ntaxa(physeq_sig) < 2 || nsamples(physeq_sig) < 2) {
    message("⚠️ Not enough data. Skipping heatmap.")
    return(NULL)
  }
  
  # Extract OTU matrix
  mat <- as(otu_table(physeq_sig), "matrix")
  if (!taxa_are_rows(physeq_sig)) mat <- t(mat)
  cat("✅ Extracted OTU matrix: ", dim(mat)[1], " taxa × ", dim(mat)[2], " samples\n")
  
  # Transform to relative abundance
  mat <- sweep(mat, 2, colSums(mat), "/")
  mat[is.na(mat)] <- 0
  cat("🔄 Transformed to relative abundance\n")
  
  # Apply transformation
  if (transform_type == "zscore") {
    mat <- t(scale(t(mat)))
    cat("🔁 Applied z-score transformation\n")
  } else if (transform_type == "clr") {
    mat <- log(mat + 1e-6)
    mat <- sweep(mat, 2, colMeans(mat), "-")
    cat("🔁 Applied CLR transformation\n")
  } else if (transform_type == "log2") {
    mat <- log2(mat + 1e-6)
    cat("🔁 Applied log2 transformation\n")
  } else if (transform_type == "arcsine") {
    mat <- asin(sqrt(mat))
    cat("🔁 Applied arcsine transformation\n")
  } else if (transform_type != "ra") {
    warning("❌ Unknown transform_type: ", transform_type)
  }
  
  # Top N variable taxa
  if (!is.null(top_n_taxa) && nrow(mat) > top_n_taxa) {
    variances <- apply(mat, 1, var)
    top_taxa <- names(sort(variances, decreasing = TRUE))[1:top_n_taxa]
    mat <- mat[top_taxa, , drop = FALSE]
    cat("🔪 Selected top ", top_n_taxa, " variable taxa\n")
  }
  
  # Taxa labels
  tax_df <- as.data.frame(tax_table(physeq_sig))
  taxa_labels <- paste0(rownames(mat), " | ", tax_df[rownames(mat), tax_level])
  rownames(mat) <- taxa_labels
  cat("🧬 Constructed taxa labels with level: ", tax_level, "\n")
  
  # Metadata
  meta <- as.data.frame(sample_data(physeq_sig))
  meta$SampleID <- rownames(meta)
  cat("📋 Sample metadata dimensions: ", dim(meta), "\n")
  
  # Group annotations
  group_annotation <- meta[, group_var]
  names(group_annotation) <- meta$SampleID
  group_annotation <- group_annotation[colnames(mat)]
  ha_col <- ComplexHeatmap::HeatmapAnnotation(
    Group = group_annotation,
    col = list(Group = structure(rainbow(length(unique(group_annotation)), end = 0.8),
                                 names = unique(group_annotation)))
  )
  
  # Asterisks for significant taxa
  sig_taxa_names <- flagged_df$TaxaID
  row_labels <- rownames(mat)
  row_labels_with_asterisk <- ifelse(
    gsub(" \\|.*", "", row_labels) %in% sig_taxa_names,
    paste0(row_labels, "*"),
    row_labels
  )
  
  # Heatmap
  cat("🎨 Drawing heatmap...\n")
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = transform_type,
    top_annotation = ha_col,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 8, rot = 45),
    row_labels = row_labels_with_asterisk,
    col = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("navy", "white", "firebrick")),
    heatmap_legend_param = list(title = transform_type)
  )
  
  # Save
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(output_dir, sprintf("heatmap_%s_vs_%s_%s_%s.png",
                                            contrast_group, control_group,
                                            filter_type, transform_type))
  png(out_path, width = 1200, height = 900, res = 150)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  message("✅ Saved ComplexHeatmap to: ", out_path)
}