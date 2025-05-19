# ───────────────────────────────────────────────────────
# 📦 Load Required Libraries
# ───────────────────────────────────────────────────────
library(phyloseq)
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(heatmaply)
library(stringr)
library(RColorBrewer)
library(htmlwidgets)

# ───────────────────────────────────────────────────────
# 🔧 Transform Abundance Matrix Function
# ───────────────────────────────────────────────────────
transform_abundance <- function(physeq_obj, method = "zscore") {
  relab <- transform_sample_counts(physeq_obj, function(x) x / sum(x))
  mat <- as(otu_table(relab), "matrix")
  if (!taxa_are_rows(relab)) mat <- t(mat)
  
  if (method == "zscore") {
    mat <- t(scale(t(mat)))
  } else if (method == "clr") {
    pseudocount <- 1e-6
    mat <- log(mat + pseudocount)
    mat <- sweep(mat, 2, colMeans(mat, na.rm = TRUE))
  } else if (method == "ra") {
    # Already relative abundance
  } else if (method == "log2") {
    pseudocount <- 1e-6
    mat <- log2(mat + pseudocount)
  } else if (method == "arcsine") {
    mat <- asin(sqrt(mat))
  } else if (method == "vst") {
    dds <- phyloseq_to_deseq2(physeq_obj, ~1)
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- estimateDispersions(dds, fitType = "local")
    vsd <- varianceStabilizingTransformation(dds)
    mat <- assay(vsd)
    if (!taxa_are_rows(physeq_obj)) mat <- t(mat)
  } else {
    stop("❌ Invalid transformation method.")
  }
  
  mat[is.na(mat)] <- 0
  return(mat)
}

# ───────────────────────────────────────────────────────
# 🧬 Extract Robust Taxa From Volcano File
# ───────────────────────────────────────────────────────
get_robust_taxa <- function(volcano_csv_path, min_methods = 2) {
  df <- read_csv(volcano_csv_path, show_col_types = FALSE)
  cat("\n📂 DEBUG: Volcano file read (", basename(volcano_csv_path), ") - rows: ", nrow(df), "\n")
  
  if ("n_methods" %in% colnames(df)) {
    robust <- df %>% filter(!is.na(n_methods), n_methods >= min_methods, is_significant)
  } else if ("method_overlap" %in% colnames(df)) {
    robust <- df %>% filter(!is.na(method_overlap), method_overlap != "DESeq2 only", is_significant)
  } else {
    stop("❌ Required column missing in: ", volcano_csv_path)
  }
  
  return(unique(robust$TaxaID))
}

# ───────────────────────────────────────────────────────
# 🔥 Main Function: Robust Taxa Heatmap Generator
# ───────────────────────────────────────────────────────
run_heatmap_from_volcano <- function(
    physeq_obj,
    volcano_csv_path,
    transformation = "zscore",
    output_dir = "./heatmap_outputs_robust_debug",
    html_dir = "./interactive_heatmaps_robust_debug",
    tax_level = "Genus"
) {
  fname <- basename(volcano_csv_path)
  parts <- stringr::str_match(fname, "deseq2_data_(.*?)_vs_(.*?)_(.*?)\\.csv")
  if (any(is.na(parts))) {
    message(sprintf("⚠️ Skipping malformed filename: %s", fname))
    return(NULL)
  }
  
  contrast_group <- parts[2]
  control_group  <- parts[3]
  filter_type    <- parts[4]
  
  cat("\n📊 RUNNING: ", contrast_group, " vs ", control_group, " (", filter_type, ")\n")
  
  physeq_sub <- subset_samples(physeq_obj, Filter_Type == filter_type)
  meta <- data.frame(sample_data(physeq_sub))
  keep_samples <- meta$Station_treatment %in% c(control_group, contrast_group)
  physeq_sub <- prune_samples(keep_samples, physeq_sub)
  physeq_sub <- prune_taxa(taxa_sums(physeq_sub) > 120, physeq_sub)
  
  if (nsamples(physeq_sub) < 4 || ntaxa(physeq_sub) < 2) {
    message(sprintf("⚠️ Not enough data for: %s vs %s (%s)", contrast_group, control_group, filter_type))
    return(NULL)
  }
  
  dds <- phyloseq_to_deseq2(physeq_sub, ~Station_treatment)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, fitType = "local")
  
  coef_name <- paste0("Station_treatment_", contrast_group, "_vs_", control_group)
  res <- if (coef_name %in% resultsNames(dds)) {
    lfcShrink(dds, coef = coef_name, type = "normal")
  } else {
    results(dds, contrast = c("Station_treatment", contrast_group, control_group))
  }
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("TaxaID") %>%
    mutate(
      adjusted_p_value = padj,
      score = log2FoldChange,
      is_significant = !is.na(padj) & padj < 0.05
    )
  
  tax_tbl <- as.data.frame(tax_table(physeq_sub))
  res_df <- left_join(res_df, rownames_to_column(tax_tbl, var = "TaxaID"), by = "TaxaID")
  
  robust_taxa <- get_robust_taxa(volcano_csv_path)
  sig_asvs <- intersect(robust_taxa, res_df$TaxaID)
  
  if (length(sig_asvs) == 0) {
    message(sprintf("❌ No robust taxa in: %s vs %s (%s)", contrast_group, control_group, filter_type))
    return(NULL)
  }
  
  physeq_sig <- prune_taxa(sig_asvs, physeq_sub)
  z_matrix <- transform_abundance(physeq_sig, method = transformation)
  z_matrix <- z_matrix[rowSums(z_matrix != 0) > 0, , drop = FALSE]
  
  if (nrow(z_matrix) == 0) {
    message(sprintf("⚠️ No nonzero robust taxa in: %s vs %s (%s)", contrast_group, control_group, filter_type))
    return(NULL)
  }
  
  tax_df <- as.data.frame(tax_table(physeq_sig))
  row_labels <- paste0(rownames(z_matrix), " | ", tax_df[rownames(z_matrix), tax_level])
  rownames(z_matrix) <- row_labels
  
  meta_df <- data.frame(sample_data(physeq_sig))
  meta_df$ID <- as.character(meta_df$ID)
  meta_df$Station_treatment <- as.character(meta_df$Station_treatment)
  meta_df$sample_id <- rownames(meta_df)
  meta_df$custom_label <- paste0(meta_df$ID, "_", meta_df$Station_treatment)
  colnames(z_matrix) <- meta_df$custom_label[match(colnames(z_matrix), meta_df$sample_id)]
  
  col_ann <- meta_df %>%
    filter(custom_label %in% colnames(z_matrix)) %>%
    select(sample_id, Station_treatment, custom_label) %>%
    rename(Group = Station_treatment)
  rownames(col_ann) <- col_ann$custom_label
  col_ann$sample_id <- col_ann$custom_label <- NULL
  
  sample_order <- rownames(col_ann)
  z_matrix <- z_matrix[, sample_order, drop = FALSE]
  
  group_colors <- setNames(brewer.pal(3, "Dark2")[1:2], unique(meta_df$Station_treatment))
  ann_colors <- list(Group = group_colors)
  
  file_prefix <- sprintf("heatmap_robust_%s_vs_%s_%s_by_%s", contrast_group, control_group, filter_type, transformation)
  png_out <- file.path(output_dir, paste0(file_prefix, ".png"))
  csv_out <- file.path(output_dir, paste0(file_prefix, "_", transformation, ".csv"))
  html_out <- file.path(html_dir, paste0(file_prefix, ".html"))
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(html_dir, showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    z_matrix_clean <- z_matrix[rowSums(is.finite(z_matrix)) == ncol(z_matrix), , drop = FALSE]
    
    if (nrow(z_matrix_clean) < 2 || ncol(z_matrix_clean) < 2) {
      message(sprintf("⚠️ Not enough valid rows/columns for: %s vs %s (%s)", contrast_group, control_group, filter_type))
      return(NULL)
    }
    
    png(png_out, width = max(1000, ncol(z_matrix_clean) * 30), height = 800, res = 150)
    pheat_out <- pheatmap(
      z_matrix_clean,
      annotation_col = col_ann,
      annotation_colors = ann_colors,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_colnames = TRUE,
      fontsize_row = 6,
      fontsize_col = 8,
      color = colorRampPalette(brewer.pal(11, "BrBG"))(100),
      main = sprintf(
        "Robust Taxa Heatmap: %s vs %s (%s)\nTransformed by: %s",
        contrast_group, control_group, filter_type, transformation
      )
    )
    dev.off()
    write.csv(z_matrix_clean, csv_out)
    
    z_html <- z_matrix_clean[pheat_out$tree_row$order, , drop = FALSE]
    html_heatmap <- heatmaply::heatmaply(
      z_html,
      col = colorRampPalette(brewer.pal(11, "BrBG"))(100),
      Rowv = TRUE,
      Colv = FALSE,
      labRow = rownames(z_html),
      labCol = colnames(z_html),
      xlab = "Sample",
      ylab = "ASV | Genus",
      fontsize_row = 8,
      fontsize_col = 10,
      main = sprintf("Interactive Robust Heatmap: %s vs %s (%s)", contrast_group, control_group, filter_type),
      height = 400 + 18 * nrow(z_html),
      width = 1000
    )
    htmlwidgets::saveWidget(html_heatmap, file = html_out, selfcontained = TRUE)
    message(sprintf("✅ Heatmaps done: %s vs %s (%s)", contrast_group, control_group, filter_type))
  }, error = function(e) {
    message(sprintf("❌ Error generating heatmap for %s vs %s (%s): %s", contrast_group, control_group, filter_type, e$message))
  })
}