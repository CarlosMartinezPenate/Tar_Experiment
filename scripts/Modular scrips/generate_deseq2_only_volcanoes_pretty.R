# ───────────────────────────────────────────────────────
# 🎨 DESeq2 Volcano + Rescued ASVs (Robust Version)
# File: generate_deseq2_only_volcanoes_pretty.R
# ───────────────────────────────────────────────────────
generate_deseq2_only_volcanoes_pretty <- function(
    results_dir,
    physeq,
    group_var,
    tax_level = "Genus",
    out_subdir = "deseq2_only_volcano_plots",
    label_top_n = 30
) {
  suppressPackageStartupMessages({
    library(phyloseq)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(ggrepel)
    library(stringr)
  })
  
  message("📂 Scanning results in: ", results_dir)
  files <- list.files(results_dir, pattern = "^results_asvs_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop("❌ No result files found.")
  
  out_dir <- file.path(results_dir, out_subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (file in files) {
    fname <- basename(file)
    parts <- str_match(fname, "^results_asvs_(.*?)_vs_(.*?)_(.*)\\.csv$")
    if (any(is.na(parts))) next
    
    contrast_group <- parts[2]
    control_group <- parts[3]
    filter_type <- parts[4]
    
    message("📊 Plotting DESeq2 Volcano (+rescued): ", contrast_group, " vs ", control_group, " (", filter_type, ")")
    
    df <- read.csv(file)
    
    # 🧠 Subset phyloseq object
    phy_sub <- prune_samples(
      sample_data(physeq)[[group_var]] %in% c(control_group, contrast_group) &
        sample_data(physeq)[["Filter_Type"]] == filter_type,
      physeq
    )
    
    # 🧮 Relative abundance
    relab <- as(otu_table(phy_sub), "matrix")
    if (!taxa_are_rows(phy_sub)) relab <- t(relab)
    relab <- sweep(relab, 2, colSums(relab), FUN = "/")
    relab <- as.data.frame(relab)
    
    meta <- data.frame(sample_data(phy_sub))
    control_samples <- rownames(meta)[meta[[group_var]] == control_group]
    contrast_samples <- rownames(meta)[meta[[group_var]] == contrast_group]
    
    relab_tbl <- tibble(
      TaxaID = rownames(relab),
      relab_control = rowMeans(relab[, control_samples, drop = FALSE]),
      relab_contrast = rowMeans(relab[, contrast_samples, drop = FALSE])
    ) %>%
      mutate(
        dominant_group = ifelse(relab_contrast > relab_control, "contrast", "control"),
        avg_relab = ifelse(dominant_group == "contrast", relab_contrast, relab_control)
      )
    
    # 🧬 Pull DESeq2, ALDEx2, LEfSe significant IDs
    aldex2_sig <- df %>% filter(source == "ALDEx2", is_significant) %>% pull(TaxaID)
    lefse_sig <- df %>% filter(source == "LEfSe", is_significant) %>% pull(TaxaID)
    
    # 🧬 Taxonomy table
    tax_df <- tax_table(physeq) %>%
      as.data.frame() %>%
      rownames_to_column("TaxaID")
    
    message("🧬 Tax Table Columns: ", paste(colnames(tax_df), collapse = ", "))
    
    # ⚠️ Remove existing tax columns to avoid duplicate renaming
    tax_cols <- colnames(tax_df)[-1]  # all except TaxaID
    df <- df %>% select(-any_of(tax_cols))
    
    # Merge DESeq2 + taxonomy + relab
    volcano_df <- df %>%
      filter(source == "DESeq2") %>%
      mutate(
        log2FoldChange = score,
        padj = adjusted_p_value,
        significant_in_aldex2 = TaxaID %in% aldex2_sig,
        significant_in_lefse = TaxaID %in% lefse_sig,
        method_overlap = case_when(
          significant_in_aldex2 & significant_in_lefse ~ "DESeq2 + ALDEx2 + LEfSe",
          significant_in_aldex2 ~ "DESeq2 + ALDEx2",
          significant_in_lefse ~ "DESeq2 + LEfSe",
          TRUE ~ "DESeq2 only"
        )
      ) %>%
      left_join(relab_tbl, by = "TaxaID") %>%
      left_join(tax_df, by = "TaxaID")
    
    # 🛑 Check for missing tax column
    if (!tax_level %in% colnames(volcano_df)) {
      stop(glue::glue("❌ Column for tax_level '{tax_level}' not found in final volcano_df.\nAvailable columns: {paste(colnames(volcano_df), collapse = ', ')}"))
    }
    
    volcano_df <- volcano_df %>%
      mutate(tax_group = .data[[tax_level]]) %>%
      filter(!is.na(padj), !is.na(log2FoldChange))
    
    volcano_df$method_overlap <- factor(
      volcano_df$method_overlap,
      levels = c("DESeq2 only", "DESeq2 + ALDEx2", "DESeq2 + LEfSe", "DESeq2 + ALDEx2 + LEfSe")
    )
    
    # 📏 Size scale training
    size_training <- volcano_df %>% filter(padj < 0.05) %>% pull(avg_relab)
    
    # 🎨 Plot
    p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = tax_group, shape = method_overlap, size = avg_relab), alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
      scale_color_brewer(palette = "Set3", na.value = "gray80") +
      scale_shape_manual(values = c(16, 17, 15, 8), drop = FALSE) +
      scale_size_continuous(range = c(1.5, 6), trans = "sqrt", limits = range(size_training, na.rm = TRUE)) +
      labs(
        title = sprintf("DESeq2 Volcano (+rescued): %s vs %s (%s)", contrast_group, control_group, filter_type),
        x = "Log₂ Fold Change",
        y = "-Log₁₀ Adjusted P",
        color = paste(tax_level, "(Tax Group)"),
        shape = "Method Overlap",
        size = "Avg Relab"
      ) +
      ggrepel::geom_text_repel(
        data = volcano_df %>%
          filter(padj < 0.05) %>%
          arrange(desc(abs(log2FoldChange))) %>%
          slice_head(n = label_top_n),
        aes(label = TaxaID),
        size = 2.5,
        max.overlaps = 20
      ) +
      theme_minimal()
    
    # 💾 Save
    out_png <- file.path(out_dir, sprintf("deseq2_volcano_%s_vs_%s_%s.png", contrast_group, control_group, filter_type))
    out_csv <- file.path(out_dir, sprintf("deseq2_data_%s_vs_%s_%s.csv", contrast_group, control_group, filter_type))
    
    ggsave(out_png, plot = p, width = 8, height = 6, dpi = 300)
    write.csv(volcano_df, out_csv, row.names = FALSE)
    message("✅ Saved: ", out_png)
  }
}