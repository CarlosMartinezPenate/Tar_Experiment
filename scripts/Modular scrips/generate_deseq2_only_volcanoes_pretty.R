generate_deseq2_only_volcanoes_pretty <- function(
    results_dir,
    physeq,
    group_var,
    tax_level = "Genus",
    out_subdir = "deseq2_only_volcano_plots",
    label_top_n = 30
) {
  message("Starting Rescued Volcano")
  suppressPackageStartupMessages({
    library(phyloseq)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(ggrepel)
    library(stringr)
    library(pals)
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
    
    # Subset phyloseq
    phy_sub <- prune_samples(
      sample_data(physeq)[[group_var]] %in% c(control_group, contrast_group) &
        sample_data(physeq)[["Filter_Type"]] == filter_type,
      physeq
    )
    
    if (nsamples(phy_sub) == 0) {
      message("⚠️ Skipping: No samples after subsetting.")
      next
    }
    
    relab <- as(otu_table(phy_sub), "matrix")
    if (!taxa_are_rows(phy_sub)) relab <- t(relab)
    relab <- sweep(relab, 2, colSums(relab), FUN = "/")
    relab <- as.data.frame(relab)
    
    meta <- data.frame(sample_data(phy_sub))
    control_samples <- rownames(meta)[meta[[group_var]] == control_group]
    contrast_samples <- rownames(meta)[meta[[group_var]] == contrast_group]
    
    avg_relab_control <- if (length(control_samples) > 0) rowMeans(relab[, control_samples, drop = FALSE]) else rep(0, nrow(relab))
    avg_relab_contrast <- if (length(contrast_samples) > 0) rowMeans(relab[, contrast_samples, drop = FALSE]) else rep(0, nrow(relab))
    
    relab_tbl <- tibble(
      TaxaID = rownames(relab),
      relab_control = avg_relab_control,
      relab_contrast = avg_relab_contrast
    ) %>%
      mutate(
        dominant_group = ifelse(relab_contrast > relab_control, "contrast", "control"),
        avg_relab = ifelse(dominant_group == "contrast", relab_contrast, relab_control)
      )
    
    aldex2_sig <- df %>% filter(source == "ALDEx2", is_significant) %>% pull(TaxaID)
    lefse_sig <- df %>% filter(source == "LEfSe", is_significant) %>% pull(TaxaID)
    
    tax_df <- tax_table(physeq) %>%
      as.data.frame() %>%
      rownames_to_column("TaxaID")
    
    tax_cols <- colnames(tax_df)[-1]
    df <- df %>% dplyr::select(-dplyr::any_of(tax_cols))
    
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
    
    if (!tax_level %in% colnames(volcano_df)) {
      stop(glue::glue("❌ Column for tax_level '{tax_level}' not found in final volcano_df."))
    }
    
    volcano_df <- volcano_df %>%
      mutate(tax_group = .data[[tax_level]]) %>%
      filter(!is.na(padj), !is.na(log2FoldChange))
    
    if (nrow(volcano_df) == 0 ||
        diff(range(volcano_df$log2FoldChange, na.rm = TRUE)) == 0 ||
        diff(range(-log10(volcano_df$padj), na.rm = TRUE)) == 0) {
      message("⚠️ Skipping plot: Insufficient variance in data.")
      next
    }
    
    volcano_df$method_overlap <- factor(
      volcano_df$method_overlap,
      levels = c("DESeq2 only", "DESeq2 + ALDEx2", "DESeq2 + LEfSe", "DESeq2 + ALDEx2 + LEfSe")
    )
    
    unique_tax_groups <- na.omit(unique(volcano_df$tax_group))
    raw_colors <- pals::polychrome(n = max(30, length(unique_tax_groups)))
    custom_colors <- setNames(raw_colors[1:length(unique_tax_groups)], unique_tax_groups)
    
    volcano_df <- volcano_df %>%
      mutate(neg_log_p = -log10(padj))
    
    size_training <- volcano_df %>% filter(padj < 0.05) %>% pull(avg_relab)
    
    p <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log_p)) +
      geom_point(aes(color = tax_group, shape = method_overlap, size = avg_relab), alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
      scale_color_manual(values = custom_colors, na.value = "gray80", name = paste(tax_level, "(Tax Group)")) +
      scale_shape_manual(values = c(16, 17, 15, 8), drop = FALSE) +
      scale_size_continuous(
        range = c(1.5, 6),
        limits = if (length(size_training) > 0 && diff(range(size_training, na.rm = TRUE)) > 0) {
          range(size_training, na.rm = TRUE)
        } else {
          range(volcano_df$avg_relab, na.rm = TRUE)
        },
        trans = "sqrt"
      ) +
      labs(
        title = sprintf("DESeq2 Volcano (+rescued): %s vs %s (%s)", contrast_group, control_group, filter_type),
        x = "Log₂ Fold Change",
        y = "-Log₁₀ Adjusted P",
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
      theme_minimal() +
      theme(
        legend.key.size = unit(0.4, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)
      )
    
    print(p)
    
    out_png <- file.path(out_dir, sprintf("deseq2_volcano_%s_vs_%s_%s.png", contrast_group, control_group, filter_type))
    out_csv <- file.path(out_dir, sprintf("deseq2_data_%s_vs_%s_%s.csv", contrast_group, control_group, filter_type))
    
    tryCatch({
      ggsave(out_png, plot = p, width = 12, height = 6, dpi = 300)
      message("✅ Saved plot: ", out_png)
    }, error = function(e) {
      warning(paste("❌ Error saving plot to", out_png, ":", e$message))
    })
    
    write.csv(volcano_df, out_csv, row.names = FALSE)
    message("✅ Saved data: ", out_csv)
  }
  
  message("Finished Rescued Volcano")
}