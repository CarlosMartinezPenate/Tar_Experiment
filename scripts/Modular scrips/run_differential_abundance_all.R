run_differential_abundance_all <- function(
    physeq_obj,
    group_var,
    contrast_group,
    control_group,
    filter_type,
    method_list,
    thresholds,
    output_dir,
    plot = TRUE,
    verbose = TRUE
) {
  message("starting run_differential_abundance_all()")
  
  is_plot_valid <- function(df, x_col = "score", y_col = "adjusted_p_value") {
    if (!all(c(x_col, y_col) %in% names(df))) return(FALSE)
    df <- df %>%
      filter(
        !is.na(.data[[x_col]]),
        !is.na(.data[[y_col]]),
        is.finite(.data[[x_col]]),
        is.finite(.data[[y_col]])
      )
    nrow(df) > 0
  }
  
  has_constant_taxa <- function(mat, group_labels) {
    apply(mat, 1, function(taxon_vals) {
      any(sapply(unique(group_labels), function(g) {
        vals <- taxon_vals[group_labels == g]
        length(unique(vals)) == 1
      }))
    })
  }
  
  # ── Prune Data ──
  meta <- data.frame(sample_data(physeq_obj))
  meta <- meta[meta[[group_var]] %in% c(control_group, contrast_group), ]
  meta[[group_var]] <- factor(meta[[group_var]], levels = c(control_group, contrast_group))
  sample_data(physeq_obj) <- sample_data(meta)
  physeq_obj <- prune_samples(rownames(meta), physeq_obj)
  physeq_obj <- prune_taxa(taxa_sums(physeq_obj) > 120, physeq_obj)
  
  if (nsamples(physeq_obj) == 0 || ntaxa(physeq_obj) == 0) {
    warning(sprintf("❌ Skipping %s vs %s (%s): no data after filtering", contrast_group, control_group, filter_type))
    return(tibble())
  }
  
  result_list <- list()
  
  # ── DESeq2 ──
  if ("DESeq2" %in% method_list) {
    dds <- phyloseq::phyloseq_to_deseq2(physeq_obj, as.formula(paste("~", group_var)))
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
    dds <- suppressMessages(DESeq2::DESeq(dds))
    
    coef_name <- paste0(group_var, "_", contrast_group, "_vs_", control_group)
    res <- if (coef_name %in% resultsNames(dds)) {
      DESeq2::lfcShrink(dds, coef = coef_name, type = "normal")
    } else {
      DESeq2::results(dds, contrast = c(group_var, contrast_group, control_group))
    }
    
    res_df <- as.data.frame(res) %>%
      rownames_to_column("TaxaID") %>%
      mutate(
        method = "DESeq2",
        adjusted_p_value = padj,
        score = log2FoldChange,
        significant = adjusted_p_value < thresholds$DESeq2
      ) %>%
      bind_cols(as.data.frame(tax_table(physeq_obj))[.$TaxaID, , drop = FALSE])
    
    write.csv(
      res_df,
      file.path(output_dir, sprintf("DESeq2_%s_vs_%s_%s.csv", contrast_group, control_group, filter_type)),
      row.names = FALSE
    )
    result_list$DESeq2 <- res_df
    if (verbose) message("✅ DESeq2 completed.")
    
    if (plot && is_plot_valid(res_df)) {
      p <- ggplot(res_df, aes(x = score, y = -log10(adjusted_p_value))) +
        geom_point(aes(color = significant), alpha = 0.6) +
        geom_hline(yintercept = -log10(thresholds$DESeq2), linetype = "dashed") +
        geom_text_repel(data = filter(res_df, significant), aes(label = TaxaID), size = 3, max.overlaps = 15) +
        scale_color_manual(values = c("black", "red")) +
        labs(
          title = sprintf("DESeq2: %s vs %s - %s", contrast_group, control_group, filter_type),
          x = "Log2 Fold Change", y = "-log10 Adjusted P"
        ) +
        theme_minimal()
      
      ggsave(
        file.path(output_dir, sprintf("DESeq2_volcano_%s_vs_%s_%s.png", contrast_group, control_group, filter_type)),
        plot = p, width = 6, height = 4
      )
    }
  }
  
  # ── ALDEx2 ──
  if ("ALDEx2" %in% method_list) {
    counts <- as.data.frame(otu_table(physeq_obj))
    if (!taxa_are_rows(physeq_obj)) counts <- t(counts)
    meta_df <- data.frame(sample_data(physeq_obj))
    condition <- as.character(meta_df[[group_var]])
    
    if (any(table(condition) < 2)) {
      warning(sprintf("⛔ Skipping ALDEx2: not enough replicates for comparison (%s vs %s)", control_group, contrast_group))
    } else {
      aldex_res <- ALDEx2::aldex(counts, conditions = condition, test = "t", effect = TRUE, denom = "all")
      
      aldex_df <- aldex_res %>%
        rownames_to_column("TaxaID") %>%
        mutate(
          method = "ALDEx2",
          adjusted_p_value = we.eBH,
          score = effect,
          significant = adjusted_p_value < thresholds$ALDEx2
        ) %>%
        bind_cols(as.data.frame(tax_table(physeq_obj))[.$TaxaID, , drop = FALSE])
      
      write.csv(
        aldex_df,
        file.path(output_dir, sprintf("ALDEx2_%s_vs_%s_%s.csv", contrast_group, control_group, filter_type)),
        row.names = FALSE
      )
      result_list$ALDEx2 <- aldex_df
      if (verbose) message("✅ ALDEx2 completed.")
      
      if (plot && is_plot_valid(aldex_df)) {
        p <- ggplot(aldex_df, aes(x = score, y = -log10(adjusted_p_value))) +
          geom_point(aes(color = significant), alpha = 0.6) +
          geom_text_repel(data = filter(aldex_df, significant), aes(label = TaxaID), size = 3, max.overlaps = 15) +
          geom_hline(yintercept = -log10(thresholds$ALDEx2), linetype = "dashed") +
          scale_color_manual(values = c("black", "red")) +
          labs(
            title = sprintf("ALDEx2: %s vs %s - %s", contrast_group, control_group, filter_type),
            x = "Effect Size", y = "-log10 Adjusted P"
          ) +
          theme_minimal()
        
        ggsave(
          file.path(output_dir, sprintf("ALDEx2_volcano_%s_vs_%s_%s.png", contrast_group, control_group, filter_type)),
          plot = p, width = 6, height = 4
        )
      } else if (verbose) {
        message("⚠️ Skipping ALDEx2 volcano plot — no valid points to plot.")
      }
    }
  }
  
  # ── LEfSe ──
  if ("LEfSe" %in% method_list) {
    otu <- as(otu_table(physeq_obj), "matrix")
    if (!taxa_are_rows(physeq_obj)) otu <- t(otu)
    tax <- as.data.frame(tax_table(physeq_obj))
    meta_df <- data.frame(sample_data(physeq_obj))
    
    se_raw <- SummarizedExperiment::SummarizedExperiment(
      assays = list(OTU = otu),
      rowData = tax,
      colData = meta_df
    )
    
    se_rel <- lefser::relativeAb(se_raw)
    group_labels <- meta_df[[group_var]]
    relab_mat <- assay(se_rel, "rel_abs")
    
    constant_taxa <- has_constant_taxa(relab_mat, group_labels)
    
    if (all(constant_taxa)) {
      message(sprintf("⛔ Skipping LEfSe: all taxa are constant in %s vs %s (%s)", contrast_group, control_group, filter_type))
      result_list$LEfSe <- tibble()
    } else {
      relab_mat <- relab_mat[!constant_taxa, , drop = FALSE]
      se_input <- SummarizedExperiment::SummarizedExperiment(
        assays = list(relab = relab_mat),
        rowData = rowData(se_rel)[!constant_taxa, , drop = FALSE],
        colData = colData(se_rel)
      )
      
      lefse_res <- tryCatch({
        lefser::lefser(
          se_input,
          groupCol = group_var,
          kruskal.threshold = thresholds$LEfSe$kruskal,
          lda.threshold = thresholds$LEfSe$lda,
          method = "fdr",
          assay = "relab"
        )
      }, error = function(e) {
        message(sprintf("❌ LEfSe failed for %s vs %s (%s): %s", contrast_group, control_group, filter_type, e$message))
        return(NULL)
      })
      
      if (!is.null(lefse_res)) {
        kw_pvals <- apply(relab_mat, 1, function(x) {
          kruskal.test(x ~ group_labels)$p.value
        })
        
        volcano_df <- as_tibble(lefse_res) %>%
          mutate(
            adjusted_p_value = p.adjust(kw_pvals[Names], method = "fdr"),
            method = "LEfSe",
            score = scores,
            TaxaID = Names,
            significant = adjusted_p_value < thresholds$LEfSe$kruskal & abs(score) >= thresholds$LEfSe$lda
          ) %>%
          left_join(tax %>% rownames_to_column("TaxaID"), by = "TaxaID") %>%
          relocate(TaxaID, method, score, adjusted_p_value, significant) %>%
          filter(abs(score) >= thresholds$LEfSe$lda)
        
        write.csv(
          volcano_df,
          file.path(output_dir, sprintf("LEfSe_%s_vs_%s_%s.csv", contrast_group, control_group, filter_type)),
          row.names = FALSE
        )
        result_list$LEfSe <- volcano_df
        if (verbose) message("✅ LEfSe completed.")
        
        if (plot && is_plot_valid(volcano_df)) {
          p <- ggplot(volcano_df, aes(x = score, y = -log10(adjusted_p_value))) +
            geom_point(aes(color = significant), alpha = 0.6) +
            geom_hline(yintercept = -log10(thresholds$LEfSe$kruskal), linetype = "dashed") +
            geom_vline(xintercept = c(-thresholds$LEfSe$lda, thresholds$LEfSe$lda), linetype = "dotted", color = "gray") +
            geom_text_repel(data = filter(volcano_df, significant), aes(label = TaxaID), size = 3, max.overlaps = 15) +
            scale_color_manual(values = c("black", "blue")) +
            labs(
              title = sprintf("LEfSe: %s vs %s - %s", contrast_group, control_group, filter_type),
              x = "LDA Score", y = "-log10 KW P-value"
            ) +
            theme_minimal()
          
          ggsave(
            file.path(output_dir, sprintf("LEfSe_volcano_%s_vs_%s_%s.png", contrast_group, control_group, filter_type)),
            plot = p, width = 6, height = 4
          )
        } else if (verbose) {
          message(sprintf("⚠️ Skipping volcano plot for LEfSe: %s vs %s (%s) — no valid data points.", contrast_group, control_group, filter_type))
        }
      } else {
        result_list$LEfSe <- tibble()
      }
    }
  }
  
  return(bind_rows(result_list, .id = "source"))
}