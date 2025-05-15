run_deseq2 <- function(physeq, group_var, control_group, contrast_group, threshold) {
  dds <- phyloseq_to_deseq2(physeq, as.formula(paste("~", group_var)))
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)
  coef_name <- paste0(group_var, "_", contrast_group, "_vs_", control_group)
  
  res <- if (coef_name %in% resultsNames(dds)) {
    lfcShrink(dds, coef = coef_name, type = "normal")
  } else {
    results(dds, contrast = c(group_var, contrast_group, control_group))
  }
  
  as.data.frame(res) %>%
    rownames_to_column("TaxaID") %>%
    mutate(
      method = "DESeq2",
      adjusted_p_value = padj,
      score = log2FoldChange,
      is_significant = !is.na(adjusted_p_value) & adjusted_p_value < threshold
    )
}