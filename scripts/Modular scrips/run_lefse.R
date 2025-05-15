run_lefse <- function(physeq, group_var, kruskal_threshold, lda_threshold) {
  otu <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) otu <- t(otu)
  
  tax <- as.data.frame(tax_table(physeq))
  meta <- data.frame(sample_data(physeq))
  
  se_raw <- SummarizedExperiment(assays = list(OTU = otu), rowData = tax, colData = meta)
  se_rel <- lefser::relativeAb(se_raw)
  se_input <- SummarizedExperiment(assays = list(relab = assay(se_rel, "rel_abs")), rowData = tax, colData = meta)
  
  lefse_res <- lefser::lefser(
    se_input,
    classCol = group_var,
    kruskal.threshold = kruskal_threshold,
    lda.threshold = lda_threshold,
    method = "fdr",
    assay = "relab"
  )
  
  kw_pvals <- apply(assay(se_input, "relab"), 1, function(x) {
    kruskal.test(x ~ colData(se_input)[[group_var]])$p.value
  })
  
  as_tibble(lefse_res) %>%
    mutate(
      adjusted_p_value = p.adjust(kw_pvals[Names], method = "fdr"),
      method = "LEfSe",
      TaxaID = Names,
      score = scores,
      is_significant = adjusted_p_value < kruskal_threshold & abs(score) >= lda_threshold
    )
}