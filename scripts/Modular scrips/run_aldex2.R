run_aldex2 <- function(physeq, group_var, control_group, contrast_group, threshold) {
  counts <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) counts <- t(counts)
  
  meta <- data.frame(sample_data(physeq))
  cond <- meta[[group_var]]
  stopifnot(ncol(counts) == length(cond))
  
  aldex_res <- ALDEx2::aldex(counts, conditions = cond, test = "t", effect = TRUE, denom = "all")
  
  aldex_res %>%
    rownames_to_column("TaxaID") %>%
    mutate(
      method = "ALDEx2",
      adjusted_p_value = we.eBH,
      score = effect,
      is_significant = !is.na(adjusted_p_value) & adjusted_p_value < threshold
    )
}