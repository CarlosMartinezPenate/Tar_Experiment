flag_significance <- function(df, thresholds) {
  df %>%
    mutate(method = tolower(source)) %>%
    mutate(is_significant = FALSE) %>%
    {
      if ("padj" %in% names(.)) {
        . <- mutate(., is_significant = ifelse(method == "deseq2", !is.na(padj) & padj < thresholds$DESeq2, is_significant))
      }
      if ("we.eBH" %in% names(.)) {
        . <- mutate(., is_significant = ifelse(method == "aldex2", !is.na(we.eBH) & we.eBH < thresholds$ALDEx2, is_significant))
      }
      mutate(., is_significant = ifelse(method == "lefse", !is.na(adjusted_p_value) & adjusted_p_value < thresholds$LEfSe$kruskal, is_significant))
    }
}
