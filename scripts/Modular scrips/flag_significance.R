flag_significance <- function(df, thresholds) {
  cat("🔍 Starting flag_significance() (base R version)...\n")
  cat("📐 Input dataframe dimensions:\n")
  print(dim(df))
  
  # Ensure method column is lowercase
  df$method <- tolower(df$source)
  
  # Initialize flags
  df$significant_in_deseq2 <- FALSE
  df$significant_in_aldex2 <- FALSE
  df$significant_in_lefse  <- FALSE
  
  # Flag DESeq2
  if ("padj" %in% colnames(df)) {
    df$significant_in_deseq2 <- with(df, method == "deseq2" & !is.na(padj) & padj < thresholds$DESeq2)
  }
  
  # Flag ALDEx2
  if ("we.eBH" %in% colnames(df)) {
    df$significant_in_aldex2 <- with(df, method == "aldex2" & !is.na(we.eBH) & we.eBH < thresholds$ALDEx2)
  }
  
  # Flag LEfSe
  if ("adjusted_p_value" %in% colnames(df)) {
    df$significant_in_lefse <- with(df, method == "lefse" & !is.na(adjusted_p_value) & adjusted_p_value < thresholds$LEfSe$kruskal)
  }
  
  # Summarize by TaxaID
  cat("📊 Summarizing significance by TaxaID...\n")
  df_summary <- aggregate(
    cbind(significant_in_deseq2, significant_in_aldex2, significant_in_lefse) ~ TaxaID,
    data = df,
    FUN = any
  )
  
  # Add count of methods
  df_summary$n_methods <- rowSums(df_summary[, c("significant_in_deseq2", "significant_in_aldex2", "significant_in_lefse")])
  df_summary$is_significant <- df_summary$n_methods > 0
  
  # Merge with original df
  df <- merge(df, df_summary, by = "TaxaID", all.x = TRUE)
  
  cat("✅ Finished flag_significance(). Output dimensions:\n")
  print(dim(df))
  
  # Print preview if columns are available
  preview_cols <- c("TaxaID", "n_methods", "is_significant")
  existing_cols <- preview_cols[preview_cols %in% colnames(df)]
  if (length(existing_cols) > 0) {
    cat("📋 Sample output (first few rows):\n")
    print(head(df[, existing_cols], 6))
  } else {
    cat("⚠️ No preview columns available in final dataframe.\n")
  }
  
  return(df)
}