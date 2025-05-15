# safe_subset_phyloseq.R
safe_subset_phyloseq <- function(physeq, group_var, control_group, contrast_group, filter_type = NULL, min_total_reads = 120) {
  cat("🔍 Entering safe_subset_phyloseq()\n")
  
  meta <- data.frame(sample_data(physeq))
  cat("📋 Metadata dimensions before filter:", dim(meta), "\n")
  
  # Optional filter
  if (!is.null(filter_type)) {
    if (!"Filter_Type" %in% colnames(meta)) {
      stop("❌ 'Filter_Type' column not found in sample_data()")
    }
    meta <- meta[meta$Filter_Type == filter_type, ]
    cat("📋 Metadata dimensions after filter_type =", filter_type, ":", dim(meta), "\n")
  }
  
  keep_samples <- rownames(meta)[meta[[group_var]] %in% c(control_group, contrast_group)]
  
  if (length(keep_samples) == 0) {
    warning("🚫 No matching samples for contrast/control groups after filtering.")
    return(NULL)
  }
  
  physeq_sub <- prune_samples(keep_samples, physeq)
  cat("🧪 Samples after pruning:", nsamples(physeq_sub), "\n")
  
  physeq_sub <- prune_taxa(taxa_sums(physeq_sub) > min_total_reads, physeq_sub)
  cat("🧬 Taxa after filtering for min_total_reads =", min_total_reads, ":", ntaxa(physeq_sub), "\n")
  
  if (nsamples(physeq_sub) == 0 || ntaxa(physeq_sub) == 0) {
    warning(sprintf("🚫 Skipping %s vs %s (%s): zero samples or taxa", contrast_group, control_group, filter_type))
    return(NULL)
  }
  
  meta_sub <- data.frame(sample_data(physeq_sub))
  if (length(unique(meta_sub[[group_var]])) < 2 || nsamples(physeq_sub) < 4) {
    warning(sprintf("🚫 Skipping %s vs %s (%s): not enough samples or diversity in group_var", contrast_group, control_group, filter_type))
    return(NULL)
  }
  
  cat("✅ Subset complete. Returning phyloseq of class:", class(physeq_sub), "\n")
  return(physeq_sub)
}