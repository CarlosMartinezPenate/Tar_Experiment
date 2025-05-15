# ───────────────────────────────────────────────────────────
# 📦 filter_phyloseq_for_comparison.R
# ───────────────────────────────────────────────────────────

filter_phyloseq_for_comparison <- function(physeq, control_group, contrast_group,
                                           group_var = "Station_treatment",
                                           super_origin_var = "Super_Origin",
                                           super_origin_control = "Tar_Experiment",
                                           super_origin_contrast = "Pollution_Events",
                                           filter_type = c("sterivex", "11_um"),
                                           filter_type_var = "Filter_Type",
                                           min_asv_count = 120) {
  message("🔍 [filter_phyloseq_for_comparison] Running filtering for: ", control_group, " vs ", contrast_group, " (", filter_type, ")")
  
  # Extract sample metadata
  meta <- data.frame(sample_data(physeq))
  
  # Subset metadata to match comparison conditions
  subset_samples <- rownames(meta[
    meta[[group_var]] %in% c(control_group, contrast_group) &
      meta[[super_origin_var]] %in% c(super_origin_control, super_origin_contrast) &
      meta[[filter_type_var]] == filter_type,
    , drop = FALSE
  ])
  
  # Filter the phyloseq object
  physeq_sub <- prune_samples(subset_samples, physeq)
  
  # Debug — print structure before matrix coercion
  message("🧼 [DEBUG] Sample data class: ", class(sample_data(physeq_sub)))
  message("🧼 [DEBUG] OTU table class: ", class(otu_table(physeq_sub)))
  
  # Coerce OTU table to matrix safely, check counts
  otu_mat <- as(otu_table(physeq_sub), "matrix")
  print(dim(otu_mat))
  
  # Remove low abundance ASVs (sum across all samples <= threshold)
  keep_taxa <- rowSums(otu_mat) > min_asv_count
  otu_mat_filtered <- otu_mat[keep_taxa, , drop = FALSE]
  
  # Debug — print ASV filter result
  message("🧼 [DEBUG] Remaining ASVs: ", sum(keep_taxa), "/", length(keep_taxa))
  
  # Rebuild phyloseq object safely
  otu_tab <- otu_table(otu_mat_filtered, taxa_are_rows = TRUE)
  new_physeq <- phyloseq::phyloseq(
    otu_tab,
    sample_data(physeq_sub),
    tax_table(physeq_sub)
  )
  
  return(new_physeq)
}