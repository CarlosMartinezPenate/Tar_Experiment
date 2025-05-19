# ───────────────────────────────────────────────────────
# 🧠 Tar Experiment: Launch All Super Comparisons
# ───────────────────────────────────────────────────────

# 📦 Load Required Libraries
suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(DESeq2)
  library(ALDEx2)
  library(lefser)
  library(pheatmap)
  library(ggrepel)
  library(ggplot2)
  library(pals)
  library(SummarizedExperiment)
  library(heatmaply)
  library(RColorBrewer)
  library(readr)
  library(htmlwidgets)
  library(stringr)
})

message("📦 Required libraries loaded.")

# 📁 Load Phyloseq Object
phy_path <- "/Users/carlosmartinez/Documents/Shewanella/contamination/Decontaminated_Prune_TarExperiment.rds"
message("🔍 Checking phyloseq path: ", phy_path)

if (!file.exists(phy_path)) {
  stop("❌ Phyloseq object not found at: ", phy_path)
}

physeq <- readRDS(phy_path)
message("✅ Phyloseq object loaded.")
print(physeq)

# 🧩 Source All Modular Functions
script_path <- "/Users/carlosmartinez/Documents/GitHub/Tar_Experiment/scripts/Modular scrips"
message("📂 Sourcing modular scripts from: ", script_path)

modular_scripts <- list.files(
  script_path,
  pattern = "\\.R$",
  full.names = TRUE
)
modular_scripts <- modular_scripts[!grepl("launch_super_comparisons\\.R$", modular_scripts)]

message("🔗 Found ", length(modular_scripts), " modular scripts to load.")
invisible(lapply(modular_scripts, function(f) {
  message("📄 Sourcing: ", basename(f))
  source(f)
}))

message("✅ All modular functions loaded.")

# ⚙️ Set Super Comparison Config
control_groups <- c("control")
comparison_groups <- c("S-2", "S-3", "DOR1", "Dor.1", "Dor.2", "Dor.3", "Dor-in", "Dor-out", "Dor.neg", "pollution_cruise")
filter_types <- c("sterivex", "11_um")

thresholds <- list(
  DESeq2 = 0.05,
  ALDEx2 = 0.25,
  LEfSe = list(kruskal = 0.15, lda = 3.5)
)

output_root <- "super_comparisons"
transform_type <- "log"

# 🚀 Run all comparisons
message("📊 Running all super comparisons...")
run_super_comparisons(
  physeq = physeq,
  control_groups = control_groups,
  comparison_groups = comparison_groups,
  filter_types = filter_types,
  thresholds = thresholds,
  methods = c("DESeq2", "ALDEx2", "LEfSe"),
  tax_level = "Genus",
  transform_type = transform_type,
  output_root = output_root
)

message("🏁 All comparisons complete.")

# ────────────────────────────────────────────────
# 🌋 Generate volcano plots once, post comparisons
# ────────────────────────────────────────────────

message("🧪 Generating DESeq2-only volcano plots...")
generate_deseq2_only_volcanoes_pretty(
  results_dir = output_root,
  physeq = physeq,
  group_var = "Station_treatment",
  tax_level = "Genus"
)

message("🌋 Generating combined volcano plots...")
generate_all_combined_volcanoes(
  results_dir = output_root,
  physeq = physeq,
  group_var = "Station_treatment",
  tax_level = "Genus"
)

message("✅ All robust volcano plots generated.")
# ───────────────────────────────────────────────────────
# 🔥 Generate Robust Heatmaps for DESeq2 Volcano Results
# ───────────────────────────────────────────────────────

message("🧬 Proceeding to generate robust heatmaps...")

volcano_dir <- file.path(output_root, "deseq2_only_volcano_plots")
heatmap_out_dir <- file.path(output_root, "heatmap_outputs")
heatmap_html_dir <- file.path(output_root, "interactive_heatmaps")
transform_methods <- c("zscore", "clr", "ra", "log2", "arcsine", "vst")

if (!dir.exists(volcano_dir)) {
  stop("❌ Volcano directory not found: ", volcano_dir)
}

volcano_files <- list.files(volcano_dir, pattern = "^deseq2_data_.*\\.csv$", full.names = TRUE)
if (length(volcano_files) == 0) {
  stop("❌ No DESeq2 volcano CSV files found in: ", volcano_dir)
}

# 🔁 Loop over transformation methods and volcano files
for (method in transform_methods) {
  message("\n🔁 Running heatmap transformation: ", method)
  
  for (f in volcano_files) {
    fname <- basename(f)
    
    # ⛏️ Extract metadata from file name
    match <- stringr::str_match(fname, "^deseq2_data_(.*)_vs_(.*)_(.*)\\.csv$")
    if (any(is.na(match))) {
      message("⚠️ Skipping malformed filename: ", fname)
      next
    }
    
    contrast_group <- match[2]
    control_group <- match[3]
    filter_type <- match[4]
    
    tryCatch({
      run_heatmap_from_volcano(
        physeq_obj = physeq,
        volcano_csv_path = f,
        transformation = method,
        output_dir = heatmap_out_dir,
        html_dir = heatmap_html_dir,
        tax_level = "Genus",
        filter_type = filter_type  # ✅ Fix was here
      )
    }, error = function(e) {
      message("❌ Failed heatmap: ", fname, " | ", e$message)
    })
  }
}

message("✅ All robust heatmaps generated.")