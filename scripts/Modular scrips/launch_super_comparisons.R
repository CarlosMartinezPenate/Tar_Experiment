# ───────────────────────────────────────────────────────
# 🧠 Tar Experiment: Launch All Super Comparisons
# ───────────────────────────────────────────────────────

# 📦 Load Required Libraries
suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
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
message("⚙️ Setting control and comparison groups...")
control_groups <- c("control")
comparison_groups <- c("S-2", "S-3", "DOR1", "Dor.1", "Dor.2", "Dor.3", "Dor-in", "Dor-out", "Dor.neg", "pollution_cruise")
filter_types <- c("sterivex", "11_um")

message("⚙️ Setting method thresholds...")
thresholds <- list(
  DESeq2 = 0.05,
  ALDEx2 = 0.25,
  LEfSe = list(kruskal = 0.15, lda = 3.5)
)

output_root <- "super_comparisons"
transform_type <- "log"  # Options: log, clr, raw

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
