summarize_asv_detection <- function(results_dir, save_dir = NULL) {
  library(dplyr)
  library(readr)
  library(stringr)
  
  message("🔍 Scanning for ASV results files in: ", results_dir)
  
  files <- list.files(results_dir, pattern = "^results_asvs_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  message("📂 Found ", length(files), " result files.")
  
  summary_list <- list()
  lost_asvs_list <- list()
  
  for (file in files) {
    fname <- basename(file)
    parts <- str_match(fname, "^results_asvs_(.*?)_vs_(.*?)_(.*)\\.csv$")
    if (any(is.na(parts))) next
    
    contrast_group <- parts[2]
    control_group <- parts[3]
    filter_type <- parts[4]
    label <- paste(contrast_group, "vs", control_group, "(", filter_type, ")")
    
    message("📄 Processing: ", label)
    df <- read_csv(file, show_col_types = FALSE)
    
    # Basic diagnostics
    cat("🧪 Preview of input data (", fname, "):\n")
    print(head(df))
    
    asv_summary <- df %>%
      group_by(TaxaID) %>%
      summarise(
        detected_in_deseq2 = any(source == "DESeq2" & is_significant),
        detected_in_aldex2 = any(source == "ALDEx2" & is_significant),
        detected_in_lefse = any(source == "LEfSe" & is_significant),
        .groups = "drop"
      ) %>%
      mutate(
        case = case_when(
          detected_in_deseq2 ~ "Detected by DESeq2",
          detected_in_aldex2 & !detected_in_deseq2 ~ "Only ALDEx2",
          detected_in_lefse & !detected_in_deseq2 ~ "Only LEfSe",
          TRUE ~ "Not significant anywhere"
        )
      )
    
    counts <- asv_summary %>%
      count(case) %>%
      mutate(contrast = label)
    
    summary_list[[label]] <- counts
    
    lost_asvs <- asv_summary %>%
      filter(case %in% c("Only ALDEx2", "Only LEfSe")) %>%
      mutate(contrast = label)
    
    lost_asvs_list[[label]] <- lost_asvs
  }
  
  full_summary <- bind_rows(summary_list)
  lost_asvs_df <- bind_rows(lost_asvs_list)
  
  # Show structure of return object
  cat("\n✅ Structure of return object:\n")
  str(list(summary = full_summary, lost_asvs = lost_asvs_df))
  
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    write_csv(full_summary, file.path(save_dir, "asv_detection_summary.csv"))
    write_csv(lost_asvs_df, file.path(save_dir, "lost_asvs_details.csv"))
    message("💾 Saved summary and lost ASVs to ", save_dir)
  }
  
  return(list(summary = full_summary, lost_asvs = lost_asvs_df))
}