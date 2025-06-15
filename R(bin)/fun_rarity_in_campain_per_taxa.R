# ' Rarity classification of MSP
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_taxo_overlap_rarity_fullsites <- function(data_and_meta_clean_fullsites){
  
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  # Load necessary libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  data <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  data <- subset(data, meta$island == "Reunion")
  meta <- subset(meta, meta$island == "Reunion")
  
  data <- data[,!colSums(data) == 0]
  
  msp_list <- names(data) 
  msp_list_filter <- msp_list[!msp_list %in% c("Bivalvia",
                                               "Calcareous_worm_tubes",
                                               "Other_foraminifera",
                                               "Cirripedia",
                                               "Soft_worm_tubes",
                                               "X_SPON",
                                               "BIOFR",
                                               "BIOFV",
                                               "CYANOB",
                                               "No_Recruitment",
                                               "Sediments",
                                               "Encrusting_Phaeophyceae_algae",
                                               "Erect_Phaeophyceae_algae",
                                               "CCA",
                                               "Encrusting_Chlorophyta_algae",
                                               "Erect_Chlorophyta_algae",
                                               "Erect_Rhodophyta_algae")]
  
  data_filtered <- data[, msp_list_filter]
  
  
  species_names <- data.frame(Species = colnames(data_filtered))
  
  pa_matrix <- vegan::decostand(data_filtered, "pa")
  
  
  corr_taxa <- species_names %>%
    mutate(
      taxa = case_when(
        grepl("bry", Species, ignore.case = TRUE) ~ "Bryozoa",
        grepl("SPON", Species, ignore.case = TRUE) ~ "Porifera",
        grepl("ASCC|ASCS", Species, ignore.case = TRUE) ~ "Ascidiacea",
        grepl("FOR", Species, ignore.case = TRUE) ~ "Foraminifera",
        grepl("BIV|pina", Species, ignore.case = TRUE) ~ "Bivalvia",
        grepl("algae|CCA", Species, ignore.case = TRUE) ~ "Algae",
        grepl("BIOF|CYAN", Species, ignore.case = TRUE) ~ "Prokaryotic biotas",
        grepl("WORM", Species, ignore.case = TRUE) ~ "Annelida",
        grepl("HYD|Cni", Species, ignore.case = TRUE) ~ "Cnidaria",
        grepl("Cirr", Species, ignore.case = TRUE) ~ "Cirripedia",
        TRUE ~ "Other"  # Pour toutes les autres espèces
      )
    )
  

  
  # Step 1: Filter to relevant taxa
  target_taxa <- c("Bryozoa", "Porifera", "Ascidiacea")
  filtered_species <- corr_taxa %>% 
    filter(taxa %in% target_taxa)
  
  # Step 2: Reshape PA matrix to long format and merge with metadata and taxa
  pa_long <- pa_matrix %>%
    rownames_to_column("arms") %>%
    pivot_longer(-arms, names_to = "Species", values_to = "Presence") %>%
    filter(Presence == 1) %>%
    left_join(meta, by = "arms")
  
  # Step 3: Keep only target taxa species
  pa_long <- pa_long %>%
    inner_join(filtered_species, by = "Species")
  
  # Step 4: Compute frequency and rarity status
  rarity_table <- pa_long %>%
    group_by(campain, taxa, Species) %>%
    summarise(freq = n_distinct(arms), .groups = "drop") %>%
    left_join(meta %>% count(campain, name = "total_samples"), by = "campain") %>%
    mutate(freq_prop = freq / total_samples,
           status = ifelse(freq_prop < 0.25, "Rare", "Not rare")) %>%
    select(campain, Species, taxa, status)
  
  # Step 5: Pivot wider on campaigns
  rarity_wide <- rarity_table %>%
    pivot_wider(names_from = campain, values_from = status)
  
  rarity_wide <- rarity_wide %>%
    mutate(
      `Status in shallow reefs` = case_when(
        !is.na(RUNARMS) & is.na(P50ARMS) ~ "Exclusive",
        TRUE ~ RUNARMS
      ),
      `Status in mesophotic reefs` = case_when(
        is.na(RUNARMS) & !is.na(P50ARMS) ~ "Exclusive",
        TRUE ~ P50ARMS
      )
    )
  
  # Step 7: Now we still have Species and taxa information -> ready for count
  summary_table <- rarity_wide %>%
    count(`Status in shallow reefs`, `Status in mesophotic reefs`, taxa) %>%
    pivot_wider(names_from = taxa, values_from = n, values_fill = 0) %>%
    mutate(`Total number of species` = rowSums(select(., any_of(target_taxa))))
  
  # Step 8: Clean final output
  colnames(summary_table)[1:2] <- c("Status in shallow reefs", "Status in mesophotic reefs")
  summary_table <- summary_table %>%
    select(
      `Status in shallow reefs`,
      `Status in mesophotic reefs`,
      `Ascidiacea` = Ascidiacea,
      `Porifera` = Porifera,
      `Bryozoa` = Bryozoa,
      `Total number of species`
    )
  
  # View the result
  print(summary_table)
  
  return(NULL)
  
}  