#' Rarity and commonness of morphospecies
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_tableau_presence <- function(data_and_meta_clean_fullsites){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites")
  
  # Charger les packages nécessaires
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(flextable)
  library(officer)
  library(scales)
  library(writexl)
  
  data <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
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
                                               "Erect_Rhodophyta_algae",
                                               "Cni_Plumulariidae")]
  
  
  data <- data[, msp_list_filter]
  
  species_names <- data.frame(Species = colnames(data))
  
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
  
  
  
  # 1. Ajouter les infos de site et triplicat depuis `meta`
  data$arms <- rownames(data)
  data_joined <- data %>%
    left_join(meta, by = "arms")
  
  # 2. Convertir en format long
  data_long <- data_joined %>%
    pivot_longer(
      cols = -c(arms, campain, triplicat, site, island),
      names_to = "taxon",
      values_to = "abundance"
    )
  
  # 3. Créer une colonne de présence/absence
  data_long <- data_long %>%
    mutate(presence = ifelse(abundance > 0, 1, 0))
  
  # 4. Résumer : présence/absence par taxon et triplicat
  presence_table <- data_long %>%
    group_by(taxon, triplicat) %>%
    summarise(presence = ifelse(sum(presence) > 0, "✓", ""), .groups = "drop")
  
  # 5. Récupérer l’ordre des triplicats triés par campagne
  triplicat_order <- meta %>%
    select(triplicat, campain) %>%
    distinct() %>%
    arrange(campain, triplicat) %>%
    pull(triplicat) %>%
    unique()
  
  # 6. Restructurer en wide format avec colonnes triées
  presence_wide <- presence_table %>%
    pivot_wider(names_from = triplicat, values_from = presence) %>%
    select(taxon, all_of(triplicat_order))
  
  # 7. Créer la ligne "campain" (en-tête)
  campain_row <- meta %>%
    select(triplicat, campain) %>%
    distinct() %>%
    filter(triplicat %in% triplicat_order) %>%
    arrange(factor(triplicat, levels = triplicat_order)) %>%
    pull(campain)
  
  # 8. Ajouter cette ligne comme première ligne
  header_row <- c("campain", campain_row)
  final_table <- rbind(header_row, presence_wide)
  names <- colnames(final_table)
  names[1] <- "Species"
  colnames(final_table) <- names
  
  final_table <- left_join(corr_taxa, final_table, by = "Species")
  
  
  
  # Exporter le data frame
  write_xlsx(final_table, "outputs/species_presence_sites.xlsx")
  
  
  final_table <- final_table %>%
    arrange(taxa)
  
  ft <- flextable(final_table)
  ft <- autofit(ft)
  ft <- align(ft, align = "center", part = "all")
  ft <- set_table_properties(ft, layout = "autofit")
  
  # Créer le document Word et y insérer le tableau
  doc <- read_docx() %>%
    body_add_par("Species Rarity Summary Table", style = "heading 1") %>%
    body_add_flextable(ft)
  
  # Exporter le fichier Word
  print(doc, target = "outputs/species_presence_sites.docx")
  
  # Renommer les c
  
  # 9. Afficher le tableau final
  print(final_table)

  return(NULL)
    
}

