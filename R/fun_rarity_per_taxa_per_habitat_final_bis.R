# ' Rarity classification of MSP
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_rarity_per_taxa_per_habitat_fullsites_bis <- function(data_and_meta_clean_fullsites){
  
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  #### Step 1 : Load necessary libraries and data ####
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(flextable)
  library(officer)
  library(scales)
  
  # data - plate faces
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)
  data <- subset(data, meta$island == "Reunion")
  meta <- subset(meta, meta$island == "Reunion")
  
  data <- data[,!colSums(data) == 0]
  
  #### Step 2 : Filtrer les taxons d'interet ####
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
  
  data_filtered_pa <- vegan::decostand(data_filtered, "pa")
  
  
  # data - ARMS
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
  data_mean <- data_mean[,!colSums(data_mean) == 0]
  
  data_mean_filtered <- data_mean[, msp_list_filter]
  data_mean_filtered_pa <- vegan::decostand(data_mean_filtered, "pa")
  
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
  
  
  #### Step 3: Separer par campagne ####
  p50arms_presence_mean <- subset(data_mean_filtered_pa, meta_mean$campain == "P50ARMS")
  runarms_presence_mean <- subset(data_mean_filtered_pa, meta_mean$campain == "RUNARMS")
  
  p50arms_meta_mean <- subset(meta_mean, meta_mean$campain == "P50ARMS")
  runarms_meta_mean <- subset(meta_mean, meta_mean$campain == "RUNARMS")
  
  #### Step 3bis: Separer par campagne ####
  p50arms_presence <- subset(data_filtered_pa, meta$campain == "P50ARMS")
  runarms_presence <- subset(data_filtered_pa, meta$campain == "RUNARMS")
  
  p50arms_meta <- subset(meta, meta$campain == "P50ARMS")
  runarms_meta <- subset(meta, meta$campain == "RUNARMS")
  
  #### Step 4: Compute Rarity threshold taxa par taxa ####  
  
  #Ascidians - p50a
  ascidiacea_species <- corr_taxa$Species[corr_taxa$taxa == "Ascidiacea"]
  data_ascidiacea_p50a <- p50arms_presence[, colnames(p50arms_presence) %in% ascidiacea_species]
  
  species_occurence <- colSums(data_ascidiacea_p50a)/nrow(data_ascidiacea_p50a)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_asc_p50a <- quantile(species_occurence$f, 0.25)
  
  #Ascidian - runa
  data_ascidiacea_runa <- runarms_presence[, colnames(runarms_presence) %in% ascidiacea_species]
  
  species_occurence <- colSums(data_ascidiacea_runa)/nrow(data_ascidiacea_runa)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_asc_runa <- quantile(species_occurence$f, 0.25)
  
  
  #Porifera
  porifera_species <- corr_taxa$Species[corr_taxa$taxa == "Porifera"]
  data_porifera_p50a <- p50arms_presence[, colnames(p50arms_presence) %in% porifera_species]
  
  species_occurence <- colSums(data_porifera_p50a)/nrow(data_porifera_p50a)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_por_p50a <- quantile(species_occurence$f, 0.25)
  
  #Porifera - runa
  data_porifera_runa <- runarms_presence[, colnames(runarms_presence) %in% porifera_species]
  
  species_occurence <- colSums(data_porifera_runa)/nrow(data_porifera_runa)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_por_runa <- quantile(species_occurence$f, 0.25)
  
  #Bryozoa
  bryozoa_species <- corr_taxa$Species[corr_taxa$taxa == "Bryozoa"]
  data_bryozoa_p50a <- p50arms_presence[, colnames(p50arms_presence) %in% bryozoa_species]
  
  species_occurence <- colSums(data_bryozoa_p50a)/nrow(data_bryozoa_p50a)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_bry_p50a <- quantile(species_occurence$f, 0.25)
  
  #Bryozoa - runa
  data_bryozoa_runa <- runarms_presence[, colnames(runarms_presence) %in% bryozoa_species]
  
  species_occurence <- colSums(data_bryozoa_runa)/nrow(data_bryozoa_runa)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_bry_runa <- quantile(species_occurence$f, 0.25)
  
  #all_species - p50a
  species_occurence <- colSums(p50arms_presence)/nrow(p50arms_presence)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_all_species_p50a <- quantile(species_occurence$f, 0.25)
  
  #all_species - runa
  species_occurence <- colSums(runarms_presence)/nrow(runarms_presence)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_all_species_runa <- quantile(species_occurence$f, 0.25)
  
  #### Step 4bis: Compute Rarity threshold taxa par taxa ####  
  
  #Ascidians - p50a
  data_mean_ascidiacea_p50a <- p50arms_presence_mean[, colnames(p50arms_presence_mean) %in% ascidiacea_species]
  
  species_occurence_mean <- colSums(data_mean_ascidiacea_p50a)/nrow(data_mean_ascidiacea_p50a)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                  Species = names(species_occurence_mean))
  threshold_asc_p50a_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #Ascidian - runa
  data_mean_ascidiacea_runa <- runarms_presence_mean[, colnames(runarms_presence_mean) %in% ascidiacea_species]
  
  species_occurence_mean <- colSums(data_mean_ascidiacea_runa)/nrow(data_mean_ascidiacea_runa)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_asc_runa_mean <- quantile(species_occurence_mean$f, 0.25)
  
  
  #Porifera - p50a
  data_mean_porifera_p50a <- p50arms_presence_mean[, colnames(p50arms_presence_mean) %in% porifera_species]
  
  species_occurence_mean <- colSums(data_mean_porifera_p50a)/nrow(data_mean_porifera_p50a)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                  Species = names(species_occurence_mean))
  threshold_por_p50a_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #Porifera - runa
  data_mean_porifera_runa <- runarms_presence_mean[, colnames(runarms_presence_mean) %in% porifera_species]
  
  species_occurence_mean <- colSums(data_mean_porifera_runa)/nrow(data_mean_porifera_runa)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_por_runa_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #Bryozoa - p50a
  data_mean_bryozoa_p50a <- p50arms_presence_mean[, colnames(p50arms_presence_mean) %in% bryozoa_species]
  
  species_occurence_mean <- colSums(data_mean_bryozoa_p50a)/nrow(data_mean_bryozoa_p50a)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_bry_p50a_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #Bryozoa - runa
  data_mean_bryozoa_runa <- runarms_presence_mean[, colnames(runarms_presence_mean) %in% bryozoa_species]
  
  species_occurence_mean <- colSums(data_mean_bryozoa_runa)/nrow(data_mean_bryozoa_runa)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_bry_runa_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #all_species - p50a
  species_occurence_mean <- colSums(p50arms_presence_mean)/nrow(p50arms_presence_mean)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                  Species = names(species_occurence_mean))
  threshold_all_species_p50a_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #all_species - runa
  species_occurence_mean <- colSums(runarms_presence_mean)/nrow(runarms_presence_mean)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                  Species = names(species_occurence_mean))
  threshold_all_species_runa_mean <- quantile(species_occurence_mean$f, 0.25)
  
  # species_list <- ascidiacea_species
  # group_name <- "Ascidiacea"
  # threshold <- threshold_asc_p50a
  # threshold_mean <- threshold_asc_p50a_mean
  
  #### Step 5: Fonctions  pour calculer la rareté pour chaque groupe ####
  compute_rarity_p50a <- function(species_list, group_name, threshold, threshold_mean) {
    df <- p50arms_presence[, colnames(p50arms_presence) %in% species_list]
    freq <- colSums(df) / nrow(df)
    
    df_mean <- p50arms_presence_mean[, colnames(p50arms_presence_mean) %in% species_list]
    freq_mean <- colSums(df_mean) / nrow(df_mean)
    
    
    df$triplicat <- p50arms_meta$triplicat
    
    result <- sapply(species_list, function(espece) {
      subset <- df[df[[espece]] == 1, ]
      n_triplicats <- length(unique(subset$triplicat))
      if (n_triplicats <= 1) {
        return("Un seul triplicat")
      } else {
        return("Plusieurs triplicats")
      }
    })
    
    out <- data.frame(Species = names(freq), f = freq, f_mean = freq_mean, n_trip = result)
    out$taxa <- group_name
    out$threshold <- threshold
    out$threshold_mean <- threshold_mean
    
    out <- out %>%
      mutate(status = case_when(
        f == 0 ~ "absent",
        f <= threshold & f_mean <= threshold_mean ~ "rare",
        TRUE ~ "common"
      ))
    
    return(out)
  }
  
  compute_rarity_runa <- function(species_list, group_name, threshold, threshold_mean) {
    df <- runarms_presence[, colnames(runarms_presence) %in% species_list]
    freq <- colSums(df) / nrow(df)
    
    df_mean <- runarms_presence_mean[, colnames(runarms_presence_mean) %in% species_list]
    freq_mean <- colSums(df_mean) / nrow(df_mean)
    
    df$triplicat <- runarms_meta$triplicat
    
    result <- sapply(species_list, function(espece) {
      subset <- df[df[[espece]] == 1, ]
      n_triplicats <- length(unique(subset$triplicat))
      if (n_triplicats <= 1) {
        return("Un seul triplicat")
      } else {
        return("Plusieurs triplicats")
      }
    })
    
    out <- data.frame(Species = names(freq), f = freq, f_mean = freq_mean, n_trip = result)
    out$taxa <- group_name
    out$threshold <- threshold
    out$threshold_mean <- threshold_mean
    
    out <- out %>%
      mutate(status = case_when(
        f == 0 ~ "absent",
        f <= threshold & f_mean <= threshold_mean ~ "rare",
        TRUE ~ "common"
      ))
    
    return(out)
  }
  
  #### Step 6: Attribuer les statuts de rareté aux espèces ####
  
  # Par groupe avec seuils spécifiques
  out_asc_runa <- compute_rarity_runa(ascidiacea_species, "Ascidiacea", threshold_asc_runa, threshold_asc_runa_mean)
  out_asc_p50a <- compute_rarity_p50a(ascidiacea_species, "Ascidiacea", threshold_asc_p50a, threshold_asc_p50a_mean)
  
  
  out_por_runa <- compute_rarity_runa(porifera_species, "Porifera", threshold_por_runa, threshold_por_runa_mean)
  out_por_p50a <- compute_rarity_p50a(porifera_species, "Porifera", threshold_por_p50a, threshold_por_p50a_mean)
  
  out_bry_runa <- compute_rarity_runa(bryozoa_species, "Bryozoa", threshold_bry_runa, threshold_bry_runa_mean)
  out_bry_p50a <- compute_rarity_p50a(bryozoa_species, "Bryozoa", threshold_bry_p50a, threshold_bry_p50a_mean)
  
  autres_groupes <- corr_taxa %>%
    filter(!taxa %in% c("Ascidiacea", "Porifera", "Bryozoa")) %>%
    pull(Species)
  
  out_others_runa <- compute_rarity_runa(autres_groupes, "Other", threshold_all_species_runa, threshold_all_species_runa_mean)
  out_others_p50a <- compute_rarity_p50a(autres_groupes, "Other", threshold_all_species_p50a, threshold_all_species_p50a_mean)
  
  recap_rarity_p50a <- bind_rows(out_asc_p50a, out_por_p50a, out_bry_p50a, out_others_p50a)
  recap_rarity_runa <- bind_rows(out_asc_runa, out_por_runa, out_bry_runa, out_others_runa)
  
  recap_rarity_p50a <- left_join(recap_rarity_p50a, corr_taxa, by = "Species") %>%
    mutate(taxa = ifelse(taxa.x == "Other", taxa.y, taxa.x)) %>%
    select(Species, taxa, f, threshold, status)
  
  recap_rarity_runa <- left_join(recap_rarity_runa, corr_taxa, by = "Species") %>%
    mutate(taxa = ifelse(taxa.x == "Other", taxa.y, taxa.x)) %>%
    select(Species, taxa, f, threshold, status)
  
  species_rarity <- data.frame(Species = recap_rarity_p50a$Species,
                               taxa = recap_rarity_p50a$taxa,
                               Status_runa = recap_rarity_runa$status, 
                               Status_p50a = recap_rarity_p50a$status)
  
  
  # First, define your categories
  species_rarity_clean <- species_rarity %>%
    mutate(
      shallow = case_when(
        Status_runa == "rare" ~ "Rare",
        Status_runa == "common" ~ "Not rare",
        Status_runa == "absent" ~ "Absent"
      ),
      mesophotic = case_when(
        Status_p50a == "rare" ~ "Rare",
        Status_p50a == "common" ~ "Not rare",
        Status_p50a == "absent" ~ "Absent"
      )
    ) %>%
    mutate(
      mesophotic = ifelse(shallow == "Absent" & mesophotic %in% c("Rare", "Not rare"), "Exclusive", mesophotic),
      shallow = ifelse(mesophotic == "Absent" & shallow %in% c("Rare", "Not rare"), "Exclusive", shallow)
    )
  
  #### Digression ####
  total_exclusive_mesophotic <- species_rarity_clean %>%
    filter(mesophotic == "Exclusive") %>%
    nrow()
  
  # Nombre d'espèces exclusives au mésophotique et "rare" dans Status_p50a
  rare_exclusive_mesophotic <- species_rarity_clean %>%
    filter(mesophotic == "Exclusive", Status_p50a == "rare") %>%
    nrow()
  
  #Pourcentage d'espèce mesophotique exclusive qui sont rare  
  percentage_mesophotic <- (rare_exclusive_mesophotic / total_exclusive_mesophotic) * 100

  total_exclusive_shallow <- species_rarity_clean %>%
    filter(shallow == "Exclusive") %>%
    nrow()
  
  # Nombre d'espèces exclusives au milieu peu profond et "rare" dans Status_runa
  rare_exclusive_shallow <- species_rarity_clean %>%
    filter(shallow == "Exclusive", Status_runa == "rare") %>%
    nrow()
  
  # Pourcentage
  percentage_shallow <- (rare_exclusive_shallow / total_exclusive_shallow) * 100
  
  

  # Étape 1 : filtrer toutes les espèces rares dans l’un ou l’autre milieu
  rare_species <- species_rarity_clean %>%
    filter(Status_runa == "rare" | Status_p50a == "rare")
  
  # Étape 2 : nombre total d'espèces rares
  total_rare_species <- nrow(rare_species)
  
  # Étape 3 : parmi elles, combien sont exclusives au peu profond
  rare_exclusive <- rare_species %>%
    filter(shallow == "Exclusive" | mesophotic == "Exclusive") %>%
    nrow()
  
  # Calcul du pourcentage d'espèce rare qui sont exclusive à un habitat ou à l'autre
  percentage_rare_exclusive <- (rare_exclusive/total_rare_species)*100 
 
  
  
  # Étape 5 : calcul des pourcentages
  percentage_exclusive_shallow <- (rare_exclusive_shallow / total_rare_species) * 100
  percentage_exclusive_mesophotic <- (rare_exclusive_mesophotic / total_rare_species) * 100
  
  
  #### Suite ####
  
  
  species_rarity_clean <- species_rarity_clean[, !names(species_rarity_clean) %in% c("Status_runa", "Status_p50a")]
  
  species_rarity_clean <- species_rarity_clean %>%
    mutate(
      taxon_group = case_when(
        taxa == "Ascidiacea" ~ "Ascidiacea",
        taxa == "Porifera" ~ "Porifera",
        taxa == "Bryozoa" ~ "Bryozoa",
        taxa %in% c("Foraminifera", "Cnidaria") ~ "Other",
        TRUE ~ "Other"  # au cas où d'autres apparaissent
      )
    )
  
  #### Step 7: décompter le nombre d'espèce rare ou pas dans les habitats ####
  
  
  # Calcul du nombre d'espèces par combinaison (shallow, mesophotic, taxon_group)
  summary_counts <- species_rarity_clean %>%
    group_by(shallow, mesophotic, taxon_group) %>%
    summarise(count = n(), .groups = "drop")
  
  # Calcul du total par groupe taxonomique (pour calculer les % par colonne)
  totals_by_group <- summary_counts %>%
    group_by(taxon_group) %>%
    summarise(total = sum(count), .groups = "drop")
  
  # Jointure et calcul des pourcentages
  summary_table <- summary_counts %>%
    left_join(totals_by_group, by = "taxon_group") %>%
    mutate(formatted = paste0(count, " (", round(100 * count / total), "%)")) %>%
    select(-count, -total) %>%
    pivot_wider(
      names_from = taxon_group,
      values_from = formatted,
      values_fill = "-"
    )
  
  # Ajout du total du nombre d'espèces pour chaque combinaison
  total_species <- species_rarity_clean %>%
    group_by(shallow, mesophotic) %>%
    summarise(Total_number_of_species = n(), .groups = "drop")
  
  # Fusion finale
  summary_table_final <- summary_table %>%
    left_join(total_species, by = c("shallow", "mesophotic")) %>%
    select(shallow, mesophotic, Ascidiacea, Porifera, Bryozoa, Other, Total_number_of_species)
  
  # Affichage
  print(summary_table_final)
  
  #### Step 8: Enregistrer le tableau ####
  
  # Reprend ton tableau résumé (déjà bien formaté avec les colonnes et valeurs type "n (x%)")
  # Ici on part du dataframe `summary_table`
  
  # Créer le flextable
  ft <- flextable(summary_table_final)
  ft <- autofit(ft)
  ft <- align(ft, align = "center", part = "all")
  ft <- set_table_properties(ft, layout = "autofit")
  
  # Créer le document Word et y insérer le tableau
  doc <- read_docx() %>%
    body_add_par("Species Rarity Summary Table", style = "heading 1") %>%
    body_add_flextable(ft)
  
  # Exporter le fichier Word
  print(doc, target = "outputs/species_rarity_summary_P50A_RUNA_04_06_2025.docx")

  # Meme comparaison entre Réunion et Rodrigues ####
  
  #### Step 1 : Load necessary libraries ####
  
  #faces
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)
  
  data <- data[meta$campain %in% c("RUNARMS", "RODARMS"), ]
  meta <- meta[meta$campain %in% c("RUNARMS", "RODARMS"), ]
  
  data <- data[,!colSums(data) == 0]
  
  #ARMS
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[meta_mean$campain %in% c("RUNARMS", "RODARMS"), ]
  meta_mean <- meta_mean[meta_mean$campain %in% c("RUNARMS", "RODARMS"), ]
  
  data_mean <- data_mean[,!colSums(data_mean) == 0]
  
  #### Step 2 : Filtrer les taxons d'interet ####
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
  
  data_filtered_pa <- vegan::decostand(data_filtered, "pa")
  
  #arms
  
  data_mean_filtered <- data_mean[, msp_list_filter]

  data_mean_filtered_pa <- vegan::decostand(data_mean_filtered, "pa")
  
  
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
  
  
  #### Step 3: Separer par campagne ####
  rodarms_presence <- subset(data_filtered_pa, meta$campain == "RODARMS")
  runarms_presence <- subset(data_filtered_pa, meta$campain == "RUNARMS")
  
  rodarms_meta <- subset(meta, meta$campain == "RODARMS")
  runarms_meta <- subset(meta, meta$campain == "RUNARMS")
  
  #### Step 3bis: Separer par campagne ####
  
  rodarms_presence_mean <- subset(data_mean_filtered_pa, meta_mean$campain == "RODARMS")
  runarms_presence_mean <- subset(data_mean_filtered_pa, meta_mean$campain == "RUNARMS")
  
  rodarms_meta_mean <- subset(meta_mean, meta_mean$campain == "RODARMS")
  runarms_meta_mean <- subset(meta_mean, meta_mean$campain == "RUNARMS")
  
  #### Step 4: Compute Rarity threshold taxa par taxa ####  
  
  #Ascidians - roda
  ascidiacea_species <- corr_taxa$Species[corr_taxa$taxa == "Ascidiacea"]
  data_ascidiacea_roda <- rodarms_presence[, colnames(rodarms_presence) %in% ascidiacea_species]
  
  species_occurence <- colSums(data_ascidiacea_roda)/nrow(data_ascidiacea_roda)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_asc_roda <- quantile(species_occurence$f, 0.25)
  
  #Ascidian - runa
  data_ascidiacea_runa <- runarms_presence[, colnames(runarms_presence) %in% ascidiacea_species]
  
  species_occurence <- colSums(data_ascidiacea_runa)/nrow(data_ascidiacea_runa)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_asc_runa <- quantile(species_occurence$f, 0.25)
  
  
  #Porifera
  porifera_species <- corr_taxa$Species[corr_taxa$taxa == "Porifera"]
  data_porifera_roda <- rodarms_presence[, colnames(rodarms_presence) %in% porifera_species]
  
  species_occurence <- colSums(data_porifera_roda)/nrow(data_porifera_roda)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_por_roda <- quantile(species_occurence$f, 0.25)
  
  #Porifera - runa
  data_porifera_runa <- runarms_presence[, colnames(runarms_presence) %in% porifera_species]
  
  species_occurence <- colSums(data_porifera_runa)/nrow(data_porifera_runa)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_por_runa <- quantile(species_occurence$f, 0.25)
  
  #Bryozoa
  bryozoa_species <- corr_taxa$Species[corr_taxa$taxa == "Bryozoa"]
  data_bryozoa_roda <- rodarms_presence[, colnames(rodarms_presence) %in% bryozoa_species]
  
  species_occurence <- colSums(data_bryozoa_roda)/nrow(data_bryozoa_roda)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_bry_roda <- quantile(species_occurence$f, 0.25)
  
  #Bryozoa - runa
  data_bryozoa_runa <- runarms_presence[, colnames(runarms_presence) %in% bryozoa_species]
  
  species_occurence <- colSums(data_bryozoa_runa)/nrow(data_bryozoa_runa)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_bry_runa <- quantile(species_occurence$f, 0.25)
  
  #all_species - roda
  species_occurence <- colSums(rodarms_presence)/nrow(rodarms_presence)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_all_species_roda <- quantile(species_occurence$f, 0.25)
  
  #all_species - runa
  species_occurence <- colSums(runarms_presence)/nrow(runarms_presence)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_all_species_runa <- quantile(species_occurence$f, 0.25)
  
  #### Step 4bis: Compute Rarity threshold taxa par taxa ####  
  
  #Ascidians - roda
  data_mean_ascidiacea_roda <- rodarms_presence_mean[, colnames(rodarms_presence_mean) %in% ascidiacea_species]
  
  species_occurence_mean <- colSums(data_mean_ascidiacea_roda)/nrow(data_mean_ascidiacea_roda)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_asc_roda_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #Ascidian - runa
  data_mean_ascidiacea_runa <- runarms_presence_mean[, colnames(runarms_presence_mean) %in% ascidiacea_species]
  
  species_occurence_mean <- colSums(data_mean_ascidiacea_runa)/nrow(data_mean_ascidiacea_runa)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_asc_runa_mean <- quantile(species_occurence_mean$f, 0.25)
  
  
  #Porifera - roda
  data_mean_porifera_roda <- rodarms_presence_mean[, colnames(rodarms_presence_mean) %in% porifera_species]
  
  species_occurence_mean <- colSums(data_mean_porifera_roda)/nrow(data_mean_porifera_roda)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_por_roda_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #Porifera - runa
  data_mean_porifera_runa <- runarms_presence_mean[, colnames(runarms_presence_mean) %in% porifera_species]
  
  species_occurence_mean <- colSums(data_mean_porifera_runa)/nrow(data_mean_porifera_runa)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_por_runa_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #Bryozoa - roda
  data_mean_bryozoa_roda <- rodarms_presence_mean[, colnames(rodarms_presence_mean) %in% bryozoa_species]
  
  species_occurence_mean <- colSums(data_mean_bryozoa_roda)/nrow(data_mean_bryozoa_roda)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_bry_roda_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #Bryozoa - runa
  data_mean_bryozoa_runa <- runarms_presence_mean[, colnames(runarms_presence_mean) %in% bryozoa_species]
  
  species_occurence_mean <- colSums(data_mean_bryozoa_runa)/nrow(data_mean_bryozoa_runa)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_bry_runa_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #all_species - roda
  species_occurence_mean <- colSums(rodarms_presence_mean)/nrow(rodarms_presence_mean)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_all_species_roda_mean <- quantile(species_occurence_mean$f, 0.25)
  
  #all_species - runa
  species_occurence_mean <- colSums(runarms_presence_mean)/nrow(runarms_presence_mean)
  species_occurence_mean <- species_occurence_mean[species_occurence_mean != 0]
  species_occurence_mean <- data.frame(f = species_occurence_mean,
                                       Species = names(species_occurence_mean))
  threshold_all_species_runa_mean <- quantile(species_occurence_mean$f, 0.25)
  
  # species_list <- ascidiacea_species
  # group_name <- "Ascidiacea"
  # threshold <- threshold_asc_roda
  # threshold_mean <- threshold_asc_roda_mean
  
  #### Step 5: Fonctions  pour calculer la rareté pour chaque groupe ####
  compute_rarity_roda <- function(species_list, group_name, threshold, threshold_mean) {
    df <- rodarms_presence[, colnames(rodarms_presence) %in% species_list]
    freq <- colSums(df) / nrow(df)
    
    df_mean <- rodarms_presence_mean[, colnames(rodarms_presence_mean) %in% species_list]
    freq_mean <- colSums(df_mean) / nrow(df_mean)
    
    
    df$triplicat <- rodarms_meta$triplicat
    
    result <- sapply(species_list, function(espece) {
      subset <- df[df[[espece]] == 1, ]
      n_triplicats <- length(unique(subset$triplicat))
      if (n_triplicats <= 1) {
        return("Un seul triplicat")
      } else {
        return("Plusieurs triplicats")
      }
    })
    
    out <- data.frame(Species = names(freq), f = freq, f_mean = freq_mean, n_trip = result)
    out$taxa <- group_name
    out$threshold <- threshold
    out$threshold_mean <- threshold_mean
    
    out <- out %>%
      mutate(status = case_when(
        f == 0 ~ "absent",
        f <= threshold & f_mean <= threshold_mean ~ "rare",
        TRUE ~ "common"
      ))
    
    return(out)
  }
  
  # species_list <- ascidiacea_species
  # group_name <- "Ascidiacea"
  # threshold <- threshold_asc_runa
  # threshold_mean <- threshold_asc_runa_mean
  
  compute_rarity_runa <- function(species_list, group_name, threshold, threshold_mean) {
    df <- runarms_presence[, colnames(runarms_presence) %in% species_list]
    freq <- colSums(df) / nrow(df)
    
    df_mean <- runarms_presence_mean[, colnames(runarms_presence_mean) %in% species_list]
    freq_mean <- colSums(df_mean) / nrow(df_mean)
    
    df$triplicat <- runarms_meta$triplicat
    
    result <- sapply(species_list, function(espece) {
      subset <- df[df[[espece]] == 1, ]
      n_triplicats <- length(unique(subset$triplicat))
      if (n_triplicats <= 1) {
        return("Un seul triplicat")
      } else {
        return("Plusieurs triplicats")
      }
    })
    
    out <- data.frame(Species = names(freq), f = freq, f_mean = freq_mean, n_trip = result)
    out$taxa <- group_name
    out$threshold <- threshold
    out$threshold_mean <- threshold_mean
    
    out <- out %>%
      mutate(status = case_when(
        f == 0 ~ "absent",
        f <= threshold & f_mean <= threshold_mean ~ "rare",
        TRUE ~ "common"
      ))
    
    return(out)
  }
  
  #### Step 6: Attribuer les statuts de rareté aux espèces ####
  
  # Par groupe avec seuils spécifiques
  out_asc_runa <- compute_rarity_runa(ascidiacea_species, "Ascidiacea", threshold_asc_runa, threshold_asc_runa_mean)
  out_asc_roda <- compute_rarity_roda(ascidiacea_species, "Ascidiacea", threshold_asc_roda, threshold_asc_roda_mean)
  
  
  out_por_runa <- compute_rarity_runa(porifera_species, "Porifera", threshold_por_runa, threshold_por_runa_mean)
  out_por_roda <- compute_rarity_roda(porifera_species, "Porifera", threshold_por_roda, threshold_por_roda_mean)
  
  out_bry_runa <- compute_rarity_runa(bryozoa_species, "Bryozoa", threshold_bry_runa, threshold_bry_runa_mean)
  out_bry_roda <- compute_rarity_roda(bryozoa_species, "Bryozoa", threshold_bry_roda, threshold_bry_roda_mean)
  
  autres_groupes <- corr_taxa %>%
    filter(!taxa %in% c("Ascidiacea", "Porifera", "Bryozoa")) %>%
    pull(Species)
  
  out_others_runa <- compute_rarity_runa(autres_groupes, "Other", threshold_all_species_runa, threshold_all_species_runa_mean)
  out_others_roda <- compute_rarity_roda(autres_groupes, "Other", threshold_all_species_roda, threshold_all_species_roda_mean)
  
  recap_rarity_roda <- bind_rows(out_asc_roda, out_por_roda, out_bry_roda, out_others_roda)
  recap_rarity_runa <- bind_rows(out_asc_runa, out_por_runa, out_bry_runa, out_others_runa)
  
  recap_rarity_roda <- left_join(recap_rarity_roda, corr_taxa, by = "Species") %>%
    mutate(taxa = ifelse(taxa.x == "Other", taxa.y, taxa.x)) %>%
    select(Species, taxa, f, threshold, status)
  
  recap_rarity_runa <- left_join(recap_rarity_runa, corr_taxa, by = "Species") %>%
    mutate(taxa = ifelse(taxa.x == "Other", taxa.y, taxa.x)) %>%
    select(Species, taxa, f, threshold, status)
  
  species_rarity <- data.frame(Species = recap_rarity_roda$Species,
                               taxa = recap_rarity_roda$taxa,
                               Status_runa = recap_rarity_runa$status, 
                               Status_roda = recap_rarity_roda$status)
  
  
  # First, define your categories
  species_rarity_clean <- species_rarity %>%
    mutate(
      Reunion = case_when(
        Status_runa == "rare" ~ "Rare",
        Status_runa == "common" ~ "Not rare",
        Status_runa == "absent" ~ "Absent"
      ),
      Rodrigues = case_when(
        Status_roda == "rare" ~ "Rare",
        Status_roda == "common" ~ "Not rare",
        Status_roda == "absent" ~ "Absent"
      )
    ) %>%
    mutate(
      Rodrigues = ifelse(Reunion == "Absent" & Rodrigues %in% c("Rare", "Not rare"), "Exclusive", Rodrigues),
      Reunion = ifelse(Rodrigues == "Absent" & Reunion %in% c("Rare", "Not rare"), "Exclusive", Reunion)
    )
  
  species_rarity_clean <- species_rarity_clean[, !names(species_rarity_clean) %in% c("Status_runa", "Status_roda")]
  
  species_rarity_clean <- species_rarity_clean %>%
    mutate(
      taxon_group = case_when(
        taxa == "Ascidiacea" ~ "Ascidiacea",
        taxa == "Porifera" ~ "Porifera",
        taxa == "Bryozoa" ~ "Bryozoa",
        taxa %in% c("Foraminifera", "Cnidaria") ~ "Other",
        TRUE ~ "Other"  # au cas où d'autres apparaissent
      )
    )
  
  #### Step 7: décompter le nombre d'espèce rare ou pas dans les habitats ####
  
  
  # Calcul du nombre d'espèces par combinaison (shallow, mesophotic, taxon_group)
  summary_counts <- species_rarity_clean %>%
    group_by(Reunion, Rodrigues, taxon_group) %>%
    summarise(count = n(), .groups = "drop")
  
  # Calcul du total par groupe taxonomique (pour calculer les % par colonne)
  totals_by_group <- summary_counts %>%
    group_by(taxon_group) %>%
    summarise(total = sum(count), .groups = "drop")
  
  # Jointure et calcul des pourcentages
  summary_table <- summary_counts %>%
    left_join(totals_by_group, by = "taxon_group") %>%
    mutate(formatted = paste0(count, " (", round(100 * count / total), "%)")) %>%
    select(-count, -total) %>%
    pivot_wider(
      names_from = taxon_group,
      values_from = formatted,
      values_fill = "-"
    )
  
  # Ajout du total du nombre d'espèces pour chaque combinaison
  total_species <- species_rarity_clean %>%
    group_by(Reunion, Rodrigues) %>%
    summarise(Total_number_of_species = n(), .groups = "drop")
  
  # Fusion finale
  summary_table_final <- summary_table %>%
    left_join(total_species, by = c("Reunion", "Rodrigues")) %>%
    select(Reunion, Rodrigues, Ascidiacea, Porifera, Bryozoa, Other, Total_number_of_species)
  
  # Affichage
  print(summary_table_final)
  
  #### Step 8: Enregistrer le tableau ####
  
  # Reprend ton tableau résumé (déjà bien formaté avec les colonnes et valeurs type "n (x%)")
  # Ici on part du dataframe `summary_table`
  
  # Créer le flextable
  ft <- flextable(summary_table_final)
  ft <- autofit(ft)
  ft <- align(ft, align = "center", part = "all")
  ft <- set_table_properties(ft, layout = "autofit")
  
  # Créer le document Word et y insérer le tableau
  doc <- read_docx() %>%
    body_add_par("Species Rarity Summary Table", style = "heading 1") %>%
    body_add_flextable(ft)
  
  # Exporter le fichier Word
  print(doc, target = "outputs/species_rarity_summary_RODA_RUNA_04_06_2025.docx")
  return(NULL)
  }  