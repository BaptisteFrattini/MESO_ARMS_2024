
#' Rarity and commonness of morphospecies
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_rarity <- function(data_and_meta_clean){
  
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Rarity)
  
  # Computing ARMS plate scale rarity index ####
  
  ## Import and format data and metadata ####
  
  data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)
  
  
  # triplicats_exclure <- c("RUNARMS2", "RUNARMS3", "RUNARMS4", "RUNARMS6", "RUNARMS7", "RUNARMS8")
  # meta <- meta[!(meta$triplicat %in% triplicats_exclure), ]
  # data <- data[rownames(data) %in% meta$name, ]
  
  data <- data[, colSums(data) != 0]
  
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
  data_pa <- vegan::decostand(data, "pa")
  
  
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
  

  
  # campagne_name = "RUNARMS"
  
  ## Fonction pour calculer recap_rarity par campagne pour chaque taxon ####
  
  compute_recap_rarity_per_campagne <- function(campagne_name) {
    

    
    # 0. Prepare data
    data_camp <- subset(data, meta$campain == campagne_name)
    data_camp <- data_camp[,!colSums(data_camp) == 0]
    data_pa_camp <- vegan::decostand(data_camp, "pa")
    meta_camp <- subset(meta, meta$campain == campagne_name)
    
    
    
    # Ascidiacea
    asc_species <- corr_taxa$Species[corr_taxa$taxa == "Ascidiacea"]
    data_pa_camp_asc <- data_pa_camp[, colnames(data_pa_camp) %in% asc_species]
    occur_asc <- colSums(data_pa_camp_asc)
    
    # 1. Calcul du poid de rareté des espèces
    rarity.weights <- rWeights(occur_asc)
    
    # 2. Transforme en vecteur
    W <- rarity.weights$W
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rownames(rarity.weights)
    
    # 4. Les données sont a fournir a la fonction avec les espèces en lignes
    data_pa_camp_asc <- t(data_pa_camp_asc)
    
    # 5. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_pa_camp_asc, rarity.weights))
    
    # 6. On fait la somme des indices par site
    res_asc <- tapply(Community_rarity_index$IsrValue, meta_camp$triplicat, sum)
    
    # Porifera
    por_species <- corr_taxa$Species[corr_taxa$taxa == "Porifera"]
    data_pa_camp_por <- data_pa_camp[, colnames(data_pa_camp) %in% por_species]
    occur_por <- colSums(data_pa_camp_por)
    
    # 1. Calcul du poid de rareté des espèces
    rarity.weights <- rWeights(occur_por)
    
    # 2. Transforme en vecteur
    W <- rarity.weights$W
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rownames(rarity.weights)
    
    # 4. Les données sont a fournir a la fonction avec les espèces en lignes
    data_pa_camp_por <- t(data_pa_camp_por)
    
    # 5. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_pa_camp_por, rarity.weights))
    
    # 6. On fait la somme des indices par site
    res_por <- tapply(Community_rarity_index$IsrValue, meta_camp$triplicat, sum)
    
    # Bryozoa
    bry_species <- corr_taxa$Species[corr_taxa$taxa == "Bryozoa"]
    data_pa_camp_bry <- data_pa_camp[, colnames(data_pa_camp) %in% bry_species]
    occur_bry <- colSums(data_pa_camp_bry)
    
    # 1. Calcul du poid de rareté des espèces
    rarity.weights <- rWeights(occur_bry)
    
    # 2. Transforme en vecteur
    W <- rarity.weights$W
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rownames(rarity.weights)
    
    # 4. Les données sont a fournir a la fonction avec les espèces en lignes
    data_pa_camp_bry <- t(data_pa_camp_bry)
    
    # 5. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_pa_camp_bry, rarity.weights))
    
    # 6. On fait la somme des indices par site
    res_bry <- tapply(Community_rarity_index$IsrValue, meta_camp$triplicat, sum)
    
    
    recap <- data.frame(bind_rows(res_asc, res_por, res_bry))
    rownames(recap) <- c("asc", "por", "bry")
    return(recap)
    
  }

  runarms_index_plate <- compute_recap_rarity_per_campagne("RUNARMS")
  p50arms_index_plate <- compute_recap_rarity_per_campagne("P50ARMS")
  rodarms_index_plate <- compute_recap_rarity_per_campagne("RODARMS")
  
  # Computing site scale rarity index ####
  
  ## Import and format data and metadata ####
  
  data <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  
  # triplicats_exclure <- c("RUNARMS2", "RUNARMS3", "RUNARMS4", "RUNARMS6", "RUNARMS7", "RUNARMS8")
  # meta <- meta[!(meta$triplicat %in% triplicats_exclure), ]
  # data <- data[rownames(data) %in% meta$name, ]
  
  data <- data[, colSums(data) != 0]
  
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
  data_pa <- vegan::decostand(data, "pa")
  
  
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
  
  
  
  # campagne_name = "P50ARMS"
  
  ## Fonction pour calculer recap_rarity par campagne pour chaque taxon ####
  
  compute_recap_rarity_per_campagne <- function(campagne_name) {
    
    
    
    # 0. Prepare data
    data_camp <- subset(data, meta$campain == campagne_name)
    data_camp <- data_camp[,!colSums(data_camp) == 0]
    data_pa_camp <- vegan::decostand(data_camp, "pa")
    meta_camp <- subset(meta, meta$campain == campagne_name)
    
    
    
    # Ascidiacea
    asc_species <- corr_taxa$Species[corr_taxa$taxa == "Ascidiacea"]
    data_pa_camp_asc <- data_pa_camp[, colnames(data_pa_camp) %in% asc_species]
    occur_asc <- colSums(data_pa_camp_asc)
    
    # 1. Calcul du poid de rareté des espèces
    rarity.weights <- rWeights(occur_asc)
    
    # 2. Transforme en vecteur
    W <- rarity.weights$W
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rownames(rarity.weights)
    
    # 4. Les données sont a fournir a la fonction avec les espèces en lignes
    data_pa_camp_asc <- t(data_pa_camp_asc)
    
    # 5. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_pa_camp_asc, rarity.weights))
    
    # 6. On fait la somme des indices par site
    res_asc <- tapply(Community_rarity_index$IsrValue, meta_camp$triplicat, sum)
    
    # Porifera
    por_species <- corr_taxa$Species[corr_taxa$taxa == "Porifera"]
    data_pa_camp_por <- data_pa_camp[, colnames(data_pa_camp) %in% por_species]
    occur_por <- colSums(data_pa_camp_por)
    
    # 1. Calcul du poid de rareté des espèces
    rarity.weights <- rWeights(occur_por)
    
    # 2. Transforme en vecteur
    W <- rarity.weights$W
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rownames(rarity.weights)
    
    # 4. Les données sont a fournir a la fonction avec les espèces en lignes
    data_pa_camp_por <- t(data_pa_camp_por)
    
    # 5. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_pa_camp_por, rarity.weights))
    
    # 6. On fait la somme des indices par site
    res_por <- tapply(Community_rarity_index$IsrValue, meta_camp$triplicat, sum)
    
    # Bryozoa
    bry_species <- corr_taxa$Species[corr_taxa$taxa == "Bryozoa"]
    data_pa_camp_bry <- data_pa_camp[, colnames(data_pa_camp) %in% bry_species]
    occur_bry <- colSums(data_pa_camp_bry)
    
    # 1. Calcul du poid de rareté des espèces
    rarity.weights <- rWeights(occur_bry)
    
    # 2. Transforme en vecteur
    W <- rarity.weights$W
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rownames(rarity.weights)
    
    # 4. Les données sont a fournir a la fonction avec les espèces en lignes
    data_pa_camp_bry <- t(data_pa_camp_bry)
    
    # 5. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_pa_camp_bry, rarity.weights))
    
    # 6. On fait la somme des indices par site
    res_bry <- tapply(Community_rarity_index$IsrValue, meta_camp$triplicat, sum)
    
    
    recap <- data.frame(bind_rows(res_asc, res_por, res_bry))
    rownames(recap) <- c("asc", "por", "bry")
    return(recap)
    
  }
  
  runarms_index_ARMS <- compute_recap_rarity_per_campagne("RUNARMS")
  p50arms_index_ARMS <- compute_recap_rarity_per_campagne("P50ARMS")
  rodarms_index_ARMS <- compute_recap_rarity_per_campagne("RODARMS")
  
  
  
  

  # 8. On fait de meme que 5. , 6. et 7. mais sur la matrice d'abondance
  data_ab_camp <- as.matrix(t(data_camp))
  
  Community_rarity_index_ab <- as.data.frame(Irr(data_ab_camp, W = W, abundance = T, Wmin = min(W), Wmax = max(W)))
  
  res_ab <- tapply(Community_rarity_index_ab$Irr, meta_camp$triplicat, function(x) sum(x, na.rm = TRUE))
  
  barplot(res_ab)
  return(NULL)  
}



