#' Rarity and commonness of morphospecies
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_rarity <- function(data_and_meta_clean_fullsites, gps_sites, runa_map, roda_map, runa_reef, roda_reef){
  
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites")
  # gps_sites = targets::tar_read("data_gps_sites")
  # runa_map = targets::tar_read("map_runa")
  # roda_map = targets::tar_read("map_roda")
  # runa_reef = targets::tar_read("reef_runa")
  # roda_reef = targets::tar_read("reef_roda")

  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Rarity)
  library(betapart)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(ggsignif)
  library(ggpubr)
  library(adespatial)
  library(sf)
  library(ggrepel)
  library(patchwork)
  library(ggspatial)
  library(vegan)
  
  
  # Computing ARMS plate scale rarity index ####
  
  ## Import and format data and metadata ####
  
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)
  
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
  
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[, msp_list_filter]
  data_mean_pa <- vegan::decostand(data_mean, "pa")
  
  # campagne_name = "RUNARMS"
  
  ## Fonction pour calculer recap_rarity par campagne pour chaque taxon ####
  
  compute_recap_rarity_per_campagne <- function(campagne_name) {
    

    
    # 0. Prepare data
    data_camp <- subset(data, meta$campain == campagne_name)
    data_camp <- data_camp[,!colSums(data_camp) == 0]
    data_pa_camp <- vegan::decostand(data_camp, "pa")
    meta_camp <- subset(meta, meta$campain == campagne_name)
    
    # 0bis. Prepare data_mean
    
    data_mean_camp <- subset(data_mean, meta_mean$campain == campagne_name)
    data_mean_camp <- data_mean_camp[,!colSums(data_mean_camp) == 0]
    data_mean_pa_camp <- vegan::decostand(data_mean_camp, "pa")
    meta_mean_camp <- subset(meta_mean, meta_mean$campain == campagne_name)
    
    
    # All species
    occur_all_ARMS <- colSums(data_mean_pa_camp)
    
    occur_all_plates <- colSums(data_pa_camp)
    
    # 1. Les données sont a fournir a la fonction avec les espèces en lignes
    data_mean_pa_camp_all <- t(data_mean_pa_camp)
    data_pa_camp_all <- t(data_pa_camp)
    data_pa_camp_all <-  data_pa_camp_all[,!colSums(data_pa_camp_all) == 0]
    
    # 2. Calcul du poid de rareté des espèces
    rarity.weights.ARMS <- rWeights(occur_all_ARMS, rCutoff = "Leroy", assemblages = data_mean_pa_camp_all)
    rarity.weights.ARMS$Species <- rownames(rarity.weights.ARMS)
    
    rarity.weights.plates <- rWeights(occur_all_plates,  rCutoff = "Leroy", assemblages = data_pa_camp_all)
    rarity.weights.plates$Species <- rownames(rarity.weights.plates)
    
    rarity.weights <- left_join(rarity.weights.ARMS, rarity.weights.plates, by = "Species", suffix = c(".ARMS", ".plates"))
      
    rarity.weights[rarity.weights$R.ARMS != rarity.weights$R.plates, "Species"]
    
    # 3. Transforme en vecteur
    W <- rarity.weights$W.ARMS + rarity.weights$W.plates
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rarity.weights$Species
    
    # 4. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_mean_pa_camp_all, W))
    
    # Community_rarity_index <- as.data.frame(Irr(data_mean_pa_camp_all, W))
    
    # 5. On fait la somme des indices par site
    res_all <- tapply(Community_rarity_index$IsrValue, meta_mean_camp$triplicat, mean)
    
    # res_all <- tapply(Community_rarity_index$Irr, meta_mean_camp$triplicat, mean)
    
    res_all <- data.frame(site = names(res_all),
                          R = res_all)
    
    
    # 6. On calcule le LCBD
    data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
    meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
    matrix.bray <- vegan::vegdist(data_mean, method = "bray")
    spe.beta <- adespatial::LCBD.comp(matrix.bray, sqrt.D = TRUE)
    names(spe.beta$LCBD) <- rownames(data_mean)
    spe.beta$Site <- meta_mean$triplicat  
    spe.beta <- data.frame(site = spe.beta$Site, LCBD = spe.beta$LCBD) 
    beta_agg <- spe.beta %>%
      group_by(site) %>%
      summarise(LCBD = mean(LCBD))
    
    res_all <- left_join(res_all, beta_agg, by = "site")
    colnames(res_all) <- c("site", "R_all", "LCBD_all")
    
    # Ascidiacea
    asc_species <- corr_taxa$Species[corr_taxa$taxa == "Ascidiacea"]
    
    data_pa_camp_asc <- data_pa_camp[, colnames(data_pa_camp) %in% asc_species]
    data_mean_pa_camp_asc <- data_mean_pa_camp[, colnames(data_mean_pa_camp) %in% asc_species]
    
    occur_asc_ARMS <- colSums(data_mean_pa_camp_asc)
    
    occur_asc_plates <- colSums(data_pa_camp_asc)
    
    # 1. Les données sont a fournir a la fonction avec les espèces en lignes
    data_mean_pa_camp_asc <- t(data_mean_pa_camp_asc)
    data_pa_camp_asc <- t(data_pa_camp_asc)
    data_pa_camp_asc <-  data_pa_camp_asc[,!colSums(data_pa_camp_asc) == 0]
    
    # 2. Calcul du poid de rareté des espèces
    rarity.weights.ARMS <- rWeights(occur_asc_ARMS, rCutoff = "Leroy", assemblages = data_mean_pa_camp_asc)
    rarity.weights.ARMS$Species <- rownames(rarity.weights.ARMS)
    
    rarity.weights.plates <- rWeights(occur_asc_plates,  rCutoff = "Leroy", assemblages = data_pa_camp_asc)
    rarity.weights.plates$Species <- rownames(rarity.weights.plates)
    
    rarity.weights <- left_join(rarity.weights.ARMS, rarity.weights.plates, by = "Species", suffix = c(".ARMS", ".plates"))
    
    rarity.weights[rarity.weights$R.ARMS != rarity.weights$R.plates, "Species"]
    
    # 3. Transforme en vecteur
    W <- rarity.weights$W.ARMS + rarity.weights$W.plates
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rarity.weights$Species
    
    # 4. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_mean_pa_camp_asc, W))
    
    # Community_rarity_index <- as.data.frame(Irr(data_mean_pa_camp_asc, W))
    
    # 5. On fait la somme des indices par site
    res_asc <- tapply(Community_rarity_index$IsrValue, meta_mean_camp$triplicat, mean)
    
    # res_asc <- tapply(Community_rarity_index$Irr, meta_mean_camp$triplicat, mean)
    
    res_asc <- data.frame(site = names(res_asc),
                          R = res_asc)
    
    # 7. On calcule le LCBD
    data_mean_asc <- data_mean[, colnames(data_mean) %in% asc_species]
    matrix.bray <- vegan::vegdist(data_mean_asc, method = "bray")
    spe.beta <- adespatial::LCBD.comp(matrix.bray, sqrt.D = TRUE)
    names(spe.beta$LCBD) <- rownames(data_mean)
    spe.beta$Site <- meta_mean$triplicat  
    spe.beta <- data.frame(site = spe.beta$Site, LCBD = spe.beta$LCBD) 
    beta_agg <- spe.beta %>%
      group_by(site) %>%
      summarise(LCBD = mean(LCBD))
    
    res_asc <- left_join(res_asc, beta_agg, by = "site")
    colnames(res_asc) <- c("site", "R_asc", "LCBD_asc")
    
    # Porifera
    por_species <- corr_taxa$Species[corr_taxa$taxa == "Porifera"]
    
    data_pa_camp_por <- data_pa_camp[, colnames(data_pa_camp) %in% por_species]
    data_mean_pa_camp_por <- data_mean_pa_camp[, colnames(data_mean_pa_camp) %in% por_species]
    
    occur_por_ARMS <- colSums(data_mean_pa_camp_por)
    
    occur_por_plates <- colSums(data_pa_camp_por)
    
    # 1. Les données sont a fournir a la fonction avec les espèces en lignes
    data_mean_pa_camp_por <- t(data_mean_pa_camp_por)
    data_pa_camp_por <- t(data_pa_camp_por)
    data_pa_camp_por <-  data_pa_camp_por[,!colSums(data_pa_camp_por) == 0]
    
    # 2. Calcul du poid de rareté des espèces
    rarity.weights.ARMS <- rWeights(occur_por_ARMS, rCutoff = "Leroy", assemblages = data_mean_pa_camp_por)
    rarity.weights.ARMS$Species <- rownames(rarity.weights.ARMS)
    
    rarity.weights.plates <- rWeights(occur_por_plates,  rCutoff = "Leroy", assemblages = data_pa_camp_por)
    rarity.weights.plates$Species <- rownames(rarity.weights.plates)
    
    rarity.weights <- left_join(rarity.weights.ARMS, rarity.weights.plates, by = "Species", suffix = c(".ARMS", ".plates"))
    
    rarity.weights[rarity.weights$R.ARMS != rarity.weights$R.plates, "Species"]
    
    # 3. Transforme en vecteur
    W <- rarity.weights$W.ARMS + rarity.weights$W.plates
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rarity.weights$Species
    
    # 4. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_mean_pa_camp_por, W))
    
    # Community_rarity_index <- as.data.frame(Irr(data_mean_pa_camp_por, W))
    
    # 5. On fait la somme des indices par site
    res_por <- tapply(Community_rarity_index$IsrValue, meta_mean_camp$triplicat, mean)
    
    # res_por <- tapply(Community_rarity_index$Irr, meta_mean_camp$triplicat, mean)
    
    res_por <- data.frame(site = names(res_por),
                          R = res_por)
    
    # 7. On calcule le LCBD
    data_mean_por <- data_mean[, colnames(data_mean) %in% por_species]
    matrix.bray <- vegan::vegdist(data_mean_por, method = "bray")
    spe.beta <- adespatial::LCBD.comp(matrix.bray, sqrt.D = TRUE)
    names(spe.beta$LCBD) <- rownames(data_mean)
    spe.beta$Site <- meta_mean$triplicat  
    spe.beta <- data.frame(site = spe.beta$Site, LCBD = spe.beta$LCBD) 
    beta_agg <- spe.beta %>%
      group_by(site) %>%
      summarise(LCBD = mean(LCBD))
    
    res_por <- left_join(res_por, beta_agg, by = "site")
    colnames(res_por) <- c("site", "R_por", "LCBD_por")
    
    
    # Bryozoa
    bry_species <- corr_taxa$Species[corr_taxa$taxa == "Bryozoa"]
    
    data_pa_camp_bry <- data_pa_camp[, colnames(data_pa_camp) %in% bry_species]
    data_mean_pa_camp_bry <- data_mean_pa_camp[, colnames(data_mean_pa_camp) %in% bry_species]
    
    occur_bry_ARMS <- colSums(data_mean_pa_camp_bry)
    
    occur_bry_plates <- colSums(data_pa_camp_bry)
    
    # 1. Les données sont a fournir a la fonction avec les espèces en lignes
    data_mean_pa_camp_bry <- t(data_mean_pa_camp_bry)
    data_pa_camp_bry <- t(data_pa_camp_bry)
    data_pa_camp_bry <-  data_pa_camp_bry[,!colSums(data_pa_camp_bry) == 0]
    
    # 2. Calcul du poid de rareté des espèces
    rarity.weights.ARMS <- rWeights(occur_bry_ARMS, rCutoff = "Leroy", assemblages = data_mean_pa_camp_bry)
    rarity.weights.ARMS$Species <- rownames(rarity.weights.ARMS)
    
    rarity.weights.plates <- rWeights(occur_bry_plates,  rCutoff = "Leroy", assemblages = data_pa_camp_bry)
    rarity.weights.plates$Species <- rownames(rarity.weights.plates)
    
    rarity.weights <- left_join(rarity.weights.ARMS, rarity.weights.plates, by = "Species", suffix = c(".ARMS", ".plates"))
    
    rarity.weights[rarity.weights$R.ARMS != rarity.weights$R.plates, "Species"]
    
    # 3. Transforme en vecteur
    W <- rarity.weights$W.ARMS + rarity.weights$W.plates
    
    # 3. Auquel on re-associe un nom d'espèce
    names(W) <- rarity.weights$Species
    
    # 4. On calcule sur la base de la matrice presence absence d'origine, l'indice de rareté de chaque face de plaque
    Community_rarity_index <- as.data.frame(Isr(data_mean_pa_camp_bry, W))
    
    # Community_rarity_index <- as.data.frame(Irr(data_mean_pa_camp_bry, W))
    
    # 5. On fait la somme des indices par site
    res_bry <- tapply(Community_rarity_index$IsrValue, meta_mean_camp$triplicat, mean)
    
    # res_bry <- tapply(Community_rarity_index$Irr, meta_mean_camp$triplicat, mean)
    
    res_bry <- data.frame(site = names(res_bry),
                          R = res_bry)
    
    # 7. On calcule le LCBD
    data_mean_bry <- data_mean[, colnames(data_mean) %in% bry_species]
    matrix.bray <- vegan::vegdist(data_mean_bry, method = "bray")
    spe.beta <- adespatial::LCBD.comp(matrix.bray, sqrt.D = TRUE)
    names(spe.beta$LCBD) <- rownames(data_mean)
    spe.beta$Site <- meta_mean$triplicat  
    spe.beta <- data.frame(site = spe.beta$Site, LCBD = spe.beta$LCBD) 
    beta_agg <- spe.beta %>%
      group_by(site) %>%
      summarise(LCBD = mean(LCBD))
    
    res_bry <- left_join(res_bry, beta_agg, by = "site")
    colnames(res_bry) <- c("site", "R_bry", "LCBD_bry")
    
    
    recap <- res_all %>%
      left_join(res_asc, by = "site") %>%
      left_join(res_por, by = "site") %>%
      left_join(res_bry, by = "site")

    # colnames(recap) <- c("all","asc", "por", "bry")
    return(recap)
    
  }

  runarms_index_plate <- compute_recap_rarity_per_campagne("RUNARMS")
  p50arms_index_plate <- compute_recap_rarity_per_campagne("P50ARMS")
  rodarms_index_plate <- compute_recap_rarity_per_campagne("RODARMS")
  
  tab <- as.data.frame(rbind(runarms_index_plate, p50arms_index_plate, rodarms_index_plate))
  # tab$site <- rownames(tab)
  
  #### Import the GPS data ####
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE)
  
  #### Import the map shapefiles ####
  runa_map <- st_read(runa_map)
  roda_map <- st_read(roda_map)
  
  
  
  data_gps <- data_gps %>%
    mutate(
      label = case_when(
        grepl("P50ARMS", Site) ~ gsub("P50ARMS", "P50A", Site),
        grepl("RUNARMS", Site) ~ gsub("RUNARMS", "RUNA", Site),
        grepl("RODARMS", Site) ~ gsub("RODARMS", "RODA", Site),
        TRUE ~ Site
      )
    )
  
  data_gps <- data_gps %>%
    mutate(
      site = case_when(
        grepl("RUNA", Site) ~ gsub("RUNA", "RUNARMS", Site),
        grepl("RODA", Site) ~ gsub("RODA", "RODARMS", Site),
        TRUE ~ Site
      )
    )
  
  map_data <- left_join(tab, data_gps, by = "site")
  
  map_data_runa <- map_data[grepl("RUNA", map_data$label),]
  
  # Convert to sf object
  map_data_runa_sf <- st_as_sf(map_data_runa, coords = c("Longitude", "Latitude"), crs = 4326)
  st_crs(map_data_runa_sf) <- 4326
  
  runa_map <- st_transform(runa_map, crs = 4326)
  
  map_data_runa_sf <- st_transform(map_data_runa_sf, st_crs(runa_map))
  st_crs(runa_map)
  st_crs(map_data_runa_sf)
  
  map_data_runa_sf <- map_data_runa_sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2])
  
  lcbd_all <- c(map_data$LCBD_all)
  min_LCBD <- min(lcbd_all)
  max_LCBD <- max(lcbd_all)
  
  scale_lcbd <- scale_size_continuous(name = "LCBD", range = c(2, 10), limits = c(min_LCBD, max_LCBD))
  
  #### Import the reef shapefiles ####
  
  runa_reef <- st_read(runa_reef)
  roda_reef <- st_read(roda_reef)
  
  runa_reef <- st_transform(runa_reef, crs = 4326)
  roda_reef <- st_transform(roda_reef, crs = 4326)
  
  #### RUNA map settings ####
  map_data_p50a <- map_data[grepl("P50A", map_data$Site),]
  map_data_p50a_sf <- st_as_sf(map_data_p50a, coords = c("Longitude", "Latitude"), crs = 4326)
  st_crs(map_data_p50a_sf) <- 4326
  map_data_p50a_sf <- st_transform(map_data_p50a_sf, st_crs(runa_map))
  st_crs(runa_map)
  st_crs(map_data_p50a_sf)
  map_data_p50a_sf <- map_data_p50a_sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2])
  
  #### RODA map settings ####
  map_data_roda <- map_data[grepl("RODA", map_data$Site),]
  map_data_roda_sf <- st_as_sf(map_data_roda, coords = c("Longitude", "Latitude"), crs = 4326)
  st_crs(roda_map)
  st_crs(map_data_roda_sf)
  map_data_roda_sf <- map_data_roda_sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2])
  
  
  # -- Créer l'objet sf principal --
  map_data_sf <- st_as_sf(map_data, coords = c("Longitude", "Latitude"), crs = 4326)
  map_data_sf <- st_transform(map_data_sf, st_crs(runa_map)) %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2])
  
  # -- Ajuster les coordonnées pour éviter les chevauchements --
  map_data_sf <- map_data_sf %>%
    mutate(
      Latitude = case_when(
        site %in% c("RUNARMS2") ~ Latitude + 0.006,
        site %in% c("RUNARMS3") ~ Latitude - 0.006,
        site %in% c("RUNARMS4") ~ Latitude + 0.005,
        site %in% c("RUNARMS5") ~ Latitude - 0.005,
        site %in% c("RUNARMS6") ~ Latitude - 0.005,
        TRUE ~ Latitude
      )
    )
  
  map_data_sf <- map_data_sf %>%
    mutate(
      Longitude = case_when(
        site %in% c("RUNARMS2") ~ Longitude - 0.005,
        site %in% c("RUNARMS3") ~ Longitude + 0.007,
        TRUE ~ Longitude
      )
    )
  
  # -- Scales communes --
  scale_r <- scale_size_continuous(
    name = "ISR", 
    range = c(1, 12),
    limits = c(min(map_data$R_all, na.rm = TRUE), max(map_data$R_all, na.rm = TRUE))
  )
  
 
  
  
  scale_lcbd <- scale_fill_gradientn(
    colours = rev(viridis::viridis(10)),
    name = "LCBD",
    limits = c(min(map_data$LCBD_all, na.rm = TRUE), max(map_data$LCBD_all, na.rm = TRUE)),
    guide = guide_colorbar(
      title.hjust = 0.5,
      barwidth = unit(4, "cm"),     # allonge la barre
      barheight = unit(0.4, "cm"),  # réduit l’épaisseur (utile pour une légende horizontale)
      label.theme = element_text(size = 8)  # taille des chiffres
    )
  )
  
  # Define common width and heights
  map_width <- 0.45
  map_height <- 0.45
  
  # Center coordinates for each map
  center_p50a <- c(x = mean(c(55.15, 55.6)), y = mean(c(-21.4, -21)))
  center_roda <- c(x = mean(c(63.20, 63.58)), y = mean(c(-19.90, -19.58)))
  # (Add RUNA similarly if needed)
  
  # Create coordinate limits
  p50a_xlim <- c(center_p50a['x'] - map_width / 2, center_p50a['x'] + map_width / 2)
  p50a_ylim <- c(center_p50a['y'] - map_height / 2, center_p50a['y'] + map_height / 2)
  
  roda_xlim <- c(center_roda['x'] - map_width / 2, center_roda['x'] + map_width / 2)
  roda_ylim <- c(center_roda['y'] - map_height / 2, center_roda['y'] + map_height / 2)
  
  # Then apply in each coord_sf:
  coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE)
  coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE)
  
  # --- 9. Création des cartes ---
  runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
  p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
  roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
  
  # -- Map the RUNA data --
  map_data_runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
  
  
  uu <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_runa_sf,
               aes(x = Longitude, y = Latitude,size = R_all, fill = LCBD_all),
               shape = 21, stroke = 0.6) +
    geom_text_repel(
      data = map_data_runa_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = "RUNA") +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Map the P50A data --
  map_data_p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
  
  vv <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_p50a_sf,
               aes(x = Longitude, y = Latitude, size = R_all, fill = LCBD_all),
               shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_p50a_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = "P50A") +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Map the RODA data --
  map_data_roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
  
  ww <- ggplot() +
    geom_sf(data = roda_map, fill = "grey85", color = "black") +
    geom_sf(data = roda_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_roda_sf,
               aes(x = Longitude, y = Latitude,size = R_all, fill = LCBD_all),
               shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_roda_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) + 
    # coord_sf(xlim = c(63.25, 63.55), ylim = c(-19.85, -19.62), expand = FALSE) +
    coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = "RODA") +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Combinaison finale --
  
  tt <- (uu | vv | ww) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Separate taxa maps ####
  
  ## Ascidiacea ####
  scale_r <- scale_size_continuous(
    name = "ISR", 
    range = c(1, 12),
    limits = c(min(map_data$R_asc, na.rm = TRUE), max(map_data$R_asc, na.rm = TRUE))
  )
  
  scale_lcbd <- scale_fill_gradientn(
    colours = rev(viridis::viridis(10)),
    name = "LCBD",
    limits = c(min(map_data$LCBD_asc, na.rm = TRUE), max(map_data$LCBD_asc, na.rm = TRUE)),
    guide = guide_colorbar(
      title.hjust = 0.5,
      barwidth = unit(4, "cm"),     # asconge la barre
      barheight = unit(0.4, "cm"),  # réduit l’épaisseur (utile pour une légende horizontale)
      label.theme = element_text(size = 8)  # taille des chiffres
    )
  )
  
  # Define common width and height
  map_width <- 0.45
  map_height <- 0.45
  
  # Center coordinates for each map
  center_p50a <- c(x = mean(c(55.15, 55.6)), y = mean(c(-21.4, -21)))
  center_roda <- c(x = mean(c(63.20, 63.58)), y = mean(c(-19.90, -19.58)))
  # (Add RUNA similarly if needed)
  
  # Create coordinate limits
  p50a_xlim <- c(center_p50a['x'] - map_width / 2, center_p50a['x'] + map_width / 2)
  p50a_ylim <- c(center_p50a['y'] - map_height / 2, center_p50a['y'] + map_height / 2)
  
  roda_xlim <- c(center_roda['x'] - map_width / 2, center_roda['x'] + map_width / 2)
  roda_ylim <- c(center_roda['y'] - map_height / 2, center_roda['y'] + map_height / 2)
  
  # Then apply in each coord_sf:
  coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE)
  coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE)
  
  # --- 9. Création des cartes ---
  runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
  p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
  roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
  
  # -- Map the RUNA data --
  map_data_runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
  
  
  uu <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_runa_sf,
               aes(x = Longitude, y = Latitude,size = R_asc, fill = LCBD_asc),
               shape = 21, stroke = 0.6) +
    geom_text_repel(
      data = map_data_runa_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("RUNA –", "Ascidiacea")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Map the P50A data --
  map_data_p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
  
  vv <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_p50a_sf,
               aes(x = Longitude, y = Latitude, size = R_asc, fill = LCBD_asc),
               shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_p50a_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("P50A –", "Ascidiacea")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Map the RODA data --
  map_data_roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
  
  ww <- ggplot() +
    geom_sf(data = roda_map, fill = "grey85", color = "black") +
    geom_sf(data = roda_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_roda_sf,
               aes(x = Longitude, y = Latitude,size = R_asc, fill = LCBD_asc),
               shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_roda_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) + 
    # coord_sf(xlim = c(63.25, 63.55), ylim = c(-19.85, -19.62), expand = FALSE) +
    coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("RODA –", "Ascidiacea")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Combinaison finale --
  
  
  xx <- (uu | vv | ww) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  
  ## Bryozoa ####
  scale_r <- scale_size_continuous(
    name = "ISR", 
    range = c(1, 12),
    limits = c(min(map_data$R_bry, na.rm = TRUE), max(map_data$R_bry, na.rm = TRUE))
  )
  
  
  scale_lcbd <- scale_fill_gradientn(
    colours = rev(viridis::viridis(10)),
    name = "LCBD",
    limits = c(min(map_data$LCBD_bry, na.rm = TRUE), max(map_data$LCBD_bry, na.rm = TRUE)),
    guide = guide_colorbar(
      title.hjust = 0.5,
      barwidth = unit(4, "cm"),     # bryonge la barre
      barheight = unit(0.4, "cm"),  # réduit l’épaisseur (utile pour une légende horizontale)
      label.theme = element_text(size = 8)  # taille des chiffres
    )
  )
  
  # Define common width and height
  map_width <- 0.45
  map_height <- 0.45
  
  # Center coordinates for each map
  center_p50a <- c(x = mean(c(55.15, 55.6)), y = mean(c(-21.4, -21)))
  center_roda <- c(x = mean(c(63.20, 63.58)), y = mean(c(-19.90, -19.58)))
  # (Add RUNA similarly if needed)
  
  # Create coordinate limits
  p50a_xlim <- c(center_p50a['x'] - map_width / 2, center_p50a['x'] + map_width / 2)
  p50a_ylim <- c(center_p50a['y'] - map_height / 2, center_p50a['y'] + map_height / 2)
  
  roda_xlim <- c(center_roda['x'] - map_width / 2, center_roda['x'] + map_width / 2)
  roda_ylim <- c(center_roda['y'] - map_height / 2, center_roda['y'] + map_height / 2)
  
  # Then apply in each coord_sf:
  coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE)
  coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE)
  
  # --- 9. Création des cartes ---
  runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
  p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
  roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
  
  # -- Map the RUNA data --
  map_data_runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
  
  
  uu <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_runa_sf,
               aes(x = Longitude, y = Latitude,size = R_bry, fill = LCBD_bry),
               shape = 21, stroke = 0.6) +
    geom_text_repel(
      data = map_data_runa_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("RUNA –", "Bryozoa")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Map the P50A data --
  map_data_p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
  
  vv <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_p50a_sf,
               aes(x = Longitude, y = Latitude, size = R_bry, fill = LCBD_bry),
               shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_p50a_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("P50A –", "Bryozoa")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Map the RODA data --
  map_data_roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
  
  ww <- ggplot() +
    geom_sf(data = roda_map, fill = "grey85", color = "black") +
    geom_sf(data = roda_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_roda_sf,
               aes(x = Longitude, y = Latitude,size = R_bry, fill = LCBD_bry),
               shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_roda_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) + 
    # coord_sf(xlim = c(63.25, 63.55), ylim = c(-19.85, -19.62), expand = FALSE) +
    coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("RODA –", "Bryozoa")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Combinaison finale --
  
  
  yy <- (uu | vv | ww) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ## Porifera ####
  

  scale_r <- scale_size_continuous(
    name = "ISR", 
    range = c(1, 12),
    limits = c(min(map_data$R_por, na.rm = TRUE), max(map_data$R_por, na.rm = TRUE))
  )
  
  
  scale_lcbd <- scale_fill_gradientn(
    colours = rev(viridis::viridis(10)),
    name = "LCBD",
    limits = c(min(map_data$LCBD_por, na.rm = TRUE), max(map_data$LCBD_por, na.rm = TRUE)),
    guide = guide_colorbar(
      title.hjust = 0.5,
      barwidth = unit(4, "cm"),     # poronge la barre
      barheight = unit(0.4, "cm"),  # réduit l’épaisseur (utile pour une légende horizontale)
      label.theme = element_text(size = 8)  # taille des chiffres
    )
  )
  
  # Define common width and height
  map_width <- 0.45
  map_height <- 0.45
  
  # Center coordinates for each map
  center_p50a <- c(x = mean(c(55.15, 55.6)), y = mean(c(-21.4, -21)))
  center_roda <- c(x = mean(c(63.20, 63.58)), y = mean(c(-19.90, -19.58)))
  # (Add RUNA similarly if needed)
  
  # Create coordinate limits
  p50a_xlim <- c(center_p50a['x'] - map_width / 2, center_p50a['x'] + map_width / 2)
  p50a_ylim <- c(center_p50a['y'] - map_height / 2, center_p50a['y'] + map_height / 2)
  
  roda_xlim <- c(center_roda['x'] - map_width / 2, center_roda['x'] + map_width / 2)
  roda_ylim <- c(center_roda['y'] - map_height / 2, center_roda['y'] + map_height / 2)
  
  # Then apply in each coord_sf:
  coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE)
  coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE)
  
  # --- 9. Création des cartes ---
  runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
  p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
  roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
  
  # -- Map the RUNA data --
  map_data_runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
  
  
  uu <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_runa_sf,
               aes(x = Longitude, y = Latitude,size = R_por, fill = LCBD_por),
               shape = 21, stroke = 0.6) +
    geom_text_repel(
      data = map_data_runa_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("RUNA –", "Porifera")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Map the P50A data --
  map_data_p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
  
  vv <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_p50a_sf,
               aes(x = Longitude, y = Latitude, size = R_por, fill = LCBD_por),
               shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_p50a_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("P50A –", "Porifera")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Map the RODA data --
  map_data_roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
  
  ww <- ggplot() +
    geom_sf(data = roda_map, fill = "grey85", color = "black") +
    geom_sf(data = roda_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_roda_sf,
               aes(x = Longitude, y = Latitude,size = R_por, fill = LCBD_por),
               shape = 21, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_roda_sf,
      aes(x = Longitude, y = Latitude, label = label),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) + 
    # coord_sf(xlim = c(63.25, 63.55), ylim = c(-19.85, -19.62), expand = FALSE) +
    coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = paste("RODA –", "Porifera")) +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Combinaison finale --
  
  zz <- (uu | vv | ww) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  final_plot <- (tt / xx / yy / zz) +
    plot_layout(guides = "collect")
  
  ggsave("outputs/Cartes - LCBD_Rarity/map_09_07_2025_IRR.pdf", plot = final_plot, width = 0.75*18, height = 18)
  
  # Computing correlation plto btw IRR and LCBD ####
  
  plot(map_data_sf$LCBD_all, map_data_sf$R_all)
  
  cor_test <- cor.test(map_data_sf$LCBD_all, map_data_sf$R_all, method = "pearson")
  
  # Extraire le coefficient et la p-value
  cor_value <- round(cor_test$estimate, 2)
  p_value <- cor_test$p.value
  method <- "Pearson"
  
  # Création du graphique
  dd <- ggplot(map_data_sf, aes(x = LCBD_all, y =R_all)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = min(map_data_sf$LCBD_all), y = max(map_data_sf$R_all),
             label = paste0(method, " r = ", cor_value, 
                            "\n(p = ", signif(p_value, 3), ")"),
             hjust = 0, vjust = 1, size = 5) +
    labs(
      x = "LCBD values",
      y = "ISR values",
      title = ""
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16),  # Titre axe X
      axis.title.y = element_text(size = 16),  # Titre axe Y
      axis.text.x  = element_text(size = 16),  # Valeurs axe X
      axis.text.y  = element_text(size = 16)
    )
    
  ggsave("outputs/ISR_LCBD_corr_plot.pdf", plot = dd, width = 12, height = 10)
  
  return(NULL)  
}


# 
# tab <- as.data.frame(rbind(runarms_index_plate))
# # tab$site <- rownames(tab)
# 
# #### Import the GPS data ####
# data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE)
# 
# #### Import the map shapefiles ####
# runa_map <- st_read(runa_map)
# roda_map <- st_read(roda_map)
# 
# 
# 
# data_gps <- data_gps %>%
#   mutate(
#     label = case_when(
#       grepl("P50ARMS", Site) ~ gsub("P50ARMS", "P50A", Site),
#       grepl("RUNARMS", Site) ~ gsub("RUNARMS", "RUNA", Site),
#       grepl("RODARMS", Site) ~ gsub("RODARMS", "RODA", Site),
#       TRUE ~ Site
#     )
#   )
# 
# data_gps <- data_gps %>%
#   mutate(
#     site = case_when(
#       grepl("RUNA", Site) ~ gsub("RUNA", "RUNARMS", Site),
#       grepl("RODA", Site) ~ gsub("RODA", "RODARMS", Site),
#       TRUE ~ Site
#     )
#   )
# 
# map_data <- left_join(tab, data_gps, by = "site")
# 
# map_data_runa <- map_data[grepl("RUNA", map_data$label),]
# 
# # Convert to sf object
# map_data_runa_sf <- st_as_sf(map_data_runa, coords = c("Longitude", "Latitude"), crs = 4326)
# st_crs(map_data_runa_sf) <- 4326
# 
# runa_map <- st_transform(runa_map, crs = 4326)
# 
# map_data_runa_sf <- st_transform(map_data_runa_sf, st_crs(runa_map))
# st_crs(runa_map)
# st_crs(map_data_runa_sf)
# 
# map_data_runa_sf <- map_data_runa_sf %>%
#   mutate(Longitude = st_coordinates(.)[,1],
#          Latitude = st_coordinates(.)[,2])
# 
# lcbd_all <- c(map_data$LCBD_all)
# min_LCBD <- min(lcbd_all)
# max_LCBD <- max(lcbd_all)
# 
# scale_lcbd <- scale_size_continuous(name = "LCBD", range = c(2, 10), limits = c(min_LCBD, max_LCBD))
# 
# #### Import the reef shapefiles ####
# 
# runa_reef <- st_read(runa_reef)
# roda_reef <- st_read(roda_reef)
# 
# runa_reef <- st_transform(runa_reef, crs = 4326)
# roda_reef <- st_transform(roda_reef, crs = 4326)
# 
# # #### RUNA map settings ####
# # map_data_p50a <- map_data[grepl("P50A", map_data$Site),]
# # map_data_p50a_sf <- st_as_sf(map_data_p50a, coords = c("Longitude", "Latitude"), crs = 4326)
# # st_crs(map_data_p50a_sf) <- 4326
# # map_data_p50a_sf <- st_transform(map_data_p50a_sf, st_crs(runa_map))
# # st_crs(runa_map)
# # st_crs(map_data_p50a_sf)
# # map_data_p50a_sf <- map_data_p50a_sf %>%
# #   mutate(Longitude = st_coordinates(.)[,1],
# #          Latitude = st_coordinates(.)[,2])
# # 
# # #### RODA map settings ####
# # map_data_roda <- map_data[grepl("RODA", map_data$Site),]
# # map_data_roda_sf <- st_as_sf(map_data_roda, coords = c("Longitude", "Latitude"), crs = 4326)
# # st_crs(roda_map)
# # st_crs(map_data_roda_sf)
# # map_data_roda_sf <- map_data_roda_sf %>%
# #   mutate(Longitude = st_coordinates(.)[,1],
# #          Latitude = st_coordinates(.)[,2])
# 
# 
# # -- Créer l'objet sf principal --
# map_data_sf <- st_as_sf(map_data, coords = c("Longitude", "Latitude"), crs = 4326)
# map_data_sf <- st_transform(map_data_sf, st_crs(runa_map)) %>%
#   mutate(Longitude = st_coordinates(.)[,1],
#          Latitude = st_coordinates(.)[,2])
# 
# # -- Ajuster les coordonnées pour éviter les chevauchements --
# map_data_sf <- map_data_sf %>%
#   mutate(
#     Latitude = case_when(
#       site %in% c("RUNARMS2") ~ Latitude + 0.006,
#       site %in% c("RUNARMS3") ~ Latitude - 0.006,
#       site %in% c("RUNARMS4") ~ Latitude + 0.005,
#       site %in% c("RUNARMS5") ~ Latitude - 0.005,
#       site %in% c("RUNARMS6") ~ Latitude - 0.005,
#       TRUE ~ Latitude
#     )
#   )
# 
# map_data_sf <- map_data_sf %>%
#   mutate(
#     Longitude = case_when(
#       site %in% c("RUNARMS2") ~ Longitude - 0.005,
#       site %in% c("RUNARMS3") ~ Longitude + 0.007,
#       TRUE ~ Longitude
#     )
#   )
# 
# # -- Scales communes --
# scale_r <- scale_size_continuous(
#   name = "ISR", 
#   range = c(1, 12),
#   limits = c(min(map_data$R_all, na.rm = TRUE), max(map_data$R_all, na.rm = TRUE))
# )
# 
# 
# 
# 
# scale_lcbd <- scale_fill_gradientn(
#   colours = rev(viridis::viridis(10)),
#   name = "LCBD",
#   limits = c(min(map_data$LCBD_all, na.rm = TRUE), max(map_data$LCBD_all, na.rm = TRUE)),
#   guide = guide_colorbar(
#     title.hjust = 0.5,
#     barwidth = unit(4, "cm"),     # allonge la barre
#     barheight = unit(0.4, "cm"),  # réduit l’épaisseur (utile pour une légende horizontale)
#     label.theme = element_text(size = 8)  # taille des chiffres
#   )
# )
# 
# # Define common width and heights
# map_width <- 0.45
# map_height <- 0.45
# 
# # Center coordinates for each map
# center_p50a <- c(x = mean(c(55.15, 55.6)), y = mean(c(-21.4, -21)))
# center_roda <- c(x = mean(c(63.20, 63.58)), y = mean(c(-19.90, -19.58)))
# # (Add RUNA similarly if needed)
# 
# # Create coordinate limits
# p50a_xlim <- c(center_p50a['x'] - map_width / 2, center_p50a['x'] + map_width / 2)
# p50a_ylim <- c(center_p50a['y'] - map_height / 2, center_p50a['y'] + map_height / 2)
# 
# roda_xlim <- c(center_roda['x'] - map_width / 2, center_roda['x'] + map_width / 2)
# roda_ylim <- c(center_roda['y'] - map_height / 2, center_roda['y'] + map_height / 2)
# 
# # Then apply in each coord_sf:
# coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE)
# coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE)
# 
# # --- 9. Création des cartes ---
# runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
# p50a_sf <- map_data_sf[grepl("P50ARMS", map_data_sf$site),]
# roda_sf <- map_data_sf[grepl("RODARMS", map_data_sf$site),]
# 
# # -- Map the RUNA data --
# map_data_runa_sf <- map_data_sf[grepl("RUNARMS", map_data_sf$site),]
# 
# 
# uu <- ggplot() +
#   geom_sf(data = runa_map, fill = "grey85", color = "black") +
#   geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
#   geom_point(data = map_data_runa_sf,
#              aes(x = Longitude, y = Latitude,size = R_all, fill = LCBD_all),
#              shape = 21, stroke = 0.6) +
#   geom_text_repel(
#     data = map_data_runa_sf,
#     aes(x = Longitude, y = Latitude, label = label),
#     size = 4,
#     fontface = "bold",
#     color = "black",
#     box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
#     max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
#     force = 1.5,               # Force de répulsion
#     segment.size = 0.2,        # Taille du segment
#     segment.color = "grey30"   # Couleur du trait
#   ) +
#   # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
#   coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
#   theme_minimal() +
#   labs(title = "RUNA") +
#   scale_r + scale_lcbd +
#   annotation_scale(location = "bl", width_hint = 0.25) +
#   annotation_north_arrow(location = "tr", which_north = "true",
#                          pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
#                          style = north_arrow_fancy_orienteering)
# 
# 
# zuzu <- (uu) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# 
# ggsave("outputs/Reunion_plot.pdf", plot = zuzu, width = 6, height = 5)
