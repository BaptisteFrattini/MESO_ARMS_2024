# ' plot values of LCBD and rarity on map by campain
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_map_final <- function(data_and_meta_clean_fullsites, gps_sites, runa_map, roda_map, runa_reef, roda_reef){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites")
  # gps_sites = targets::tar_read("data_gps_sites")
  # runa_map = targets::tar_read("map_runa")
  # roda_map = targets::tar_read("map_roda")
  # runa_reef = targets::tar_read("reef_runa")
  # roda_reef = targets::tar_read("reef_roda")

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
  
  # All taxa map ####
  ## Computing indices ####
  ### Singularity/uniqueness/rarity index ####
  
  data_LCBD <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_LCBD <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
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
                                               "Erect_Rhodophyta_algae")]
  
  data_filtered <- data[, msp_list_filter]
  data_filtered_pa <- vegan::decostand(data_filtered, "pa")
  
  
  species_names <- data.frame(Species = colnames(data_filtered))
  
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
  
  campagne_name = "P50ARMS"
  # 2. Fonction pour calculer recap_rarity par campagne
  compute_recap_rarity_per_campagne <- function(campagne_name) {
    
    # Filtrer les données pour la campagne
    
    data_campagne <- subset(data_filtered_pa, meta$campain == campagne_name)
    meta_campagne <- subset(meta, meta$campain == campagne_name)
    
    # Ascidiacea
    asc_species <- corr_taxa$Species[corr_taxa$taxa == "Ascidiacea"]
    data_asc <- data_campagne[, colnames(data_campagne) %in% asc_species]
    occ_asc <- colSums(data_asc)
    occ_asc <- occ_asc[occ_asc != 0] / nrow(data_campagne)
    threshold_asc <- quantile(occ_asc, 0.25)
    
    asc_species <- asc_species[colSums(data_asc) > 0]
    data_asc$triplicat <- meta_campagne$triplicat
    
    result <- sapply(asc_species, function(espece) {
      subset <- data_asc[data_asc[[espece]] == 1, ]
      n_triplicats <- length(unique(subset$triplicat))
      if (n_triplicats <= 1) {
        return("Un seul triplicat")
      } else {
        return("Plusieurs triplicats")
      }
    })
    
    out_asc <- data.frame(Species = names(occ_asc),
                          f = occ_asc, taxa = "Ascidiacea",
                          threshold = threshold_asc,
                          repartition = result
                          )
    
    # Porifera
    por_species <- corr_taxa$Species[corr_taxa$taxa == "Porifera"]
    data_por <- data_campagne[, colnames(data_campagne) %in% por_species]
    occ_por <- colSums(data_por)
    occ_por <- occ_por[occ_por != 0] / nrow(data_campagne)
    threshold_por <- quantile(occ_por, 0.25)
    
    por_species <- por_species[colSums(data_por) > 0]
    data_por$triplicat <- meta_campagne$triplicat
    
    result <- sapply(por_species, function(espece) {
      subset <- data_por[data_por[[espece]] == 1, ]
      n_triplicats <- length(unique(subset$triplicat))
      if (n_triplicats <= 1) {
        return("Un seul triplicat")
      } else {
        return("Plusieurs triplicats")
      }
    })
    
    out_por <- data.frame(Species = names(occ_por),
                          f = occ_por, taxa = "Porifera",
                          threshold = threshold_por,
                          repartition = result
    )
    
    # Bryozoa
    bry_species <- corr_taxa$Species[corr_taxa$taxa == "Bryozoa"]
    data_bry <- data_campagne[, colnames(data_campagne) %in% bry_species]
    occ_bry <- colSums(data_bry)
    occ_bry <- occ_bry[occ_bry != 0] / nrow(data_campagne)
    threshold_bry <- quantile(occ_bry, 0.25)
    
    bry_species <- bry_species[colSums(data_bry) > 0]
    data_bry$triplicat <- meta_campagne$triplicat
    
    result <- sapply(bry_species, function(espece) {
      subset <- data_bry[data_bry[[espece]] == 1, ]
      n_triplicats <- length(unique(subset$triplicat))
      if (n_triplicats <= 1) {
        return("Un seul triplicat")
      } else {
        return("Plusieurs triplicats")
      }
    })
    
    out_bry <- data.frame(Species = names(occ_bry),
                          f = occ_bry, taxa = "Bryozoa",
                          threshold = threshold_bry,
                          repartition = result
    )
    
    # Autres (We compute a threshold based on every msp to use with the other)
    oth_species <- corr_taxa$Species
    data_oth <- data_campagne[, colnames(data_campagne) %in% oth_species]
    occ_oth <- colSums(data_oth)
    occ_oth <- occ_oth[occ_oth != 0] / nrow(data_campagne)
    threshold_oth <- quantile(occ_oth, 0.25)
    
    oth_species <- oth_species[colSums(data_oth) > 0]
    data_oth$triplicat <- meta_campagne$triplicat
    
    result <- sapply(oth_species, function(espece) {
      subset <- data_oth[data_oth[[espece]] == 1, ]
      n_triplicats <- length(unique(subset$triplicat))
      if (n_triplicats <= 1) {
        return("Un seul triplicat")
      } else {
        return("Plusieurs triplicats")
      }
    })
    
    out_oth <- data.frame(Species = names(occ_oth),
                          f = occ_oth, taxa = "Other",
                          threshold = threshold_oth,
                          repartition = result
    )
    # Fusionner
    recap <- bind_rows(out_asc, out_por, out_bry, out_oth) 
    
    # Créer la nouvelle colonne "rarete"
    recap$rarity <- ifelse(recap$f <= recap$threshold & recap$repartition == "Un seul triplicat", "rare", "common")
    
    # Ajouter campagne
    recap$campagne <- campagne_name
    return(recap)
  }
  
  # 3. Appliquer pour chaque campagne
  campagnes <- unique(meta$campain)
  
  recap_rarity_all <- lapply(campagnes, compute_recap_rarity_per_campagne) %>%
    bind_rows()
  
  
  # 4. Optionnel : split en 3 objets pour RUNA, P50A, RODA
  recap_rarity_runa <- recap_rarity_all %>% filter(campagne == "RUNARMS")
  recap_rarity_p50a <- recap_rarity_all %>% filter(campagne == "P50ARMS")
  recap_rarity_roda <- recap_rarity_all %>% filter(campagne == "RODARMS")
  
  # Espèces rares
  rare_species_runa <- recap_rarity_runa %>%
    filter(rarity == "rare") %>%
    pull(Species)
  
  rare_species_p50a <- recap_rarity_p50a %>%
    filter(rarity == "rare") %>%
    pull(Species)
  
  rare_species_roda <- recap_rarity_roda %>%
    filter(rarity == "rare") %>%
    pull(Species)
  
  
  # Richness
  
  richesse <- specnumber(data_filtered_pa, groups = meta$triplicat)
  
  # Calcul du nombre d'espèces rares par triplicat
  
  data_with_meta <- cbind(data_filtered_pa, meta)
  
  data_with_meta_runa <- subset(data_with_meta, campain ==  "RUNARMS")
  
  nb_rare_df_runa <- data_with_meta_runa %>%
    select(all_of(rare_species_runa), triplicat) %>%
    group_by(triplicat) %>%
    summarise(nb_rare = sum(colSums(across(everything())) > 0), .groups = "drop")
  
  data_with_meta_p50a <- subset(data_with_meta, campain ==  "P50ARMS")
  
  nb_rare_df_p50a <- data_with_meta_p50a %>%
    select(all_of(rare_species_p50a), triplicat) %>%
    group_by(triplicat) %>%
    summarise(nb_rare = sum(colSums(across(everything())) > 0), .groups = "drop")
  
  data_with_meta_roda <- subset(data_with_meta, campain ==  "RODARMS")
  
  nb_rare_df_roda <- data_with_meta_roda %>%
    select(all_of(rare_species_roda), triplicat) %>%
    group_by(triplicat) %>%
    summarise(nb_rare = sum(colSums(across(everything())) > 0), .groups = "drop")
  
  nb_rare_df <- data.frame(rbind(nb_rare_df_runa, nb_rare_df_p50a, nb_rare_df_roda))
  
  # Mise en tableau final
  indice_rareté <- data.frame(
    triplicat = names(richesse),
    richesse = richesse
  ) %>%
    left_join(nb_rare_df, by = "triplicat") %>%
    mutate(R = nb_rare / richesse)
  # mutate(R = nb_rare)
  
  #### Compute LCBD ####
  
  
  matrix.bray <- vegan::vegdist(data_LCBD, method = "bray")
  spe.beta <- adespatial::LCBD.comp( matrix.bray, sqrt.D = TRUE)
  spe.beta$LCBD
  names(spe.beta$LCBD) <- rownames(data_LCBD)
  
  
  
  # Aggregate beta diversity by site (mean of the replicates)
  
  spe.beta$Site <- substr(names(spe.beta$LCBD), 1, 8)  # Extract main site name
  spe.beta <- data.frame(site = spe.beta$Site, LCBD = spe.beta$LCBD) 
  beta_agg <- spe.beta %>%
    group_by(site) %>%
    summarise(LCBD = mean(LCBD))
  
  #### Import the GPS data ####
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE)
  
  # Merge with site coordinates
  map_data <- as.data.frame(cbind(data_gps, beta_agg))
  
  
  
  ## Mapping the indices ####
  
  #### Import the map shapefiles ####
  runa_map <- st_read(runa_map)
  roda_map <- st_read(roda_map)
  
  map_data_runa <- map_data[grepl("RUNA", map_data$Site),]
  
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
  
  
  
  lcbd_all <- c(map_data$LCBD)
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
  
  
  
  # -- Harmoniser les noms de colonnes
  indice_rareté <- indice_rareté %>% rename(site = triplicat)
  
  # -- Modifier les noms de sites dans data_gps
  data_gps$site <- gsub("^RODA", "RODARMS", data_gps$Site)
  data_gps$site <- gsub("^RUNA", "RUNARMS", data_gps$site)
  
  # -- Creer le tableau des informations à génerer sur la carte
  rarity_lcbd_map_data <- left_join(indice_rareté, beta_agg, by = "site") %>%
    left_join(data_gps, by = "site")
  
  # -- Créer l'objet sf principal --
  map_data_sf <- st_as_sf(rarity_lcbd_map_data, coords = c("Longitude", "Latitude"), crs = 4326)
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
  # -- Renommer certain site
  map_data_sf <- map_data_sf %>%
    mutate(
      label = case_when(
        grepl("P50ARMS", site) ~ gsub("P50ARMS", "P50A", site),
        grepl("RUNARMS", site) ~ gsub("RUNARMS", "RUNA", site),
        grepl("RODARMS", site) ~ gsub("RODARMS", "RODA", site),
        TRUE ~ site
      )
    )
  
  
  # -- Scales communes --
  scale_r <- scale_size_continuous(
    name = "Indice de rareté (R)", 
    range = c(2, 10),
    limits = c(min(rarity_lcbd_map_data$R, na.rm = TRUE), max(rarity_lcbd_map_data$R, na.rm = TRUE))
  )
  
  scale_lcbd <- scale_fill_gradientn(
    colours = rev(viridis::viridis(10)),
    name = "LCBD",
    limits = c(min(rarity_lcbd_map_data$LCBD, na.rm = TRUE), max(rarity_lcbd_map_data$LCBD, na.rm = TRUE)),
    guide = guide_colorbar(
      title.hjust = 0.5,
      barwidth = unit(4, "cm"),     # allonge la barre
      barheight = unit(0.4, "cm"),  # réduit l’épaisseur (utile pour une légende horizontale)
      label.theme = element_text(size = 8)  # taille des chiffres
    )
  )
  
  # Define common width and height
  map_width <- 0.45
  map_height <- 0.4
  
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
               aes(x = Longitude, y = Latitude, size = R, fill = LCBD),
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
               aes(x = Longitude, y = Latitude, size = R, fill = LCBD),
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
               aes(x = Longitude, y = Latitude, size = R, fill = LCBD),
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
  
  
  # ########################################################################## #
  
  # Separate taxa maps ####
  # group_name <- "Ascidiacea"
  # rare_species_df <- recap_rarity_all
  
  make_maps_for_group <- function(group_name, rare_species_df) {
    # --- 1. Sélection des espèces du groupe ---
    group_species <- corr_taxa %>%
      filter(taxa == group_name) %>%
      pull(Species)
    
    # --- 2. Recalculer le LCBD pour le groupe uniquement ---
    data_group <- data_LCBD[, colnames(data_LCBD) %in% group_species]
    # data_mean_group <- vegan::decostand(data_mean_group, "pa")
    matrix.bray <- vegan::vegdist(data_group, method = "bray")
    spe.beta <- adespatial::LCBD.comp(matrix.bray, sqrt.D = TRUE)
    names(spe.beta$LCBD) <- rownames(data_LCBD)
    spe.beta$Site <- substr(names(spe.beta$LCBD), 1, 8)
    beta_agg <- data.frame(site = spe.beta$Site, LCBD = spe.beta$LCBD) %>%
      group_by(site) %>%
      summarise(LCBD = mean(LCBD), .groups = "drop")
    
    # --- 3. Sélection des espèces rares du groupe ---
    rare_species_runa <- rare_species_df %>%
      filter(rarity == "rare", taxa == group_name, campagne == "RUNARMS") %>%
      pull(Species)
    rare_species_p50a <- rare_species_df %>%
      filter(rarity == "rare", taxa == group_name, campagne == "P50ARMS") %>%
      pull(Species)
    rare_species_roda <- rare_species_df %>%
      filter(rarity == "rare", taxa == group_name, campagne == "RODARMS") %>%
      pull(Species)
    
    # --- 4. Calcul richesse spécifique du groupe ---
    mat_group <- data_with_meta[, colnames(data_with_meta) %in% group_species]
    richesse <- specnumber(mat_group, groups = data_with_meta$triplicat)
    
    # --- 5. Calcul du nombre d'espèces rares présentes par site ---
    data_with_meta_runa <- subset(data_with_meta, campain ==  "RUNARMS")
    
    nb_rare_df_runa <- data_with_meta_runa %>%
      select(all_of(rare_species_runa), triplicat) %>%
      group_by(triplicat) %>%
      summarise(nb_rare = sum(colSums(across(everything())) > 0), .groups = "drop")
    
    data_with_meta_p50a <- subset(data_with_meta, campain ==  "P50ARMS")
    
    nb_rare_df_p50a <- data_with_meta_p50a %>%
      select(all_of(rare_species_p50a), triplicat) %>%
      group_by(triplicat) %>%
      summarise(nb_rare = sum(colSums(across(everything())) > 0), .groups = "drop")
    
    data_with_meta_roda <- subset(data_with_meta, campain ==  "RODARMS")
    
    nb_rare_df_roda <- data_with_meta_roda %>%
      select(all_of(rare_species_roda), triplicat) %>%
      group_by(triplicat) %>%
      summarise(nb_rare = sum(colSums(across(everything())) > 0), .groups = "drop")
    
    nb_rare_df <- data.frame(rbind(nb_rare_df_runa, nb_rare_df_p50a, nb_rare_df_roda))
    
    # --- 6. Calcul de l'indice R ---
    indice_rareté <- data.frame(
      triplicat = names(richesse),
      richesse = richesse
    ) %>%
      left_join(nb_rare_df, by = "triplicat") %>%
      mutate(R = nb_rare / richesse,
             site = triplicat)
    # mutate(R = nb_rare,
    #      site = triplicat)
    
    # --- 7. Préparation des données de cartographie ---
    rarity_lcbd_map_data <- indice_rareté %>%
      left_join(beta_agg, by = "site") %>%
      left_join(data_gps, by = "site")
    
    rarity_lcbd_map_data <- rarity_lcbd_map_data[,-7]
    
    
    map_data_sf <- st_as_sf(rarity_lcbd_map_data, coords = c("Longitude", "Latitude"), crs = 4326) %>%
      st_transform(st_crs(runa_map)) %>%
      mutate(
        Longitude = st_coordinates(.)[,1],
        Latitude = st_coordinates(.)[,2],
        Latitude = case_when(
          site %in% c("RUNARMS2") ~ Latitude + 0.006,
          site %in% c("RUNARMS3") ~ Latitude - 0.006,
          site %in% c("RUNARMS4") ~ Latitude + 0.005,
          site %in% c("RUNARMS5") ~ Latitude - 0.005,
          site %in% c("RUNARMS6") ~ Latitude - 0.005,
          TRUE ~ Latitude
        ),
        Longitude = case_when(
          site %in% c("RUNARMS2") ~ Longitude - 0.005,
          site %in% c("RUNARMS3") ~ Longitude + 0.007,
          TRUE ~ Longitude
        ),
        label = case_when(
          grepl("P50ARMS", site) ~ gsub("P50ARMS", "P50A", site),
          grepl("RUNARMS", site) ~ gsub("RUNARMS", "RUNA", site),
          grepl("RODARMS", site) ~ gsub("RODARMS", "RODA", site),
          TRUE ~ site
        )
      )
    
    # --- 8. Scales ---
    scale_r <- scale_size_continuous(
      name = "Indice de rareté (R)", 
      range = c(2, 10),
      limits = c(min(rarity_lcbd_map_data$R, na.rm = TRUE), max(rarity_lcbd_map_data$R, na.rm = TRUE))
    )
    
    scale_lcbd <- scale_fill_gradientn(
      colours = rev(viridis::viridis(10)),
      name = "LCBD",
      limits = c(min(rarity_lcbd_map_data$LCBD, na.rm = TRUE), max(rarity_lcbd_map_data$LCBD, na.rm = TRUE)),
      guide = guide_colorbar(
        title.hjust = 0.5,
        barwidth = unit(4, "cm"),     # allonge la barre
        barheight = unit(0.4, "cm"),  # réduit l’épaisseur (utile pour une légende horizontale)
        label.theme = element_text(size = 8)  # taille des chiffres
      )
    )
    
    # Define common width and height
    map_width <- 0.45
    map_height <- 0.4
    
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
    
    # RUNA
    aa <- ggplot() +
      geom_sf(data = runa_map, fill = "grey85", color = "black") +
      geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
      geom_point(data = runa_sf,
                 aes(x = Longitude, y = Latitude, size = R, fill = LCBD),
                 shape = 21, stroke = 0.6) +
      geom_text_repel(data = runa_sf,
                      aes(x = Longitude, y = Latitude, label = label),
                      size = 4, fontface = "bold", color = "black",
                      box.padding = 0.8, force = 1.5,
                      segment.size = 0.2, segment.color = "grey30",
                      max.overlaps = Inf) +
      # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
      coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
      theme_minimal() +
      labs(title = paste("RUNA –", group_name)) +
      scale_lcbd + scale_r +
      annotation_scale(location = "bl", width_hint = 0.25) +
      annotation_north_arrow(location = "tr", which_north = "true",
                             pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                             style = north_arrow_fancy_orienteering)
    
    # P50A
    bb <- ggplot() +
      geom_sf(data = runa_map, fill = "grey85", color = "black") +
      geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
      geom_point(data = p50a_sf,
                 aes(x = Longitude, y = Latitude, size = R, fill = LCBD),
                 shape = 21, stroke = 0.6, color = "black") +
      geom_text_repel(data = p50a_sf,
                      aes(x = Longitude, y = Latitude, label = label),
                      size = 4, fontface = "bold", color = "black",
                      box.padding = 0.8, force = 1.5,
                      segment.size = 0.2, segment.color = "grey30",
                      max.overlaps = Inf) +
      # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
      coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
      theme_minimal() +
      labs(title = paste("P50A –", group_name)) +
      scale_lcbd + scale_r +
      annotation_scale(location = "bl", width_hint = 0.25) +
      annotation_north_arrow(location = "tr", which_north = "true",
                             pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                             style = north_arrow_fancy_orienteering)
    
    # RODA
    cc <- ggplot() +
      geom_sf(data = roda_map, fill = "grey85", color = "black") +
      geom_sf(data = roda_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
      geom_point(data = roda_sf,
                 aes(x = Longitude, y = Latitude, size = R, fill = LCBD),
                 shape = 21, stroke = 0.6, color = "black") +
      geom_text_repel(data = roda_sf,
                      aes(x = Longitude, y = Latitude, label = label),
                      size = 4, fontface = "bold", color = "black",
                      box.padding = 0.8, force = 1.5,
                      segment.size = 0.2, segment.color = "grey30",
                      max.overlaps = Inf) +
      # coord_sf(xlim = c(63.25, 63.55), ylim = c(-19.85, -19.62), expand = FALSE) +
      coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE) +
      theme_minimal() +
      labs(title = paste("RODA –", group_name)) +
      scale_lcbd + scale_r +
      annotation_scale(location = "bl", width_hint = 0.25) +
      annotation_north_arrow(location = "tr", which_north = "true",
                             pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                             style = north_arrow_fancy_orienteering)
    
    # --- 10. Combinaison finale ---
    return((aa | bb | cc) + plot_layout(guides = "collect") & theme(legend.position = "bottom"))
  }
  
  ### Make the maps ####
  
  xx <- make_maps_for_group("Ascidiacea", recap_rarity_all)
  yy <- make_maps_for_group("Porifera", recap_rarity_all)
  zz <-make_maps_for_group("Bryozoa", recap_rarity_all)
  
  
  final_plot <- (tt / xx / yy / zz) +
    plot_layout(guides = "collect")
  
  ggsave("outputs/Cartes - LCBD_Rarity/final_map.pdf", plot = final_plot, width = 0.75*18, height = 18)
  
  
  
  return(NULL) 
}
