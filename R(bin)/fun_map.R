# ' plot values of LCBD and rarity on map
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_map <- function(data_and_meta_clean_fullsites, gps_sites, runa_map, roda_map, runa_reef, roda_reef){
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

  #### Compute Rarity threshold ####  
  #Ascidians
  ascidiacea_species <- corr_taxa$Species[corr_taxa$taxa == "Ascidiacea"]
  data_ascidiacea <- data_filtered_pa[, colnames(data_filtered_pa) %in% ascidiacea_species]
  
  species_occurence <- colSums(data_ascidiacea)/nrow(data_ascidiacea)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_asc <- quantile(species_occurence$f, 0.25)
  
  #Porifera
  porifera_species <- corr_taxa$Species[corr_taxa$taxa == "Porifera"]
  data_porifera <- data_filtered_pa[, colnames(data_filtered_pa) %in% porifera_species]
  
  species_occurence <- colSums(data_porifera)/nrow(data_porifera)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_por <- quantile(species_occurence$f, 0.25)
  
  #Bryozoa
  bryozoa_species <- corr_taxa$Species[corr_taxa$taxa == "Bryozoa"]
  data_bryozoa <- data_filtered_pa[, colnames(data_filtered_pa) %in% bryozoa_species]
  
  species_occurence <- colSums(data_bryozoa)/nrow(data_bryozoa)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_bry <- quantile(species_occurence$f, 0.25)
  
  #all_species
  species_occurence <- colSums(data_filtered_pa)/nrow(data_filtered_pa)
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold_all_species <- quantile(species_occurence$f, 0.25)
  
  
  # Fonction pour calculer la rareté pour un groupe
  compute_rarity <- function(species_list, group_name, threshold) {
    df <- data_filtered_pa[, colnames(data_filtered_pa) %in% species_list]
    freq <- colSums(df) / nrow(df)
    freq <- freq[freq != 0]
    out <- data.frame(Species = names(freq), f = freq)
    out$taxa <- group_name
    out$threshold <- threshold
    out$rarity <- ifelse(out$f <= threshold, "rare", "common")
    return(out)
  }
  
  # Par groupe avec seuils spécifiques
  out_asc <- compute_rarity(ascidiacea_species, "Ascidiacea", threshold_asc)
  out_por <- compute_rarity(porifera_species, "Porifera", threshold_por)
  out_bry <- compute_rarity(bryozoa_species, "Bryozoa", threshold_bry)
  
  # Pour les autres groupes → on applique le seuil général
  autres_groupes <- corr_taxa %>%
    filter(!taxa %in% c("Ascidiacea", "Porifera", "Bryozoa")) %>%
    pull(Species)
  
  out_others <- compute_rarity(autres_groupes, "Other", threshold_all_species)
  
  # Combiner tous les résultats
  recap_rarity <- bind_rows(out_asc, out_por, out_bry, out_others)
  
  recap_rarity <- left_join(recap_rarity, corr_taxa, by = "Species") %>%
    mutate(taxa = ifelse(taxa.x == "Other", taxa.y, taxa.x)) %>%
    select(Species, taxa, f, threshold, rarity)
  
  # S'assurer que les noms de lignes sont bien les mêmes
  data_filtered_pa$name <- rownames(data_filtered_pa)
  
  # Fusion avec les métadonnées
  data_with_meta <- left_join(data_filtered_pa, meta, by = "name")

  # Espèces rares
  rare_species <- recap_rarity %>%
    filter(rarity == "rare") %>%
    pull(Species)
  
  # Ajouter le nom des échantillons
  data_filtered_pa$name <- rownames(data_filtered_pa)
  
  # Fusion avec les métadonnées
  data_with_meta <- left_join(data_filtered_pa, meta, by = "name")
  
  # Liste des colonnes espèces uniquement
  species_cols <- colnames(data_filtered_pa)[!colnames(data_filtered_pa) %in% c("name")]
  
  # Calcul avec vegan::specnumber()
  # On commence par isoler la matrice espèces seule
  mat_species <- data_with_meta[, species_cols]
  
  # Calcul de la richesse par triplicat
  richesse <- specnumber(mat_species, groups = data_with_meta$triplicat)
  
  # Calcul du nombre d'espèces rares par triplicat
  nb_rare_df <- data_with_meta %>%
    select(all_of(rare_species), triplicat) %>%
    group_by(triplicat) %>%
    summarise(nb_rare = sum(colSums(across(everything())) > 0), .groups = "drop")
  
  # Mise en tableau final
  indice_rareté <- data.frame(
    triplicat = names(richesse),
    richesse = richesse
  ) %>%
    left_join(nb_rare_df, by = "triplicat") %>%
    mutate(R = nb_rare / richesse)
  
  #### Compute LCBD ####
  
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  data_mean_filtered <- data_mean[, msp_list_filter]
  
  matrix.jacc <- vegan::vegdist(data_mean_filtered, method = "jaccard")
  spe.beta <- adespatial::LCBD.comp( matrix.jacc, sqrt.D = TRUE)

  spe.beta$LCBD
  
  names(spe.beta$LCBD) <- rownames(data_mean_filtered)
  
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE)
  
  
  # Read the shapefile
  runa_map <- st_read(runa_map)
  roda_map <- st_read(roda_map)
  
  # Aggregate beta diversity by site (mean of the replicates)
  
  spe.beta$Site <- substr(names(spe.beta$LCBD), 1, 8)  # Extract main site name
  spe.beta <- data.frame(site = spe.beta$Site, LCBD = spe.beta$LCBD) 
  beta_agg <- spe.beta %>%
    group_by(site) %>%
    summarise(LCBD = mean(LCBD))
  
  # Merge with site coordinates
  map_data <- as.data.frame(cbind(data_gps, beta_agg))
  
  ## Mapping the indices ####
  
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
  
  # importer les reef
  
  runa_reef <- st_read(runa_reef)
  roda_reef <- st_read(roda_reef)
  
  runa_reef <- st_transform(runa_reef, crs = 4326)
  roda_reef <- st_transform(roda_reef, crs = 4326)
  
  map_data_p50a <- map_data[grepl("P50A", map_data$Site),]
  map_data_p50a_sf <- st_as_sf(map_data_p50a, coords = c("Longitude", "Latitude"), crs = 4326)
  st_crs(map_data_p50a_sf) <- 4326  
  map_data_p50a_sf <- st_transform(map_data_p50a_sf, st_crs(runa_map))
  st_crs(runa_map)
  st_crs(map_data_p50a_sf)
  map_data_p50a_sf <- map_data_p50a_sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2])
  
  map_data_roda <- map_data[grepl("RODA", map_data$Site),]
  map_data_roda_sf <- st_as_sf(map_data_roda, coords = c("Longitude", "Latitude"), crs = 4326)
  st_crs(roda_map)
  st_crs(map_data_roda_sf)
  map_data_roda_sf <- map_data_roda_sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2])
  
  
  # Harmoniser les noms de colonnes
  indice_rareté <- indice_rareté %>% rename(site = triplicat)
  # Modifier les noms de sites dans data_gps
  data_gps$site <- gsub("^RODA", "RODARMS", data_gps$Site)
  data_gps$site <- gsub("^RUNA", "RUNARMS", data_gps$site)
  rarity_lcbd_map_data <- left_join(indice_rareté, beta_agg, by = "site") %>%
    left_join(data_gps, by = "site")
  
  # -- Créer l'objet sf principal --
  map_data_sf <- st_as_sf(rarity_lcbd_map_data, coords = c("Longitude", "Latitude"), crs = 4326)
  map_data_sf <- st_transform(map_data_sf, st_crs(runa_map)) %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2])
  
  
  # -- Ajuster les coordonnées pour éviter les chevauchements --

  
  # Modifier les noms de sites dans data_gps
  data_gps$site <- gsub("^RODA", "RODARMS", data_gps$Site)
  data_gps$site <- gsub("^RUNA", "RUNARMS", data_gps$site)
  
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
    colours = c("darkblue", "coral"),
    name = "LCBD",
    limits = c(min(rarity_lcbd_map_data$LCBD, na.rm = TRUE), max(rarity_lcbd_map_data$LCBD, na.rm = TRUE))
  )
  
  # -- RUNA --
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
    coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    theme_minimal() +
    labs(title = "RUNA") +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- P50A --
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
    coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    theme_minimal() +
    labs(title = "P50A") +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- RODA --
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
    coord_sf(xlim = c(63.25, 63.55), ylim = c(-19.85, -19.62), expand = FALSE) +
    theme_minimal() +
    labs(title = "RODA") +
    scale_r + scale_lcbd +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  # -- Combinaison finale avec guides partagés --

  
  tt <- (uu | vv | ww) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  
# Separate taxa maps ####
 
  group_name <- "Ascidiacea"
  rare_species_df <- recap_rarity
  
  make_maps_for_group <- function(group_name, rare_species_df) {
    # --- 1. Sélection des espèces du groupe ---
    group_species <- corr_taxa %>%
      filter(taxa == group_name) %>%
      pull(Species)
    
    # --- 2. Recalculer le LCBD pour le groupe uniquement ---
    data_mean_group <- data_mean[, colnames(data_mean) %in% group_species]
    data_mean_group <- vegan::decostand(data_mean_group, "pa")
    matrix.jacc <- vegan::vegdist(data_mean_group, method = "jaccard")
    spe.beta <- adespatial::LCBD.comp(matrix.jacc, sqrt.D = TRUE)
    names(spe.beta$LCBD) <- rownames(data_mean_group)
    spe.beta$Site <- substr(names(spe.beta$LCBD), 1, 8)
    beta_agg <- data.frame(site = spe.beta$Site, LCBD = spe.beta$LCBD) %>%
      group_by(site) %>%
      summarise(LCBD = mean(LCBD), .groups = "drop")
    
    # --- 3. Sélection des espèces rares du groupe ---
    rare_species <- rare_species_df %>%
      filter(rarity == "rare", taxa == group_name) %>%
      pull(Species)
    
    # --- 4. Calcul richesse spécifique du groupe ---
    mat_group <- data_with_meta[, colnames(data_with_meta) %in% group_species]
    richesse <- specnumber(mat_group, groups = data_with_meta$triplicat)
    
    # --- 5. Calcul du nombre d'espèces rares présentes par site ---
    nb_rare_df <- data_with_meta %>%
      select(all_of(rare_species), triplicat) %>%
      group_by(triplicat) %>%
      summarise(nb_rare = sum(colSums(across(everything())) > 0), .groups = "drop")
    
    # --- 6. Calcul de l'indice R ---
    indice_rareté <- data.frame(
      triplicat = names(richesse),
      richesse = richesse
    ) %>%
      left_join(nb_rare_df, by = "triplicat") %>%
      mutate(R = nb_rare / richesse,
             site = triplicat)
    
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
      colours = c("darkblue", "coral"),
      name = "LCBD",
      limits = c(min(rarity_lcbd_map_data$LCBD, na.rm = TRUE), max(rarity_lcbd_map_data$LCBD, na.rm = TRUE))
    )
    
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
      coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
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
      coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
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
      coord_sf(xlim = c(63.25, 63.55), ylim = c(-19.85, -19.62), expand = FALSE) +
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
  
  xx <- make_maps_for_group("Ascidiacea", recap_rarity)
  yy <- make_maps_for_group("Porifera", recap_rarity)
  zz <-make_maps_for_group("Bryozoa", recap_rarity)
  
  final_plot <- (tt / xx / yy / zz) +
    plot_layout(guides = "collect")
  
  
  ggsave("outputs/Cartes - LCBD_Rarity/final_map_plot.pdf", plot = final_plot, width = 0.75*18, height = 18)
  
  
  return(NULL)
}

  