# ' plot values of LCBD and rarity on map
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_map_alpha <- function(data_and_meta_clean_fullsites, gps_sites, runa_map, roda_map, runa_reef, roda_reef){
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
  library(dplyr)
  library(iNEXT)
  
  # All taxa map ####
  ## Computing indices ####
  ### richness ####
  
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)

  data_pa <- vegan::decostand(data, "pa")
  
  S <- vegan::specnumber(data_pa, groups = meta$triplicat)
  
  
  triplicat_name = "RODARMS2"
  
  fun_plot_iNEXT_triplicat <- function(triplicat_name, data, meta) {

    
    # Filtrage
    meta_sel <- subset(meta, triplicat == triplicat_name)
    data_sel <- subset(data_pa, meta$triplicat == triplicat_name)
    
    
    # Créer la liste pour iNEXT
    data_list <- list(sel = t(data_sel))
    
    # Appliquer iNEXT
    out <- iNEXT(data_list, q = 0, datatype = "incidence_raw")
    
    # Plot
    ggiNEXT(out, type = 1, color.var = "Order.q") +
      ggtitle(paste("Courbe d'accumulation pour", triplicat_name))
  
    estimates <- out$iNextEst
    
    
    # Filtrer pour 48 unités d'effort
    richness_48 <- lapply(estimates, function(df) {
      df %>%
        filter(t == 48) %>%
        select(Method, Order.q, t, qD, qD.LCL, qD.UCL)
    })
    
    # Affichage
    richness_48 
    
    #Rchesse estimé à 48 faces de plaque pour RODARMS3 : qD = 54.10615
    
    }
  
  # Utilisation
  fun_plot_iNEXT_triplicat("RODARMS3", data, meta)
  
  
  S["RODARMS3"] <- 54
  

  
  #specpool# Mapping the indices ####
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE)
  
  map_data <- data.frame(site = names(S),
                         S = S)
  
  map_data <- as.data.frame(cbind(data_gps, map_data))
  
  # Read the shapefile
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
  
  

  # Modifier les noms de sites dans data_gps
  data_gps$site <- gsub("^RODA", "RODARMS", data_gps$Site)
  data_gps$site <- gsub("^RUNA", "RUNARMS", data_gps$site)
  
  # -- Ajuster les coordonnées pour éviter les chevauchements --
  map_data_runa_sf <- map_data_runa_sf %>%
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
  
  map_data_runa_sf <- map_data_runa_sf %>%
    mutate(
      Longitude = case_when(
        site %in% c("RUNARMS2") ~ Longitude - 0.005,
        site %in% c("RUNARMS3") ~ Longitude + 0.007,
        TRUE ~ Longitude
      )
    )
  
  map_data_runa_sf <- map_data_runa_sf %>%
    mutate(
      label = case_when(
        grepl("P50ARMS", site) ~ gsub("P50ARMS", "P50A", site),
        grepl("RUNARMS", site) ~ gsub("RUNARMS", "RUNA", site),
        grepl("RODARMS", site) ~ gsub("RODARMS", "RODA", site),
        TRUE ~ site
      )
    )
  
  map_data_roda_sf <- map_data_roda_sf %>%
    mutate(
      label = case_when(
        grepl("P50ARMS", site) ~ gsub("P50ARMS", "P50A", site),
        grepl("RUNARMS", site) ~ gsub("RUNARMS", "RUNA", site),
        grepl("RODARMS", site) ~ gsub("RODARMS", "RODA", site),
        TRUE ~ site
      )
    )
  
  map_data_p50a_sf <- map_data_p50a_sf %>%
    mutate(
      label = case_when(
        grepl("P50ARMS", site) ~ gsub("P50ARMS", "P50A", site),
        grepl("RUNARMS", site) ~ gsub("RUNARMS", "RUNA", site),
        grepl("RODARMS", site) ~ gsub("RODARMS", "RODA", site),
        TRUE ~ site
      )
    )
  
  
  # -- Scales communes --
  scale_S <- scale_size_continuous(
    name = "Indice de Richesse spécifique (S)", 
    range = c(2, 10),
    limits = c(min(map_data$S, na.rm = TRUE), max(map_data$S, na.rm = TRUE))
  )

  uu <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_runa_sf,
               aes(x = Longitude, y = Latitude, size = S, fill = "coral"),
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
    scale_S +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)  
  
  vv <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_p50a_sf,
               aes(x = Longitude, y = Latitude, size = S, fill = "coral"),
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
    scale_S +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  ww <- ggplot() +
    geom_sf(data = roda_map, fill = "grey85", color = "black") +
    geom_sf(data = roda_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_roda_sf,
               aes(x = Longitude, y = Latitude, size = S, fill = "coral"),
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
    scale_S +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
  
  
  tt <- (uu | vv | ww) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  
  ggsave("outputs/Cartes - Richness.png", plot = tt, width = 12, height = 7)
  
  
  return(NULL)
}
  