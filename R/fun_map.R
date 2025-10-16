# ' plot values of LCBD and rarity on map
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_map <- function(data_and_meta_clean_fullsites, gps_sites, runa_map, roda_map, runa_reef, roda_reef, world_map){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites")
  # gps_sites = targets::tar_read("data_gps_sites")
  # runa_map = targets::tar_read("map_runa")
  # roda_map = targets::tar_read("map_roda")
  # runa_reef = targets::tar_read("reef_runa")
  # roda_reef = targets::tar_read("reef_roda")
  # world_map = targets::tar_read("map_world")
  
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
  
  # map ####
  
  #specpool# Mapping the indices ####
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE)
  
  map_data <- as.data.frame(data_gps)
  
  # Read the shapefile
  runa_map <- st_read(runa_map)
  roda_map <- st_read(roda_map)
  world_map <- st_read(world_map)
  
  
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
  
  
  
  # importer les reef et la world map
  runa_reef <- st_read(runa_reef)
  roda_reef <- st_read(roda_reef)
  
  runa_reef <- st_transform(runa_reef, crs = 4326)
  roda_reef <- st_transform(roda_reef, crs = 4326)
  world_map <- st_transform(world_map, crs = 4326)
  
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
        Site %in% c("RUNA2") ~ Latitude + 0.002,
        Site %in% c("RUNA3") ~ Latitude - 0.002,
        Site %in% c("RUNA4") ~ Latitude + 0.001,
        Site %in% c("RUNA5") ~ Latitude - 0.001,
        Site %in% c("RUNA6") ~ Latitude - 0.001,
        TRUE ~ Latitude
      )
    )
  
  # -- Ajuster les coordonnées pour éviter les chevauchements --
  map_data_cina_sf <- subset(map_data_runa_sf, map_data_runa_sf$Site == "RUNA2")
  
  map_data_cina_sf <- map_data_cina_sf %>%
    mutate(
      Site = case_when(
        grepl("RUNA2", Site) ~ gsub("RUNA2", "CINA1-4", Site),
        TRUE ~ Site
      )
    )
    
  

  
  # Define common width and height
  map_width <- 0.45
  map_height <- 0.4
  
  # Center coordinates for each map
  center_p50a <- c(x = mean(c(55.15, 55.6)), y = mean(c(-21.4, -21)))
  center_roda <- c(x = mean(c(63.20, 63.58)), y = mean(c(-19.90, -19.58)))
  center_world <- c(x = mean(c(42.42, 67.4)), y = mean(c(-30.1, -10.17)))
  # (Add RUNA similarly if needed)
  
  # Create coordinate limits
  p50a_xlim <- c(center_p50a['x'] - map_width / 2, center_p50a['x'] + map_width / 2)
  p50a_ylim <- c(center_p50a['y'] - map_height / 2, center_p50a['y'] + map_height / 2)
  
  roda_xlim <- c(center_roda['x'] - map_width / 2, center_roda['x'] + map_width / 2)
  roda_ylim <- c(center_roda['y'] - map_height / 2, center_roda['y'] + map_height / 2)
  
  world_xlim <- c(42.42, 67.4)
  world_ylim <- c(-27,-10.17)
  
  # Then apply in each coord_sf:
  coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE)
  coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE)
  coord_sf(xlim = world_xlim, ylim = world_ylim, expand = FALSE)
  
  ## just a map for ppt ####
  
  masca <- ggplot() +
    geom_sf(data = world_map, fill = "grey85", color = "black") +
    coord_sf(xlim = world_xlim, ylim = world_ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = "Mascarene archipelago") +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)  
  
  ggsave("outputs/Cartes_des_masca.png", plot = masca, width = 12/2, height = 7/2)
  
  xx <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_runa_sf,
               aes(x = Longitude, y = Latitude),
               fill = "black",
               shape = 20, stroke = 0.6) +
    geom_text_repel(
      data = map_data_runa_sf,
      aes(x = Longitude, y = Latitude, label = Site),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    theme_minimal() +
    labs(title = "RUNA") +
    annotation_scale(location = "bl", width_hint = 0.25) 
    # annotation_north_arrow(location = "tr", which_north = "true",
    #                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
    #                        style = north_arrow_fancy_orienteering)  
  
  xx_bis <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_cina_sf,
               aes(x = Longitude, y = Latitude),
               fill = "black",
               shape = 20, stroke = 0.6) +
    geom_text_repel(
      data = map_data_cina_sf,
      aes(x = Longitude, y = Latitude, label = Site),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    theme_minimal() +
    labs(title = "CINA") +
    annotation_scale(location = "bl", width_hint = 0.25) 
    # + annotation_north_arrow(location = "tr", which_north = "true",
                           # pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                           # style = north_arrow_fancy_orienteering) 
  
  map_data_p50a_sf$Site <- c("P50A1", "P50A2", "P50A3")
  
  yy <- ggplot() +
    geom_sf(data = runa_map, fill = "grey85", color = "black") +
    geom_sf(data = runa_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_p50a_sf,
               aes(x = Longitude, y = Latitude), 
               fill = "black",
               shape = 20, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_p50a_sf,
      aes(x = Longitude, y = Latitude, label = Site),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    coord_sf(xlim = p50a_xlim, ylim = p50a_ylim, expand = FALSE) +
    # coord_sf(xlim = c(55.15, 55.6), ylim = c(-21.4, -21), expand = FALSE) +
    theme_minimal() +
    labs(title = "P50A") +
    annotation_scale(location = "bl", width_hint = 0.25) 
    # annotation_north_arrow(location = "tr", which_north = "true",
    #                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
    #                        style = north_arrow_fancy_orienteering)
  
  zz <- ggplot() +
    geom_sf(data = roda_map, fill = "grey85", color = "black") +
    geom_sf(data = roda_reef, fill = "darkcyan", color = "darkcyan", alpha = 0.5) +
    geom_point(data = map_data_roda_sf,
               aes(x = Longitude, y = Latitude),
               fill = "black",
               shape = 20, color = "black", stroke = 0.6) +
    geom_text_repel(
      data = map_data_roda_sf,
      aes(x = Longitude, y = Latitude, label = Site),
      size = 4,
      fontface = "bold",
      color = "black",
      box.padding = 0.8,         # Distance entre point et étiquette (en "lines")
      max.overlaps = Inf,        # Ne jamais supprimer d’étiquettes
      force = 1.5,               # Force de répulsion
      segment.size = 0.2,        # Taille du segment
      segment.color = "grey30"   # Couleur du trait
    ) +
    coord_sf(xlim = roda_xlim, ylim = roda_ylim, expand = FALSE) +
    # coord_sf(xlim = c(63.20, 63.58), ylim = c(-19.90, -19.58), expand = FALSE) +
    theme_minimal() +
    labs(title = "RODA") +
    annotation_scale(location = "bl", width_hint = 0.25) 
    # annotation_north_arrow(location = "tr", which_north = "true",
    #                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
    #                        style = north_arrow_fancy_orienteering)
    # 
  aa <- (xx | xx_bis | yy | zz) +
    plot_layout(guides = "collect") 
  
  ggsave("outputs/Cartes_des_sites.png", plot = aa, width = 15, height = 9)
  
  return(NULL)
}
