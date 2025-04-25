# ' 
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_meta_explo <- function(data_and_meta_clean_fullsites){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites")

  library(vegan)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  
  
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)
  
  meta_roda <- subset(meta, meta$campain == "RODARMS" )
  meta_p50a <- subset(meta, meta$campain == "P50ARMS" )
  meta <- as.data.frame(rbind(meta_p50a, meta_roda))

  
  # Créer une variable combinée pour l'ouverture et l'orientation
  meta <- meta %>%
    mutate(orientation_type = paste0(orientation, "_", open_close))
  
  # Compter le nombre de faces par île, site, orientation et ouverture
  face_counts <- meta %>%
    group_by(island, site_names, orientation, orientation_type, open_close) %>%
    summarise(n_faces = n(), .groups = "drop")
  
  face_counts_p50a <- subset(face_counts, island == "Reunion" )
  
  
  # Graphique : nombre de faces par orientation et ouverture, par île et site
  aa <- ggplot(face_counts_p50a, aes(x = interaction(orientation, open_close), y = n_faces, fill = orientation_type)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_grid(island ~ site_names, scales = "free_x", space = "free_x") +
        labs(
          x = "Orientation et Ouverture",
          y = "Nombre de faces",
          fill = "Ouverture",
          title = "Déséquilibre d'échantillonnage des ARMS par orientation et ouverture"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 20)
  
  
  face_counts_roda <- subset(face_counts, island == "Rodrigues" )
  
  
  # Graphique : nombre de faces par orientation et ouverture, par île et site
  bb <- ggplot(face_counts_roda, aes(x = interaction(orientation, open_close), y = n_faces, fill = orientation_type)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_grid(island ~ site_names, scales = "free_x", space = "free_x") +
        labs(
          x = "Orientation et Ouverture",
          y = "Nombre de faces",
          fill = "Ouverture",
          title = "Déséquilibre d'échantillonnage des ARMS par orientation et ouverture"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 20)
      
      
  # Identifier les lignes B_o et T_o à Rodrigues
  rodrigues_to_filter <- meta %>%
    filter(island == "Rodrigues",
           (orientation == "B" & open_close == "o") |
             (orientation == "T" & open_close == "o"))
  
  
  # Pour chaque ARMS, retirer 2 B_o et 1 T_o avec les plus petits plate_number
  rows_to_remove <- rodrigues_to_filter %>%
    group_by(arms) %>%
    arrange(plate_number) %>%
    group_modify(~ {
      b_o <- .x %>% filter(orientation == "B", open_close == "o") %>% slice_head(n = 1)
    }) %>%
    ungroup()
  
  
  
  # Retirer les lignes sélectionnées du jeu de données original
  meta_filtered <- anti_join(meta, rows_to_remove)
  
  # Vérifier les lignes supprimées
  rows_to_remove %>%
    select(name, arms, orientation, open_close, plate_number, orientation_type)
  
  face_counts <- meta_filtered %>%
    group_by(island, site_names, orientation, open_close, orientation_type) %>%
    summarise(n_faces = n(), .groups = "drop")
  
  face_counts_p50a <- subset(face_counts, island == "Reunion" )
  
  cc <- ggplot(face_counts_p50a, aes(x = interaction(orientation, open_close), y = n_faces, fill = orientation_type)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_grid(island ~ site_names, scales = "free_x", space = "free_x") +
        labs(
          x = "Orientation et Ouverture",
          y = "Nombre de faces",
          fill = "Ouverture",
          title = "Déséquilibre des P50ARMS par orientation et ouverture"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 20)
  
  
  face_counts_roda <- subset(face_counts, island == "Rodrigues" )

  dd <- ggplot(face_counts_roda, aes(x = interaction(orientation, open_close), y = n_faces, fill = orientation_type)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_grid(island ~ site_names, scales = "free_x", space = "free_x") +
        labs(
          x = "Orientation et Ouverture",
          y = "Nombre de faces",
          fill = "Ouverture",
          title = "Déséquilibre d'échantillonnage des ARMS par orientation et ouverture"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 20)
  
  face_totals <- meta_filtered %>%
    group_by(campain) %>%
    summarise(total_faces = n(), .groups = "drop")
  
  # Afficher le tableau
  print(face_totals)
  
  
  rows_to_remove <- rodrigues_to_filter %>%
    group_by(arms) %>%
    arrange(plate_number) %>%
    group_modify(~ {
      b_o <- .x %>% filter(orientation == "B", open_close == "o") %>% slice_head(n = 2)
      t_o <- .x %>% filter(orientation == "T", open_close == "o") %>% slice_head(n = 1)
      bind_rows(b_o, t_o)
    }) %>%
    ungroup()
  
  meta_filtered_2 <- anti_join(meta, rows_to_remove)
  
  # Vérifier les lignes supprimées
  rows_to_remove %>%
    select(name, arms, orientation, open_close, plate_number, orientation_type)
  
  face_counts <- meta_filtered_2 %>%
    group_by(island, site_names, orientation, open_close, orientation_type) %>%
    summarise(n_faces = n(), .groups = "drop")
  
  face_counts_p50a <- subset(face_counts, island == "Reunion" )
  
  ee <- ggplot(face_counts_p50a, aes(x = interaction(orientation, open_close), y = n_faces, fill = orientation_type)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_grid(island ~ site_names, scales = "free_x", space = "free_x") +
        labs(
          x = "Orientation et Ouverture",
          y = "Nombre de faces",
          fill = "Ouverture",
          title = "Déséquilibre des P50ARMS par orientation et ouverture"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 20)
  
  
  face_counts_roda <- subset(face_counts, island == "Rodrigues" )
  
  ff <- ggplot(face_counts_roda, aes(x = interaction(orientation, open_close), y = n_faces, fill = orientation_type)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        facet_grid(island ~ site_names, scales = "free_x", space = "free_x") +
        labs(
          x = "Orientation et Ouverture",
          y = "Nombre de faces",
          fill = "Ouverture",
          title = "Déséquilibre d'échantillonnage des ARMS par orientation et ouverture"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 20)
  
  face_totals <- meta_filtered_2 %>%
    group_by(campain) %>%
    summarise(total_faces = n(), .groups = "drop")
  
  # Afficher le tableau
  print(face_totals)
  
  rows_to_remove <- rodrigues_to_filter %>%
    filter(triplicat %in% c("RODARMS1", "RODARMS2")) %>%  # <-- Ajout ici
    group_by(arms) %>%
    arrange(plate_number) %>%
    group_modify(~ {
      b_o <- .x %>% filter(orientation == "B", open_close == "o") %>% slice_head(n = 2)
      t_o <- .x %>% filter(orientation == "T", open_close == "o") %>% slice_head(n = 1)
      bind_rows(b_o, t_o)
    }) %>%
    ungroup()
  
  meta_filtered_3 <- anti_join(meta, rows_to_remove)

  
  rows_to_remove %>%
    select(name, arms, orientation, open_close, plate_number, orientation_type)
  
  face_counts <- meta_filtered_3 %>%
    group_by(island, site_names, orientation, open_close, orientation_type) %>%
    summarise(n_faces = n(), .groups = "drop")
  
  face_counts_p50a <- subset(face_counts, island == "Reunion" )
  
  gg <- ggplot(face_counts_p50a, aes(x = interaction(orientation, open_close), y = n_faces, fill = orientation_type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(island ~ site_names, scales = "free_x", space = "free_x") +
    labs(
      x = "Orientation et Ouverture",
      y = "Nombre de faces",
      fill = "Ouverture",
      title = "Déséquilibre des P50ARMS par orientation et ouverture"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 20)
  
  
  face_counts_roda <- subset(face_counts, island == "Rodrigues" )
  
  hh <- ggplot(face_counts_roda, aes(x = interaction(orientation, open_close), y = n_faces, fill = orientation_type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(island ~ site_names, scales = "free_x", space = "free_x") +
    labs(
      x = "Orientation et Ouverture",
      y = "Nombre de faces",
      fill = "Ouverture",
      title = "Déséquilibre d'échantillonnage des ARMS par orientation et ouverture"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 20)
  
  face_totals <- meta_filtered_3 %>%
    group_by(campain) %>%
    summarise(total_faces = n(), .groups = "drop")
  
  # Afficher le tableau
  print(face_totals)
  
  return(NULL)
  
  }

  
    