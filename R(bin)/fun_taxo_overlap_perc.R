#' Aoccurence classification of MSP
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_taxo_overlap_perc <- function(data_and_meta_clean){
  
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  data <- read.csv(data_and_meta_clean["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean["path_meta"], row.names = 1)
  data <- subset(data, meta$island == "Reunion")
  meta <- subset(meta, meta$island == "Reunion")
  
  # Load necessary libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  
  # Step 3: Filter the rows in `data_mean` corresponding to each campaign
  p50arms_data <- subset(data, meta$campain == "P50ARMS")
  runarms_data <- subset(data, meta$campain == "RUNARMS")
  
  # Step 4: Compute the cumulative abundances for each species
  # p50arms_abundance <- colSums(p50arms_data)
  # runarms_abundance <- colSums(runarms_data)
  
  p50arms_presence <- vegan::decostand(p50arms_data, "pa")
  runarms_presence <- vegan::decostand(runarms_data, "pa")
  
  p50arms_perc <- (colSums(p50arms_presence)*100)/length(rownames(p50arms_presence))
  runarms_perc <- (colSums(runarms_presence)*100)/length(rownames(runarms_presence))
  
  # Step 5: Combine into a single data frame
  species_abundance <- data.frame(
    Species = colnames(data),
    P50ARMS = p50arms_perc,
    RUNARMS = runarms_perc,
    row.names = NULL
  )
  species_abundance <- species_abundance[!species_abundance$Species=="No_Recruitment",]
  species_abundance <- species_abundance[!species_abundance$Species=="Sediments",]
  
  # species_abundance$ab_total <- species_abundance$P50ARMS + species_abundance$RUNARMS
  # 
  # species_abundance <- subset(species_abundance, ab_total != 0)
  # 
  # species_abundance_run_sorted <- species_abundance %>%
  #   arrange(RUNARMS)
  # species_abundance_run_sorted <- subset(species_abundance_run_sorted, RUNARMS != 0)
  # 
  # thresh_run <- species_abundance_run_sorted$ab_total[0.25*length(species_abundance_run_sorted$Species)]
  # 
  # species_abundance_p50_sorted <- species_abundance %>%
  #   arrange(P50ARMS)
  # species_abundance_p50_sorted <- subset(species_abundance_p50_sorted, P50ARMS != 0)
  # 
  # thresh_p50 <- species_abundance_p50_sorted$ab_total[0.25*length(species_abundance_p50_sorted$Species)]

  
  # Ajouter la colonne "Status" avec la nouvelle condition
  
  species_abundance$Status <- ifelse(
    species_abundance$P50ARMS == 0 & species_abundance$RUNARMS == 0, 
    "absent", 
    ifelse(
      species_abundance$P50ARMS > 25 & species_abundance$RUNARMS == 0, 
      "Common deep exclusive", 
      ifelse(
        species_abundance$RUNARMS > 25 & species_abundance$P50ARMS == 0, 
        "Common shallow exclusive", 
        ifelse(
          species_abundance$P50ARMS > 25 & species_abundance$RUNARMS > 25, 
          "Common generalist", 
          ifelse(
            species_abundance$P50ARMS < 25 & species_abundance$RUNARMS == 0, 
            "Rare deep exclusive", 
            ifelse(
              species_abundance$RUNARMS < 25 & species_abundance$P50ARMS == 0, 
              "Rare shallow exclusive", 
              ifelse(
                species_abundance$P50ARMS < 25 & species_abundance$RUNARMS < 25, 
                "Rare generalist", 
                ifelse(
                  species_abundance$P50ARMS > 25 & species_abundance$RUNARMS < 25, 
                  "Deep specialist",
                  ifelse(
                    species_abundance$P50ARMS < 25 & species_abundance$RUNARMS > 25, 
                    "Shallow specialist",
                    NA  # Si aucune condition n'est remplie
                  )
                )           
              )
            )
          )
        )
      )
    )
  )
  #delete species that are absent
  species_abundance <- subset(species_abundance, Status != "absent")

  # species_abundance <- species_abundance[!species_abundance$Species=="X_SPON",]
  # Visualiser le tableau mis à jour
  print(species_abundance)
  
  
  library(ggplot2)
  
  species_long <- species_abundance %>%
    pivot_longer(
      cols = c(P50ARMS, RUNARMS),  # Columns for each campaign
      names_to = "Campaign",
      values_to = "Percentage"
    )
  
  
  # species_long <- species_long %>%
  #   mutate(Percentage = log1p(Percentage))
  
  species_long <- species_long %>%
    mutate(Status = factor(Status, levels = c(
      "Common deep exclusive",
      "Common shallow exclusive",
      "Common generalist",
      "Rare deep exclusive",
      "Rare shallow exclusive",
      "Rare generalist",
      "Deep specialist",
      "Shallow specialist"
    )))
  
  g1 <- ggplot(species_long, aes(x = reorder(Species, -Percentage), y = Percentage, fill = Status)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Campaign, scales = "fixed") +  # Uniforme entre les campagnes
    coord_flip() +
    scale_fill_manual(values = c(
      "Common deep exclusive" = "darkblue",
      "Common shallow exclusive" = "darkgreen",
      "Common generalist" = "grey41",
      "Deep specialist" = "slateblue",
      "Shallow specialist" = "mediumseagreen",
      "Rare deep exclusive" = "skyblue",
      "Rare shallow exclusive" = "olivedrab1",
      "Rare generalist" = "grey"
    ), drop = FALSE) +
    labs(
      title = "Species Percent Presence by depth",
      x = "Species",
      y = "Morpho-species occurence percentage"
    ) +
    theme_minimal()
  

  
  ggsave("outputs/taxo_overlap_barplot_perc.pdf", g1, width = 9, height = 15 )
  

  # Ajouter la colonne "taxa"
  species_abundance <- species_abundance %>%
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
  
  
  
  status_colors <- c(
    "Common deep exclusive" = "darkblue",
    "Common shallow exclusive" = "darkgreen",
    "Common generalist" = "grey41",
    "Deep specialist" = "slateblue",
    "Shallow specialist" = "mediumseagreen",
    "Rare deep exclusive" = "skyblue",
    "Rare shallow exclusive" = "olivedrab1",
    "Rare generalist" = "grey"
  )
  
  
  histogram_data <- species_abundance %>%
    group_by(taxa, Status) %>%
    summarise(count = n(), .groups = "drop")  
  
  
  # Étape 1 : Ajouter les catégories absentes à Status
  # Étape 1 : Définir les niveaux souhaités pour Status
  levels_status <- c(
    "Common deep exclusive",
    "Common shallow exclusive",
    "Common generalist",
    "Rare deep exclusive",
    "Rare shallow exclusive",
    "Rare generalist",
    "Deep specialist",
    "Shallow specialist"
  )
  
  # Étape 2 : Convertir Status en facteur avec les niveaux définis
  histogram_data$Status <- factor(histogram_data$Status, levels = levels_status)
  
  # Étape 3 : Ajouter des lignes fictives pour les catégories manquantes
  for (status in levels_status) {
    if (!status %in% histogram_data$Status) {
      histogram_data <- rbind(histogram_data, data.frame(
        taxa = NA,  # Pas de taxa pour les lignes fictives
        count = 0,  # Pas de contribution réelle
        Status = status
      ))
    }
  }
  
  # Étape 4 : Créer le graphique avec les couleurs personnalisées
  cumulative_histogram <- ggplot(histogram_data, aes(x = taxa, y = count, fill = Status)) +
    geom_bar(stat = "identity", position = "stack", color = "white") +
    scale_fill_manual(
      values = c(
        "Common deep exclusive" = "darkblue",
        "Common shallow exclusive" = "darkgreen",
        "Common generalist" = "grey41",
        "Deep specialist" = "slateblue",
        "Shallow specialist" = "mediumseagreen",
        "Rare deep exclusive" = "skyblue",
        "Rare shallow exclusive" = "olivedrab1",
        "Rare generalist" = "grey"
      )
    ) +
    labs(
      title = "Cumulative Histogram of Species Status by Taxa",
      x = "Taxa",
      y = "Number of Species",
      fill = "Status"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Afficher le graphique
  print(cumulative_histogram)
  
  ggsave("outputs/Cumul_hist_perc.pdf", cumulative_histogram, width = 8, height = 5 )
  
  
  return(NULL)
}
