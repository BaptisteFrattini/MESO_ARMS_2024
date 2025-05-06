#' Beta diversity decomposing
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_LCBD <- function(data_and_meta_clean){
  
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  library(betapart)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(ggsignif)
  library(ggpubr)
  library(adespatial)
  library(tidyr)
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  msp_list <- names(data_mean) 
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
  
  data_mean_filtered <- data_mean[, msp_list_filter]
  
  # data_mean_filtered_pa <- vegan::decostand(data_mean_filtered, "pa")
  # data_mean_pa <- vegan::decostand(data_mean, "pa")
  
  # spe.beta <- adespatial::beta.div(data_mean_filtered_pa, method = "hellinger", nperm = 9999)
  # ?adespatial::beta.div()
  
  # LCBD ####
  
  matrix.bray <- vegan::vegdist(data_mean_filtered, method = "bray")
  spe.beta <- adespatial::LCBD.comp( matrix.bray, sqrt.D = TRUE)
  
  names(spe.beta$LCBD) <- rownames(data_mean_filtered)
  
  # spe.beta_bis <- adespatial::beta.div(data_mean_pa, method = "hellinger", nperm = 9999)
  
  mean(spe.beta$LCBD)
  
  couleurs <- rep(c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
    "#bcbd22", "#17becf", "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173"
  ), each = 3)
  
  
  couleurs <- couleurs[-18]
  
  mean(spe.beta$LCBD)
  plot( as.factor(meta_mean$arms), spe.beta$LCBD, col=couleurs, pch=20, cex=4, yaxt="n", ylab="", main="Contribution moyenne de chaque site à la diversité Beta - LCBD \n (la ligne rouge représente la LCBD moyenne des différents sites)", xlab="LCBD")                 
  abline(v=0.038, col="red")
  # legend(0.06,20.5, legend=c("1A","1B","1C","2A","2B","2C","3A","3B","3C","4A","4B","4C","5A","5B","5C","6A","6B","6C","7A","7B","7C","8A","8B","8C","9A","9B","9C"), fill=colo, cex=0.5, lty=)
  # ytick=c("RUNARMS1","RUNARMS2","RUNARMS3","RUNARMS4","RUNARMS5","RUNARMS6","RUNARMS7","RUNARMS8","RUNARMS9")
  
  lcbd_data <- data.frame(Site = names(spe.beta$LCBD),
                          LCBD = as.numeric(spe.beta$LCBD),
                          Campain = meta_mean$campain)
  
  
  
  
  ggplot(lcbd_data, aes(x = reorder(Site, -LCBD), y = LCBD, color = Campain)) +
    geom_point(size = 4) +  # Ajoute des points de taille 4
    scale_color_manual(values = c("P50ARMS" = "darkblue", "RODARMS" = "#DD8D29", "RUNARMS" = "#46ACC8")) +  # Applique les couleurs définies
    geom_hline(yintercept = mean(lcbd_data$LCBD), linetype = "dashed", color = "red", size = 1) +  # Ligne moyenne
    coord_flip() +  # Inverse les axes pour une meilleure lisibilité
    labs(title = "LCBD Values by Site", x = "Site", y = "LCBD") +
    theme_minimal()
  
  # SCBD ####
  
  data_mean_filtered_pa <- vegan::decostand(data_mean_filtered, "pa")
  
  spe.beta <- adespatial::beta.div(data_mean_filtered_pa, method = "hellinger", nperm = 9999)
  
  data_SCBD <- data.frame(Species = names(spe.beta$SCBD),
                          SCBD = spe.beta$SCBD)
  
  data_SCBD <- data_SCBD %>%
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
  
  
  
  # Ensure rownames of data_mean_filtered_pa are a column
  data_mean_filtered_pa <- data_mean_filtered_pa %>%
    tibble::rownames_to_column(var = "arms")
  
  # Merge with metadata
  merged_data <- data_mean_filtered_pa %>%
    left_join(meta_mean, by = "arms")
  
  # Summarize species presence by campaign
  species_presence <- merged_data %>%
    select(-arms, -triplicat, -site, -island) %>% 
    group_by(campain) %>%
    summarise(across(everything(), ~ any(. == 1), .names = "{.col}")) %>%
    pivot_longer(-campain, names_to = "species", values_to = "presence")
  
  # Compute presence status across campaigns
  species_status <- species_presence %>%
    pivot_wider(names_from = campain, values_from = presence, values_fill = FALSE) %>%
    mutate(status = case_when(
      P50ARMS & !RUNARMS & !RODARMS ~ "only in P50ARMS",
      !P50ARMS & RUNARMS & !RODARMS ~ "only in RUNARMS",
      !P50ARMS & !RUNARMS & RODARMS ~ "only in RODARMS",
      P50ARMS & !RUNARMS & RODARMS ~ "in both P50ARMS and RODARMS",
      !P50ARMS & RUNARMS & RODARMS ~ "in both RUNARMS and RODARMS",
      P50ARMS & RUNARMS & !RODARMS ~ "in both P50ARMS and RUNARMS",
      P50ARMS & RUNARMS & RODARMS ~ "in all campains"
    ))
  
  # View result
  species_status %>% select(species, status)
  
  
  data_SCBD$status <- species_status$status
  rownames(data_SCBD) <- NULL
  
  
  ggplot(data_SCBD, aes(x = status, y = SCBD, fill = taxa)) +
    geom_bar(stat = "identity", position = "dodge") +  # Bar plot
    labs(x = "Status", y = "SCBD", fill = "Taxa") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels for readability
    scale_fill_brewer(palette = "Set3")  # Nice color scheme
  
  
  
  
  ggplot(data_SCBD, aes(x = taxa, y = SCBD, fill = taxa)) +
    geom_boxplot(alpha = 0.6) +  # Boxplot to show distribution
    geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  # Add points for visibility
    labs(x = "Taxa", y = "SCBD", fill = "Taxa") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
  
  
  # Filter out unwanted taxa
  filtered_data <- data_SCBD %>%
    filter(!taxa %in% c("Cnidaria", "Foraminifera", "Bivalvia"))
  
  # Create the plot
  ggplot(filtered_data, aes(x = taxa, y = SCBD, color = status)) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.8) +  # Jitter plot with status-based colors
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +  # Mean as black diamond
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +  # SD as error bars
    labs(x = "Taxa", y = "SCBD", color = "Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_brewer(palette = "Set3")  # Use a nice color scheme for status
  
  
  # Compute occurrence frequency
  species_occurence <- colSums(data_mean_filtered_pa[,-1])/nrow(data_mean_filtered_pa)
  
  species_occurence <- species_occurence[species_occurence != 0]
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  
  
  # Merge with data_SCBD to keep SCBD values
  data_SCBD <- left_join(data_SCBD, species_occurence, by = "Species")
  
  filtered_data <- data_SCBD %>%
    filter(!taxa %in% c("Cnidaria", "Foraminifera", "Bivalvia"))
  
  # Plot the results
  ggplot(filtered_data, aes(x = taxa, y = SCBD, color = f)) +
    geom_point(size = 3, alpha = 0.8) +  # Scatter plot with color-coded points
    scale_color_gradient(low = "blue", high = "red") +  # Color scale for occurrence frequency
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "darkgrey") +  # SD as error bars
    theme_minimal()
  
  
  
  return(NULL)   
}