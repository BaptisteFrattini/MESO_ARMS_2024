#' occ freq of species in plates arms and sites
#'
#' @param data_and_meta_clean_fullsites the path to the clean data
#' 
#'
#' @return NULL
#' @export
#'
fun_f_occ_taxa <- function(data_and_meta_clean_fullsites){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(mgcv)
  library(cowplot)
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)
  # data <- subset(data, meta$island == "Reunion")
  # meta <- subset(meta, meta$island == "Reunion")
  # 
  
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
  
  
  # fréquence d'occurence à l'echele des faces de plaques ####
  f <- colSums(data_pa)/nrow(data_pa)
  
  f <- data.frame(f = f,
                  Species = names(f))
  
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
  

  f_taxa <- left_join(corr_taxa, f, by = "Species")
  
  f_taxa$log_f <- log(f_taxa$f)
    
  f_taxa_filtered <- f_taxa %>%
    filter(taxa %in% c("Ascidiacea", "Bryozoa", "Porifera"))
  
  vv <- ggplot(f_taxa_filtered, aes(x = taxa, y = log_f)) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  # Boîte vide et noire, sans outliers
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +  # Points légèrement décalés horizontalement
    theme_minimal() +
    labs(
      x = "Taxa",
      y = "Log(occurrence frequency)",
      title = ""
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16),  # Titre axe X
      axis.title.y = element_text(size = 16),  # Titre axe Y
      axis.text.x  = element_text(size = 16),  # Valeurs axe X
      axis.text.y  = element_text(size = 16)
    )
  
  
  tt <- ggplot(f_taxa_filtered, aes(x = taxa, y = f)) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  # Boîte vide et noire, sans outliers
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +  # Points légèrement décalés horizontalement
    theme_minimal() +
    labs(
      x = "Taxa",
      y = "Occurrence frequency",
      title = ""
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 16),  # Titre axe X
      axis.title.y = element_text(size = 16),  # Titre axe Y
      axis.text.x  = element_text(size = 16),  # Valeurs axe X
      axis.text.y  = element_text(size = 16)
    )
  
  fig1 <- cowplot::plot_grid(tt, vv, ncol = 2, nrow = 1)
  
  path_to_stack_c_chart <-  here::here("outputs/plot_f_occ_taxa.pdf")
  ggsave(path_to_stack_c_chart, plot = fig1, width =12, height = 9)
  
  return(NULL)
  
  
}