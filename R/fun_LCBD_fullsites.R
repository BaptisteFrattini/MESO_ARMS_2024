#' Beta diversity decomposing
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_LCBD_fullsites <- function(data_and_meta_clean_fullsites){

  library(betapart)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(ggsignif)
  library(ggpubr)
  library(adespatial)
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
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
  
  data_mean_filtered_pa <- vegan::decostand(data_mean_filtered, "pa")
  data_mean_pa <- vegan::decostand(data_mean, "pa")
  
  spe.beta <- adespatial::beta.div(data_mean_filtered_pa, method = "hellinger", nperm = 9999)
  
  spe.beta_bis <- adespatial::beta.div(data_mean_pa, method = "hellinger", nperm = 9999)
  
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
                          LCBD = as.numeric(spe.beta$LCBD))
  
  mean(lcbd_data$LCBD)
  ggplot(lcbd_data, aes(x = reorder(Site, -LCBD), y = LCBD)) +
    geom_point(size = 4) +  # Ajoute des points de taille 4
    scale_color_manual(values = couleurs) +  # Applique les couleurs définies
    geom_hline(yintercept = 0.022, linetype = "dashed", color = "red", size = 1) +  # Ligne moyenne
    coord_flip() +  # Inverse les axes pour une meilleure lisibilité
    labs(title = "LCBD Values by Site", x = "Site", y = "LCBD") +
    theme_minimal()
 return(NULL)   
}
  