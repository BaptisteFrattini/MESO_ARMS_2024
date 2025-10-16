#' Stack column chart 
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to the chart
#' @export
#'
fun_stack_c_chart_fullsites <- function(data_and_meta_clean_fullsites){
  #  data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 

  
  data_pool_mean <- read.csv(data_and_meta_clean_fullsites["path_data_pool_mean"], row.names = 1)
  meta_pool_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  Sediments <- data_pool_mean$Sediments
  Bare_plate <- data_pool_mean$No_Recruitment
  data_pool_mean$Ascidiacea <- rowSums(data.frame(data_pool_mean$ascidiacea_s, data_pool_mean$ascidiacea_c))
  data_pool_mean <- data_pool_mean[, !(colnames(data_pool_mean) %in% c("ascidiacea_c", "ascidiacea_s", "No_Recruitment", "Sediments"))]
  
  
  sums <- colSums(data_pool_mean)
  
  # Ordonner les indices des colonnes en fonction de leur somme
  ordered_indices <- order(sums, decreasing = FALSE)
  
  # Réarranger les colonnes de votre tableau en utilisant l'ordre obtenu
  data_pool_mean <- data_pool_mean[, ordered_indices]
  
  data_pool_mean <- data.frame(cbind(Sediments, Bare_plate, data_pool_mean))
  
  
  spe.T <- colnames(data_pool_mean)
  spe.T <- make.names(spe.T)
  spe.T <- c("Sediment",
             "Bare plate",
             "Cirripedia",
             "Other algae",
             "Hydrozoa",
             "foraminifera",
             "Ascidiacea",
             "Prokariotic biotas",
             "Bivalvia",
             "Porifera",
             "Bryozoa",
             "Crustose coralline algae",
             "Annelida") 
  
  data_pool <- data_pool_mean
  
  data_pool_t <- t(data_pool)
  row.names(data_pool_t) <- spe.T 
  data_pool <- t(data_pool_t)
  
  data_pool <- (data_pool/rowSums(data_pool))*100
  rowSums(data_pool)
  
  pcm <- reshape2::melt(data_pool, id = rownames(as.matrix(data_pool)))
  
  ncol(pcm)
  ####colours####
  col<-NULL
  col[1] <-  "#FAEFD1" 
  col[2] <-  "#899DA4"
  col[3] <-  "navyblue"
  col[4] <-  "#5BBCD6"
  col[5] <-  "darkgreen"
  col[6] <-  "#F1BB7B"
  col[7] <-  "darkorchid1"
  col[8] <-  "#00A08A"
  col[9] <-  "darkolivegreen2"
  col[10] <- "#F2AD00"
  col[11] <- "coral"
  col[12] <- "#FF0000"
  col[13] <- "#446455"
  
  names(col)<- colnames(data_pool)
 
  # pcm_P50A <- subset(pcm, grepl("P50A", Var1))
  # pcm_P50A <- pcm_P50A[order(pcm_P50A$Var1), ]
  # 
  # pcm_RODA <- subset(pcm, grepl("RODA", Var1))
  # pcm_RODA <- pcm_RODA[order(pcm_RODA$Var1), ]
  # 
  # pcm_RUNA <- subset(pcm, grepl("RUNA", Var1))
  # pcm_RUNA <- pcm_RUNA[order(pcm_RUNA$Var1), ]
  # 
  # pcm <- as.data.frame(rbind(pcm_RUNA, pcm_P50A, pcm_RODA))
  
  # Reorder Var1 alphabetically

  
  
  library(ggplot2)
  mx <- ggplot(pcm, aes(x =Var1, fill =Var2, y = value)) + 
    geom_bar(stat = "identity", colour = "black") + 
    theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
          axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
          legend.text = element_text(size = 12, face = "bold", colour = "black"), 
          axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
    scale_y_continuous(expand = c(0,0)) + 
    labs(x = "", y = "Average covers of ARMS plate faces (%)", fill = "Categories :") + 
    scale_fill_manual(values = col)
  mx = mx  + coord_flip()
  
  
  mx = mx + guides(fill = guide_legend(reverse = TRUE))
  
  # library(forcats)
  # library(dplyr)
  # pcm <- pcm %>%
  #   mutate(Var1 = fct_relevel(Var1, sort(unique(as.character(Var1)))))
  
  
  path_to_stack_c_chart <-  here::here("outputs/stack_c_chart_fullsites.pdf")
  ggsave(path_to_stack_c_chart, plot = mx, width = 8, height = 9.5)
  
  return(path_to_stack_c_chart)
}

