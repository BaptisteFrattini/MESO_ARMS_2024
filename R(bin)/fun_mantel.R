#' Mantel test
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_mantel <- function(data_and_meta_clean, gps_sites){
  # data_and_meta_clean = targets::tar_read("clean_data_metadata")
  # gps_sites = targets::tar_read("data_gps_sites")
  library(betapart)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE)
  data_gps <-  data_gps %>%
    slice(rep(1:n(), each = 3))
  data_gps <- data_gps[-18,]
  
  data_gps
  
  coords <- data_gps[, c("Longitude", "Latitude")]

  matrix.dist <- geosphere::distm(coords)
  matrix.dist <- as.data.frame(matrix.dist)
  matrix.dist <- as.matrix(matrix.dist)
  matrix.dist <- as.dist(matrix.dist)
  class(matrix.dist)
  # row.names(matrix.dist) <- data_gps$Site
  # colnames(matrix.dist) <- data_gps$Site
  
  ## representation of the set of distances 
  a <- as.numeric(levels(factor(as.vector(as.matrix((matrix.dist))))))
  b <- c(1:length(a))
  plot(x = a, y = b, yaxt="n", frame = FALSE, xlab = "Distance modalities between sites (km)", ylab = " ", pch = 19, col = "black", cex = 1.2, cex.axis = 2, cex.lab = 1.65)
  dev.off()
  
  mat.bc <- vegan::vegdist(data_mean, dist = "bray")
  data_mean_pa <- vegan::decostand(data_mean, method = "pa")
  B.pair.pa <- betapart::beta.pair(data_mean_pa, index.family = "jaccard")
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- B.pair.pa$beta.jac

  ##### jacc #####

  # aa = as.vector(mat.bc)
  # tt = as.vector(matrix.dist)
  # #new data frame with vectorized distance matrices
  # mat = data.frame(aa,tt)
  # 
  # mm1 = ggplot2::ggplot(mat, aes(y = aa, x = tt)) + 
  #   geom_point(size = 3, alpha = 0.5, color = "black") + 
  #   labs(x = NULL,
  #        y = "Jaccard dissimilarity") +
  #   geom_smooth(method = "gam", 
  #               colour = "red", 
  #               alpha = 0.2, 
  #               fill = "red") +
  #   theme( axis.text.x = element_text(colour = "black",
  #                                     size = 12), 
  #          axis.text.y = element_text(size = 11, 
  #                                     colour = "black"), 
  #          axis.title = element_text(size = 14, 
  #                                    colour = "black"), 
  #          panel.background = element_blank(), 
  #          panel.border = element_rect(fill = NA,
  #                                      colour = "black")) +
  #   scale_x_continuous(breaks = seq(0, 100, 10), minor_breaks = seq(0, 100, 1)) +
  #   ylim(0,0.7)
  
  
  return(NULL)
  
}

  
  
  