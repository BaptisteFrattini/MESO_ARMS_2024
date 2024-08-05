#' Explore the cleaned data
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to the subseted raw data file
#' @export
#'
fun_data_exploring <- function(data_and_meta_clean){
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  
  data_and_meta_clean
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  
  #### alpha diversity ####
  library(ggpubr)
  library(forcats)
  library(ggplot2)
  S <- vegan::specnumber(data_mean)
  
  data_div <- data.frame(s = S,
                         name = meta_mean$arms,
                         campain = meta_mean$campain,
                         triplicat = meta_mean$triplicat)

  ff <- ggplot(data_div, aes(x = fct_relevel(campain, "RUNARMS", "RODARMS", "P50ARMS"), y = S)) +
    geom_boxplot(fill =  c("lightblue", "lightblue3","darkblue")) +
    labs(title = "",
         x = "Depth",
         y = "Average species richness in an ARMS") +
    theme_classic()
  
  gg <- ggplot(data_div, aes(x = fct_relevel(triplicat, "RUNARMS1","RUNARMS5","RUNARMS9", "RODARMS1", "RODARMS2","RODARMS3", "P50ARMS1", "P50ARMS2", "P50ARMS3"), y = S)) +
    geom_boxplot(fill =  c(rep("lightblue",3),rep("lightblue3",3), rep("darkblue",3))) +
    labs(title = "",
         x = "Depth",
         y = "Average species richness in an ARMS") +
    theme_classic()
 
  ggsave("outputs/alpha_div(1).pdf", plot = ff, width = 9, height = 8)
  ggsave("outputs/alpha_div(2).pdf", plot = gg, width = 9, height = 8)
  
  ####NMDS####
  #bray
  library(wesanderson)
  fit <-  vegan::metaMDS(data_mean, distance = "bray")
  col1 <- c(wesanderson::wes_palette("FantasticFox1", 5), wesanderson::wes_palette("Royal2", 1))
  col1 <- c("#DD8D29", "#E2D200", "#46ACC8", "#B40F20","darkblue", "#9A8822")
  mds_name <- paste0("NMDS_MESO_bay.pdf")
  mds_path <- here::here("outputs/", mds_name)
  pdf(file =  mds_path)
  plot(fit)
  vegan::ordihull(fit, meta_mean$triplicat, col=col1)
  vegan::ordiellipse(fit, meta_mean$triplicat, col=col1, kind = "ehull", lwd=2)
  vegan::ordispider(fit, meta_mean$triplicat, col=col1, cex = 0.5)
  vegan::ordiellipse(fit, meta_mean$triplicat, col=col1, draw="polygon", label = TRUE)
  # vegan::ordilabel(fit, display = "species", choices = c(1, 2), cex = 0.3, border = NA)
  a <- paste0("stress = ",round(fit$stress, 3))
  text(0.8,0.95, a)
  dev.off()
  #Jac
  data_mean_pa <- vegan::decostand(data_mean, "pa")
  fit2 <-  vegan::metaMDS(data_mean_pa, distance = "jaccard")
  ?metaMDS
  mds_name <- paste0("NMDS_MESO_jac.pdf")
  mds_path <- here::here("outputs/", mds_name)
  pdf(file =  mds_path)
  plot(fit2)
  vegan::ordihull(fit2, meta_mean$triplicat, col=col1)
  vegan::ordiellipse(fit2, meta_mean$triplicat, col=col1, kind = "ehull", lwd=2)
  vegan::ordispider(fit2, meta_mean$triplicat, col=col1, cex = 0.5)
  vegan::ordiellipse(fit2, meta_mean$triplicat, col=col1, draw="polygon", label = TRUE)
  # vegan::ordilabel(fit2, display = "species", choices = c(1, 2), cex = 0.3, border = NA)
  a <- paste0("stress = ",round(fit2$stress, 3))
  text(0.8,0.8, a)
  dev.off()
  
  #### species accumulation curves #### 
  
  data <- read.csv(data_and_meta_clean["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean["path_meta"], row.names = 1)
  
  #### All MESO ARMS combined ####
  
  div_alpha_name <- paste0("div_alpha_plate.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 15, height = 5)
  
  par(mfrow = c(1, 3))
  
  data_P50 <- subset(data, meta$campain == "P50ARMS")
   
  s <- vegan::specaccum(data_P50, method = "random", permutations = 999,
                        conditioned =TRUE)
  pool <- vegan::specpool(data_P50)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in all 9 mesophotic ARMS",
       ylim=c(1,100))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(80,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  #### RODARMS ####
  
  data_RODARMS <- subset(data, meta$campain == "RODARMS")
  
  s <- vegan::specaccum(data_RODARMS, method = "random", permutations = 999,
                        conditioned =TRUE)
  
  pool <- vegan::specpool(data_RODARMS)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in all 8 ARMS from Rodrigues",
       ylim=c(1,100))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(100,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  #### RUNARMS ####
  
  data_RUNARMS <- subset(data, meta$campain == "RUNARMS")
  
  s <- vegan::specaccum(data_RUNARMS, method = "random", permutations = 999,
                        conditioned =TRUE)
  
  pool <- vegan::specpool(data_RUNARMS)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in all 8 ARMS from Reunion",
       ylim=c(1,100))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(80,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()
  
  # #### P50ARMS1 ####
  # 
  # data_P50ARMS1 <- subset(data, meta$triplicat == "P50ARMS1")
  # 
  # s <- vegan::specaccum(data_P50ARMS1, method = "random", permutations = 999,
  #                       conditioned =TRUE)
  # 
  # pool <- vegan::specpool(data_P50ARMS1)
  # 
  # plot(s, 
  #      ci.type="poly",
  #      col="blue", 
  #      lwd=2,
  #      ci.lty=0,
  #      ci.col="lightblue",
  #      xlab="Number of plates analysed",
  #      ylab="# morpho-species detected in all 3 mesophotic ARMS of P50ARMS1",
  #      ylim=c(1,100))
  # 
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  # 
  # 
  # text(26,
  #      10,
  #      paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
  #      cex = 0.85)
  # 
  # 
  # 
  # #### P50ARMS1 ####
  # 
  # data_P50ARMS2 <- subset(data, meta$triplicat == "P50ARMS2")
  # 
  # s <- vegan::specaccum(data_P50ARMS2, method = "random", permutations = 999,
  #                       conditioned =TRUE)
  # 
  # pool <- vegan::specpool(data_P50ARMS2)
  # 
  # plot(s, 
  #      ci.type="poly",
  #      col="blue", 
  #      lwd=2,
  #      ci.lty=0,
  #      ci.col="lightblue",
  #      xlab="Number of plates analysed",
  #      ylab="# morpho-species detected in all 3 mesophotic ARMS of P50ARMS2",
  #      ylim=c(1,100))
  # 
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  # 
  # 
  # text(25,
  #      10,
  #      paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
  #      cex = 0.85)
  # 
  # #### P50ARMS3 ####
  # 
  # data_P50ARMS3 <- subset(data, meta$triplicat == "P50ARMS3")
  # 
  # s <- vegan::specaccum(data_P50ARMS3, method = "random", permutations = 999,
  #                       conditioned =TRUE)
  # 
  # pool <- vegan::specpool(data_P50ARMS3)
  # 
  # plot(s, 
  #      ci.type="poly",
  #      col="blue", 
  #      lwd=2,
  #      ci.lty=0,
  #      ci.col="lightblue",
  #      xlab="Number of plates analysed",
  #      ylab="# morpho-species detected in all 3 mesophotic ARMS of P50ARMS3",
  #      ylim=c(1,100))
  # 
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  # 
  # 
  # text(26,
  #      10,
  #      paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
  #      cex = 0.85)
  # 
  # dev.off()
  # #### Spec curve by ARMS (instead of plates) ####
  # 
  # 
  # div_alpha_name <- paste0("div_alpha_ARMS.pdf")
  # div_alpha_path <- here::here("outputs/", div_alpha_name)
  # pdf(file =  div_alpha_path, width = 11, height = 8.5)
  # 
  # data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  # meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  # 
  # #### All MESO ARMS combined ####
  # 
  # data_mean_P50ARMS <- subset(data_mean, meta_mean$campain == "P50ARMS")
  # 
  # s <- vegan::specaccum(data_mean_P50ARMS, method = "random", permutations = 999,
  #                       conditioned =TRUE)
  # pool <- vegan::specpool(data_mean_P50ARMS)
  # 
  # plot(s, 
  #      ci.type="poly",
  #      col="blue", 
  #      lwd=2,
  #      ci.lty=0,
  #      ci.col="lightblue",
  #      xlab="Number of plates analysed",
  #      ylab="# morpho-species detected in all 9 mesophotic ARMS",
  #      ylim=c(1,100))
  # 
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  # 
  # 
  # text(5,
  #      10,
  #      paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
  #      cex = 0.85)
  # 
  # dev.off()
  #### 
  
  return(div_alpha_path)
}