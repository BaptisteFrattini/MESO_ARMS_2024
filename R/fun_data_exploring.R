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

  ff <- ggplot(data_div, aes(x = fct_relevel(campain, "RUNARMS", "P50ARMS"), y = S)) +
    geom_boxplot(fill =  c("lightblue","darkblue")) +
    labs(title = "",
         x = "Depth",
         y = "Average species richness in an ARMS") +
    theme_classic()
  
  gg <- ggplot(data_div, aes(x = fct_relevel(triplicat, "RUNARMS1","RUNARMS5","RUNARMS9", "P50ARMS1", "P50ARMS2", "P50ARMS3"), y = S)) +
    geom_boxplot(fill =  c(rep("lightblue",3),rep("darkblue",3))) +
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
  a <- paste0("stress = ",round(fit$stress, 3))
  text(0.8,1.1, a)
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
  a <- paste0("stress = ",round(fit2$stress, 3))
  text(0.8,0.8, a)
  dev.off()
  
  
  return(NULL)
}