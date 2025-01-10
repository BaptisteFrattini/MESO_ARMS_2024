#' PERMANOVA and SIMPER
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_permanova <- function(data_and_meta_clean){
  
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  # PERMANOVA ####
  ## Effet de la profondeur ####
  
  data_mean_RUN <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean_RUN <- subset(meta_mean, meta_mean$island == "Reunion")
  
  matrix_bc <- vegan::vegdist(data_mean_RUN, "bray")
  
  perm <- vegan::adonis2(data_mean_RUN ~ campain, strata = meta_mean_RUN$site, data = meta_mean_RUN , method = "bray", permutations = 99999)
  
  sp_contrib <- summary(vegan::simper(data_mean_RUN, meta_mean_RUN $campain))
  
  ## Effet de l'ile ####
  
  data_mean_shallow <- data_mean[-c(1:9),]
  meta_mean_shallow <- meta_mean[-c(1:9),]
  
  perm <- vegan::adonis2(data_mean_shallow ~ island, data = meta_mean_shallow, method = "bray", permutations = 99999)
  
  sp_contrib <- summary(vegan::simper(data_mean_shallow, meta_mean_shallow$campain))
  
  # PERMDISP ####
  
  dis.bray <- vegan::vegdist(data_mean, "bray")
  # dis.jacc <- vegan::vegdist(data_mean_pa, "jaccard")
  
  a <- vegan::betadisper(
    dis.bray,
    meta_mean$triplicat,
    type = "median",
    bias.adjust = FALSE,
    sqrt.dist = FALSE,
    add = TRUE
  )
  boxplot(a)
  anova(a)
  plot(a)
  
  ?vegan::betadisper
  
  b <- vegan::betadisper(
    dis.bray,
    meta_mean$campain,
    type = "median",
    bias.adjust = FALSE,
    sqrt.dist = FALSE,
    add = TRUE
  )
  boxplot(b)
  
  anova(b)
  
  return(NULL)
}