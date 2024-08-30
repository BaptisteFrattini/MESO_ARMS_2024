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
  
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
  matrix_bc <- vegan::vegdist(data_mean, "bray")
  
  perm <- vegan::adonis2(data_mean ~ campain, strata = meta_mean$site, data = meta_mean, method = "bray", permutations = 99999)
  
  sp_contrib <- summary(vegan::simper(data_mean, meta_mean$campain))

  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[-c(1:9),]
  meta_mean <- meta_mean[-c(1:9),]
  
  perm <- vegan::adonis2(data_mean ~ island, data = meta_mean, method = "bray", permutations = 99999)
  
  sp_contrib <- summary(vegan::simper(data_mean, meta_mean$campain))
  
  
  return(NULL)
}