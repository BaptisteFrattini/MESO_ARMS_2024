#' Venn diagrams
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to the venn diag
#' @export
#'
fun_venn_diag <- function(data_and_meta_clean){
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  library(dplyr)
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_depth <- data_mean %>% 
    group_by(meta_mean$campain) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_depth <- as.data.frame(data_depth)
  row.names(data_depth) <- data_depth$`meta_mean$campain`
  data_depth <- data_depth[,-1]
  
  data_depth_pa <- vegan::decostand(data_depth, "pa")
  
  data_depth_pa <- data.frame(t(data_depth_pa))
  data_depth_pa <- data.frame(msp = rownames(data_depth_pa),
                              deep = data_depth_pa$P50ARMS,
                              shallow = data_depth_pa$RUNARMS)
  
  deep_msp <- data_depth_pa$msp[data_depth_pa["deep"] == 1]
  length(deep_msp)
  shallow_msp <- data_depth_pa$msp[data_depth_pa["shallow"] == 1]
  length(shallow_msp)
  
  # Trouver les MSP présentes dans les deux colonnes (intersection)
  msp_deep_shallow <- intersect(deep_msp, shallow_msp)

  
  only_deep_msp <- setdiff(deep_msp, msp_deep_shallow)
  
  only_shallow_msp <- setdiff(shallow_msp, msp_deep_shallow)
  
  devtools::install_github("yanlinlin82/ggvenn")
  library(ggvenn)
  
  x <- list(Deep = deep_msp,
            Shallow = shallow_msp)
  
  Venn_depth_path <- here::here("outputs/Venn_depth.pdf")
  
  venn <- ggvenn(x, fill_color = c("navy", "lightblue"))
  
  ggsave(filename =  Venn_depth_path, plot = venn , width = 6, height = 6)
  
  
  #### Pour Saint Leu uniquement ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[grepl("RUNARMS5|P50ARMS2", rownames(data_mean)), ]
  meta_mean <-meta_mean[grepl("RUNARMS5|P50ARMS2", meta_mean$triplicat), ]
  
  data_depth <- data_mean %>% 
    group_by(meta_mean$campain) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_depth <- as.data.frame(data_depth)
  row.names(data_depth) <- data_depth$`meta_mean$campain`
  data_depth <- data_depth[,-1]
  
  data_depth_pa <- vegan::decostand(data_depth, "pa")
  
  data_depth_pa <- data.frame(t(data_depth_pa))
  data_depth_pa <- data.frame(msp = rownames(data_depth_pa),
                              deep = data_depth_pa$P50ARMS,
                              shallow = data_depth_pa$RUNARMS)
  
  deep_msp <- data_depth_pa$msp[data_depth_pa["deep"] == 1]
  shallow_msp <- data_depth_pa$msp[data_depth_pa["shallow"] == 1]
  
  # Trouver les MSP présentes dans les deux colonnes (intersection)
  msp_deep_shallow <- intersect(deep_msp, shallow_msp)
  
  
  only_deep_msp <- setdiff(deep_msp, msp_deep_shallow)
  
  only_shallow_msp <- setdiff(shallow_msp, msp_deep_shallow)
  
  devtools::install_github("yanlinlin82/ggvenn")
  library(ggvenn)
  
  x <- list(Deep = deep_msp,
            Shallow = shallow_msp)
  
  Venn_depth_path <- here::here("outputs/Venn_depth_Saint_Leu.pdf")
  
  venn <- ggvenn(x, fill_color = c("navy", "lightblue"))
  
  ggsave(filename =  Venn_depth_path, plot = venn , width = 6, height = 6)
  
  
  
  #### Pour cap la houssaye uniquement ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[grepl("RUNARMS1|P50ARMS1", rownames(data_mean)), ]
  meta_mean <-meta_mean[grepl("RUNARMS1|P50ARMS1", meta_mean$triplicat), ]
  
  data_depth <- data_mean %>% 
    group_by(meta_mean$campain) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_depth <- as.data.frame(data_depth)
  row.names(data_depth) <- data_depth$`meta_mean$campain`
  data_depth <- data_depth[,-1]
  
  data_depth_pa <- vegan::decostand(data_depth, "pa")
  
  data_depth_pa <- data.frame(t(data_depth_pa))
  data_depth_pa <- data.frame(msp = rownames(data_depth_pa),
                              deep = data_depth_pa$P50ARMS,
                              shallow = data_depth_pa$RUNARMS)
  
  deep_msp <- data_depth_pa$msp[data_depth_pa["deep"] == 1]
  shallow_msp <- data_depth_pa$msp[data_depth_pa["shallow"] == 1]
  
  # Trouver les MSP présentes dans les deux colonnes (intersection)
  msp_deep_shallow <- intersect(deep_msp, shallow_msp)
  
  
  only_deep_msp <- setdiff(deep_msp, msp_deep_shallow)
  
  only_shallow_msp <- setdiff(shallow_msp, msp_deep_shallow)
 
  x <- list(Deep = deep_msp,
            Shallow = shallow_msp)
  
  Venn_depth_path <- here::here("outputs/Venn_depth_Cap.pdf")
  
  venn <- ggvenn(x, fill_color = c("navy", "lightblue"))
  
  ggsave(filename =  Venn_depth_path, plot = venn , width = 6, height = 6)
  
  #### Pour grand bois uniquement ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[grepl("RUNARMS9|P50ARMS3", rownames(data_mean)), ]
  meta_mean <-meta_mean[grepl("RUNARMS9|P50ARMS3", meta_mean$triplicat), ]
  
  data_depth <- data_mean %>% 
    group_by(meta_mean$campain) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_depth <- as.data.frame(data_depth)
  row.names(data_depth) <- data_depth$`meta_mean$campain`
  data_depth <- data_depth[,-1]
  
  data_depth_pa <- vegan::decostand(data_depth, "pa")
  
  data_depth_pa <- data.frame(t(data_depth_pa))
  data_depth_pa <- data.frame(msp = rownames(data_depth_pa),
                              deep = data_depth_pa$P50ARMS,
                              shallow = data_depth_pa$RUNARMS)
  
  deep_msp <- data_depth_pa$msp[data_depth_pa["deep"] == 1]
  shallow_msp <- data_depth_pa$msp[data_depth_pa["shallow"] == 1]
  
  # Trouver les MSP présentes dans les deux colonnes (intersection)
  msp_deep_shallow <- intersect(deep_msp, shallow_msp)
  
  
  only_deep_msp <- setdiff(deep_msp, msp_deep_shallow)
  
  only_shallow_msp <- setdiff(shallow_msp, msp_deep_shallow)
  
  
  x <- list(Deep = deep_msp,
            Shallow = shallow_msp)
  
  Venn_depth_path <- here::here("outputs/Venn_depth_Grd_Bois.pdf")
  
  
  venn <- ggvenn(x, fill_color = c("navy", "lightblue"))
  
  ggsave(filename =  Venn_depth_path, plot = venn , width = 6, height = 6)
  
  
  
  return(Venn_depth_path)
  
}