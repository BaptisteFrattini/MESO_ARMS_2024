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
  
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
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
  
  #### Pour tous les sites ####
  
  # Trouver les MSP présentes dans les deux colonnes (intersection)
  msp_deep_shallow <- intersect(deep_msp, shallow_msp)

  
  only_deep_msp <- setdiff(deep_msp, msp_deep_shallow)
  
  only_shallow_msp <- setdiff(shallow_msp, msp_deep_shallow)
  
  devtools::install_github("yanlinlin82/ggvenn")
  library(ggvenn)
  
  x <- list(Deep = deep_msp,
            Shallow = shallow_msp)
  
  Venn_depth_path <- here::here("outputs/Venn/Venn_depth.png")
  
  venn <- ggvenn(x, fill_color = c("navy", "lightblue"), text_size = 5.5)
  
  ggsave(filename =  Venn_depth_path, plot = venn , width = 4, height = 4)
  
  
  #### Pour Saint Leu uniquement ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
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
  
  Venn_leu_path <- here::here("outputs/Venn/Venn_depth_Saint_Leu.png")
  
  venn_leu <- ggvenn(x, fill_color = c("navy", "lightblue"), text_size = 5.5)
  
  ggsave(filename =  Venn_leu_path, plot = venn_leu , width = 4, height = 4)
  
  
  
  #### Pour cap la houssaye uniquement ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
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
  
  Venn_cap_path <- here::here("outputs/Venn/Venn_depth_Cap.png")
  
  venn_cap <- ggvenn(x, fill_color = c("navy", "lightblue"), text_size = 5.5)
  
  ggsave(filename =  Venn_cap_path, plot = venn_cap , width = 4, height = 4)
  
  #### Pour grand bois uniquement ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
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
  
  Venn_gd_bois_path <- here::here("outputs/Venn/Venn_depth_Grd_Bois.png")
  
  
  venn_gd_bois <- ggvenn(x, fill_color = c("navy", "lightblue"), text_size = 5.5)
  
  ggsave(filename =  Venn_gd_bois_path, plot = venn_gd_bois , width = 4, height = 4)
  
  
  #### Entre Reunion et Rodrigues ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean_masc <- subset(data_mean, meta_mean$campain != "P50ARMS")
  meta_mean_masc <- subset(meta_mean, meta_mean$campain != "P50ARMS")
  
  data_island <- data_mean_masc %>% 
    group_by(meta_mean_masc$island) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_island <- as.data.frame(data_island)
  row.names(data_island) <- data_island$`meta_mean$island`
  data_island <- data_island[,-1]
  
  data_island_pa <- vegan::decostand(data_island, "pa")
  
  data_island_pa <- data.frame(t(data_island_pa))
  data_island_pa <- data.frame(msp = rownames(data_island_pa),
                              reunion = data_island_pa$X1,
                              rodrigues = data_island_pa$X2)
  
  reunion_msp <- data_island_pa$msp[data_island_pa["reunion"] == 1]
  rodrigues_msp <- data_island_pa$msp[data_island_pa["rodrigues"] == 1]
  
  # Trouver les MSP présentes dans les deux colonnes (intersection)
  msp_reunion_rodrigues <- intersect(reunion_msp, rodrigues_msp)
  
  
  only_reunion_msp <- setdiff(reunion_msp, msp_reunion_rodrigues)
  
  only_rodrigues_msp <- setdiff(rodrigues_msp, msp_reunion_rodrigues)
  
  
  x <- list(Reunion = reunion_msp,
            Rodrigues = rodrigues_msp)
  
  Venn_Reu_Rod_path <- here::here("outputs/Venn/Venn_island.png")
  
  library(ggvenn)
  Venn_Reu_Rod <- ggvenn(x, fill_color = c("navy", "lightblue"), text_size = 5.5)
  
  ggsave(filename =  Venn_Reu_Rod_path, plot = Venn_Reu_Rod , width = 4, height = 4)
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean_masc <- subset(data_mean, meta_mean$campain != "P50ARMS")
  meta_mean_masc <- subset(meta_mean, meta_mean$campain != "P50ARMS")
  
  data_island <- data_mean_masc %>% 
    group_by(meta_mean_masc$island) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_island <- as.data.frame(data_island)
  row.names(data_island) <- data_island$`meta_mean$island`
  data_island <- data_island[,-1]
  
  data_island_pa <- vegan::decostand(data_island, "pa")
  
  data_island_pa <- data.frame(t(data_island_pa))
  data_island_pa <- data.frame(msp = rownames(data_island_pa),
                              reunion = data_island_pa$X1,
                              rodrigues = data_island_pa$X2)
  
  reunion_msp <- data_island_pa$msp[data_island_pa["reunion"] == 1]
  rodrigues_msp <- data_island_pa$msp[data_island_pa["rodrigues"] == 1]
  
  # Trouver les MSP présentes dans les deux colonnes (intersection)
  msp_reunion_rodrigues <- intersect(reunion_msp, rodrigues_msp)
  
  
  only_reunion_msp <- setdiff(reunion_msp, msp_reunion_rodrigues)
  
  only_rodrigues_msp <- setdiff(rodrigues_msp, msp_reunion_rodrigues)
  
  
  x <- list(Reunion = reunion_msp,
            Rodrigues = rodrigues_msp)
  
  Venn_Reu_Rod_path <- here::here("outputs/Venn/Venn_island.png")
  
  library(ggvenn)
  Venn_Reu_Rod <- ggvenn(x, fill_color = c("navy", "lightblue"), text_size = 5.5)
  
  ggsave(filename =  Venn_Reu_Rod_path, plot = Venn_Reu_Rod , width = 4, height = 4)
  
  
  #### Entre les trois ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  
  data <- data_mean %>% 
    group_by(meta_mean$campain) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data <- as.data.frame(data)
  row.names(data) <- data$`meta_mean$campain`
  data <- data[,-1]
  
  data_pa <- vegan::decostand(data, "pa")
  
  data_pa <- data.frame(t(data_pa))
  data <- data_pa
  
  
 
  
  library(eulerr)

  
  # Create the Euler diagram
  fit <- euler(c(
    "P50ARMS" = sum(data$P50ARMS == 1 & data$RODARMS == 0 & data$RUNARMS == 0),
    "RODARMS" = sum(data$P50ARMS == 0 & data$RODARMS == 1 & data$RUNARMS == 0),
    "RUNARMS" = sum(data$P50ARMS == 0 & data$RODARMS == 0 & data$RUNARMS == 1),
    "P50ARMS&RODARMS" = sum(data$P50ARMS == 1 & data$RODARMS == 1 & data$RUNARMS == 0),
    "P50ARMS&RUNARMS" = sum(data$P50ARMS == 1 & data$RODARMS == 0 & data$RUNARMS == 1),
    "RODARMS&RUNARMS" = sum(data$P50ARMS == 0 & data$RODARMS == 1 & data$RUNARMS == 1),
    "P50ARMS&RODARMS&RUNARMS" = sum(data$P50ARMS == 1 & data$RODARMS == 1 & data$RUNARMS == 1)
  ))
  
  count <- fit$original.values
  
  total = sum(count)
  
  percentages <- (count/total)*100 
  
  sum(percentages)
  
  labels <- paste0(count, " (", round(percentages, 1), "%)")
  labels[1] <- paste0("P50ARMS \n",labels[1])
  labels[2] <- paste0("RODARMS \n",labels[2])
  labels[3] <- paste0("RUNARMS \n",labels[3])
  
  venn_name <- paste0("Venn_full.pdf")
  venn_path <- paste0("outputs/Venn/", venn_name)
  pdf(file =  venn_path)
  plot(fit, fills = list(fill = c("dodgerblue2", "darkorange", "lightgreen")), labels = labels, main = "")
  dev.off()
  
  data$msp <- rownames(data)
  
  
  reunion_msp <- data$msp[data["RUNARMS"] == 1]
  rodrigues_msp <- data$msp[data["RODARMS"] == 1]
  p50_msp <- data$msp[data["P50ARMS"] == 1]
  
  # Calcul des différents ensembles
  only_reunion <- setdiff(reunion_msp, union(rodrigues_msp, p50_msp))
  only_rodrigues <- setdiff(rodrigues_msp, union(reunion_msp, p50_msp))
  only_p50 <- setdiff(p50_msp, union(reunion_msp, rodrigues_msp))
  
  reunion_rodrigues_only <- setdiff(intersect(reunion_msp, rodrigues_msp), p50_msp)
  reunion_p50_only <- setdiff(intersect(reunion_msp, p50_msp), rodrigues_msp)
  rodrigues_p50_only <- setdiff(intersect(rodrigues_msp, p50_msp), reunion_msp)
  
  intersection_all <- intersect(intersect(reunion_msp, rodrigues_msp), p50_msp)
  
  # Chargement de la librairie pour écrire dans un fichier Excel
  install.packages("writexl")
  library(writexl)
  
  # Création d'une liste de data frames pour chaque colonne
  results <- list(
    "Only Reunion" = only_reunion,
    "Only Rodrigues" = only_rodrigues,
    "Only P50" = only_p50,
    "Reunion & Rodrigues (sans P50)" = reunion_rodrigues_only,
    "Reunion & P50 (sans Rodrigues)" = reunion_p50_only,
    "Rodrigues & P50 (sans Reunion)" = rodrigues_p50_only,
    "Intersection All" = intersection_all
  )
  
  # Créer un data frame avec des colonnes de longueur égale en remplissant avec des NA
  max_length <- max(sapply(results, length))
  df <- data.frame(lapply(results, function(x) { c(x, rep(NA, max_length - length(x))) }))
  
  # Écriture du data frame dans un fichier Excel
  write_xlsx(df, "outputs/Venn/combinations_msp.xlsx")
  
  
  
  
  
  # MSP présente à la reunion ET à rodrigues
  msp_reunion_rodrigues <- intersect(reunion_msp, rodrigues_msp)
  # MSp présente à la reunion, a rodrigues, mais pas a P50
  msp_reunion_rodrigues_seulement <- setdiff(msp_reunion_rodrigues, p50_msp)
  
  
  # MSP présente à la reunion , à rodrigues et à P50
  msp_all <- intersect(msp_reunion_rodrigues, p50_msp)
  
  
  
  
  only_rodrigues_msp <- setdiff(rodrigues_msp, msp_reunion_rodrigues)

  #### return ####
  
  
  
  fin <- cowplot::plot_grid(venn, venn_cap, venn_leu, venn_gd_bois, Venn_Reu_Rod,
                            ncol = 2,
                            nrow = 3,
                            labels = c("All sites", "Cap La Houssaye", "Saint-Leu", "Grand Bois", "Rodrigues"),
                            rel_widths = c(0.6, 0.6),  # Adjust as needed
                            rel_heights = c(0.6, 0.6)) # Adjust as needed
  
  
  full_venn_path <- here::here("outputs/Venn/full_venn.png")
  
  ggsave(filename =  full_venn_path, plot = fin , width = 11, height = 11)
  
  Venn_tot <- c(Venn_depth_path,Venn_cap_path, Venn_leu_path, Venn_gd_bois_path, full_venn_path)
  
  return(Venn_tot)
  
}
