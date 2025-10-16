#' Mantel test
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_mantel <- function(data_and_meta_clean_fullsites, gps_sites){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites")
  # gps_sites = targets::tar_read("data_gps_sites")
  library(betapart)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[meta_mean$campain %in% c("RODARMS", "RUNARMS"), ]
  meta_mean <- meta_mean[meta_mean$campain %in% c("RODARMS", "RUNARMS"), ]
  
  data_mean <- data_mean[, colSums(data_mean, na.rm = TRUE) != 0]
  
  msp_list <- names(data_mean)
  # msp_list_filter <- msp_list[!msp_list %in% c("No_Recruitment",
  #                                              "Sediments")]

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



  data_mean <- data_mean[, msp_list_filter]
  
  
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE)
  data_gps <- data_gps[-c(1:3),]
  
  data_gps <-  data_gps %>%
    slice(rep(1:n(), each = 3))
  data_gps <- data_gps[-9,]
  
  data_gps$arms <- meta_mean$arms

  
  coords <- data_gps[, c("Longitude", "Latitude")]
  rownames(coords) <- data_gps$arms
  
  matrix.dist <- geosphere::distm(coords)
  row.names(matrix.dist) <- data_gps$arms
  colnames(matrix.dist) <- data_gps$arms
  matrix.dist <- matrix.dist
  matrix.dist <- as.data.frame(matrix.dist)
  matrix.dist <- as.matrix(matrix.dist)
  matrix.dist <- as.dist(matrix.dist)
  matrix.dist[matrix.dist == 0] <- 3


  class(matrix.dist)

  
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

  ##### bray #####

  aa = as.vector(mat.bc)
  tt = log10(as.vector(matrix.dist))
  #new data frame with vectorized distance matrices
  mat = data.frame(aa,tt)

  mat$ile <- ifelse(mat$tt > 5.0, "btw_islands", "within_islands")
  
  mm1 = ggplot2::ggplot(mat, aes(y = aa, x = tt)) +
    geom_point(
      data = subset(mat, ile == "within_islands"),  # points à triangle uniquement pour btw_islands
      shape = 16,  # 
      size = 3,
      alpha = 0.5
    ) +
    geom_point(
      data = subset(mat, ile == "btw_islands"),  # points à triangle uniquement pour btw_islands
      shape = 17,  # triangle
      size = 3,
      alpha = 0.5
    ) +
    labs(x = "log10(Distance géographique)",
         y = "Dissimilarité de Bray-Curtis") +
    geom_smooth(method = "lm",
                colour = "red",
                alpha = 0.2,
                fill = "red", se = FALSE) +
    theme( axis.text.x = element_text(colour = "black",
                                      size = 12),
           axis.text.y = element_text(size = 11,
                                      colour = "black"),
           axis.title = element_text(size = 14,
                                     colour = "black"),
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA,
                                       colour = "black")) + 
    theme(legend.position = "none")
  
  mm1
  
  
  log_mat_dist <- as.dist(log10(matrix.dist))
  
  correl = vegan::mantel(log_mat_dist, as.dist(mat.bc), method = "pearson")
  p <- correl$signif
  if (p < 0.001) {
    significativite <- "***"
  } else if (p >= 0.001 && p < 0.01) {
    significativite <- "**"
  } else if (p >= 0.01 && p < 0.05) {
    significativite <- "*"
  } else {
    significativite <- " NS"
  }
  
  mm1 = mm1 + annotate(geom = "text",  x = min(mat$tt), y = 0.9, label = paste0("Mantel R = ",  round(correl$statistic, 3), "; p = ", correl$signif,significativite),
                       color = "black", size = 3.5, hjust = 0, vjust = 1)
    # scale_x_continuous(breaks = seq(0, 100, 10), minor_breaks = seq(0, 100, 1)) +
    # ylim(0,0.7)
  # Modèle linéaire simple : area en fonction de cpc
  mod <- lm(aa ~ tt, data = mat)
  coef(mod)
  cor <- cor(mat$aa, mat$tt, method = "spearman")
  
  # Afficher le résumé du modèle
  summary(mod)
  # Extraire R²
  r2 <- summary(mod)$r.squared
  
  print(paste("R² du modèle linéaire :", round(r2, 3)))
  
  shapiro.test(residuals(mod))
  
  qqnorm(residuals(mod))
  qqline(residuals(mod), col = "red")
  # plot(mod)
  
  
  ####Jacc ####
  
  bb = as.vector(mat.jacc)
  tt = log10(as.vector(matrix.dist))
  #new data frame with vectorized distance matrices
  mat = data.frame(bb,tt)
  
  mat$ile <- ifelse(mat$tt > 5.0, "btw_islands", "within_islands")
  
  mm1_bis = ggplot2::ggplot(mat, aes(y = bb, x = tt)) +
    geom_point(
      data = subset(mat, ile == "within_islands"),  # points à triangle uniquement pour btw_islands
      shape = 16,  # 
      size = 3,
      alpha = 0.5
    ) +
    geom_point(
      data = subset(mat, ile == "btw_islands"),  # points à triangle uniquement pour btw_islands
      shape = 17,  # triangle
      size = 3,
      alpha = 0.5
    ) +
    labs(x = "log10(Distance géographique)",
         y = "Dissimilarité de Jaccard") +
    geom_smooth(method = "lm",
                colour = "red",
                alpha = 0.2,
                fill = "red", se = FALSE) +
    theme( axis.text.x = element_text(colour = "black",
                                      size = 12),
           axis.text.y = element_text(size = 11,
                                      colour = "black"),
           axis.title = element_text(size = 14,
                                     colour = "black"),
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA,
                                       colour = "black")) + 
    theme(legend.position = "none")
  
  mm1_bis
  
  correl = vegan::mantel(log_mat_dist, as.dist(mat.jacc), method = "pearson")
  p <- correl$signif
  if (p < 0.001) {
    significativite <- "***"
  } else if (p >= 0.001 && p < 0.01) {
    significativite <- "**"
  } else if (p >= 0.01 && p < 0.05) {
    significativite <- "*"
  } else {
    significativite <- " NS"
  }
  
  mm1_bis = mm1_bis + annotate(geom = "text",  x = min(mat$tt), y = 0.90, label = paste0("Mantel R = ",  round(correl$statistic, 3), "; p = ", correl$signif,significativite),
                       color = "black", size = 3.5, hjust = 0, vjust = 1)
  
  # Modèle linéaire simple : area en fonction de cpc
  mod <- lm(bb ~ tt, data = mat)
  coef(mod)
  
  # Afficher le résumé du modèle
  summary(mod)
  # Extraire R²
  r2 <- summary(mod)$r.squared
  
  print(paste("R² du modèle linéaire :", round(r2, 3)))
  
  shapiro.test(residuals(mod))
  
  #### A tester avec une separation espèce rare et espèces communes ####
  
  
  species_occurence <- colSums(data_mean)/nrow(data_mean)
  species_occurence <- data.frame(f = species_occurence,
                                  Species = names(species_occurence))
  threshold <- quantile(species_occurence$f, 0.25)
  
  rare_species <- species_occurence[species_occurence$f <= threshold,]
  length(rare_species$Species)
  common_species <- species_occurence[species_occurence$f >= threshold,]
  length(common_species$Species)
  
  data_mean_rare <- data_mean[,rare_species$Species]
  data_mean_common <- data_mean[,common_species$Species]
  
  
  mat.bc.rare <- vegan::vegdist(data_mean_rare, dist = "bray")
  data_mean_rare_pa <- vegan::decostand(data_mean_rare, method = "pa")
  B.pair.pa.rare <- betapart::beta.pair(data_mean_rare_pa, index.family = "jaccard")
  mat.turn.rare <- B.pair.pa.rare$beta.jtu
  mat.nest.rare <- B.pair.pa.rare$beta.jne
  mat.jacc.rare <- B.pair.pa.rare$beta.jac
  
  mat.bc.common <- vegan::vegdist(data_mean_common, dist = "bray")
  data_mean_common_pa <- vegan::decostand(data_mean_common, method = "pa")
  B.pair.pa.common <- betapart::beta.pair(data_mean_common_pa, index.family = "jaccard")
  mat.turn.common <- B.pair.pa.common$beta.jtu
  mat.nest.common <- B.pair.pa.common$beta.jne
  mat.jacc.common <- B.pair.pa.common$beta.jac
  
  
  ##### jacc #####
  
  # aa = log10(as.vector(mat.bc.common))
  # bb = log10(as.vector(mat.bc.rare))
  # cc = log10(as.vector(mat.bc))
  # tt = log10(as.vector(matrix.dist))
  
  aa = as.vector(mat.bc.common)
  bb = as.vector(mat.bc.rare)
  tt = log10(as.vector(matrix.dist))
  
  #new data frame with vectorized distance matrices
  mat <- data.frame(
    value = c(aa, bb),
    distance = rep(tt, 2),
    group = rep(c("common", "rare"), each = length(tt))
  )
  
  mat$ile <- ifelse(mat$distance > 5.0, "btw_islands", "within_islands")
  
  
  
  # Plot
  library(ggplot2)
  
  mm2 <- ggplot(mat, aes(x = distance, y = value, color = group)) +
    geom_point(
      data = subset(mat, ile == "within_islands"),  # points à triangle uniquement pour btw_islands
      shape = 16,  # 
      size = 3,
      alpha = 0.5
    ) +
    geom_point(
      data = subset(mat, ile == "btw_islands"),  # points à triangle uniquement pour btw_islands
      shape = 17,  # triangle
      size = 3,
      alpha = 0.5
    ) +
    geom_smooth(method = "lm", alpha = 0.2, aes(fill = group), se = FALSE) +
    labs(x = "log10(Distance géographique)",
         y = "Dissimilarité de Bray-Curtis",
         color = "Group",
         fill = "Group") +
    theme(
      axis.text.x = element_text(colour = "black", size = 12),
      axis.text.y = element_text(colour = "black", size = 11),
      axis.title = element_text(colour = "black", size = 14),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black")
    )+
    scale_color_manual(values = c("common" = "darkblue", "rare" = "darkorange"))+ 
    theme(legend.position = "none")
  
  mm2
  
  
  correl = vegan::mantel(log_mat_dist, as.dist(mat.bc.common), method = "pearson")
  p <- correl$signif
  if (p < 0.001) {
    significativite <- "***"
  } else if (p >= 0.001 && p < 0.01) {
    significativite <- "**"
  } else if (p >= 0.01 && p < 0.05) {
    significativite <- "*"
  } else {
    significativite <- " NS"
  }
  
  mm2 = mm2 + annotate(geom = "text",  x = min(mat$distance), y = 0.9, label = paste0("Mantel R = ",  round(correl$statistic, 3), "; p = ", correl$signif,significativite),
                               color = "darkblue", size = 3.5, hjust = 0, vjust = 1)
  
  mod <- lm(aa ~ tt, data = mat)
  
  mod_rare <- lm(bb~ tt, data = mat)
  
  coef(mod)
  coef(mod_rare)
  
  r2 <- summary(mod)$r.squared
  r2
  
  r2_rare <- summary(mod_rare)$r.squared
  r2_rare
  
  summary(mod)
 
  #### Jacc ####
  
  cc = as.vector(mat.jacc.common)
  dd = as.vector(mat.jacc.rare)
  tt = log10(as.vector(matrix.dist))
  
  #new data frame with vectorized distance matrices
  mat_jacc <- data.frame(
    value = c(cc, dd),
    distance = rep(tt, 2),
    group = rep(c("common", "rare"), each = length(tt))
  )
  
  mat_jacc <- mat_jacc %>%
    mutate(ile = ifelse(distance > 5.0, "btw_islands", "within_islands"))
  
  mm3 <- ggplot(mat_jacc, aes(x = distance, y = value, color = group)) +
    geom_point(
      data = subset(mat_jacc, ile == "within_islands"),  # points à triangle uniquement pour btw_islands
      shape = 16,  # 
      size = 3,
      alpha = 0.5
    ) +
    geom_point(
      data = subset(mat_jacc, ile == "btw_islands"),  # points à triangle uniquement pour btw_islands
      shape = 17,  # triangle
      size = 3,
      alpha = 0.5
    ) +
    geom_smooth(method = "lm", aes(fill = group), se = FALSE) +
    labs(x = "log10(Distance géographique)",
         y = "Dissimilarité de Jaccard",
         color = "Group",
         fill = "Group") +
    theme(
      axis.text.x = element_text(colour = "black", size = 12),
      axis.text.y = element_text(colour = "black", size = 11),
      axis.title = element_text(colour = "black", size = 14),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"))+
    scale_color_manual(values = c("common" = "darkblue", "rare" = "darkorange"))+ 
    theme(legend.position = "none")
  
  
  mm3
  
  correl = vegan::mantel(log_mat_dist, as.dist(mat.jacc.common), method = "pearson")
  p <- correl$signif
  if (p < 0.001) {
    significativite <- "***"
  } else if (p >= 0.001 && p < 0.01) {
    significativite <- "**"
  } else if (p >= 0.01 && p < 0.05) {
    significativite <- "*"
  } else {
    significativite <- " NS"
  }
  
  mm3 = mm3 + annotate(geom = "text",  x = min(mat$distance), y = 0.90, label = paste0("Mantel R = ",  round(correl$statistic, 3), "; p = ", correl$signif,significativite),
                       color = "darkblue", size = 3.5, hjust = 0, vjust = 1)
  
  mod <- lm(cc ~ tt)
  
  mod_rare <- lm(bb~ tt)
  
  coef(mod)
  coef(mod_rare)
  
  r2 <- summary(mod)$r.squared
  r2
  
  r2_rare <- summary(mod_rare)$r.squared
  r2_rare
  
  summary(mod)
  
  fig1 <- cowplot::plot_grid(mm1, mm1_bis, mm2, mm3, ncol = 2, nrow = 2)
  
  
  fin_1 <- cowplot::ggdraw(fig1) +
    theme(plot.margin = margin(t = 1, r = 3, b = 1, l = 3, unit = "cm"))
  
  path_to_mantel_island <- paste0("outputs/Graph_mantel_island_mantel.pdf")
  ggsave(filename =  path_to_mantel_island , plot = fin_1, width = 10.5, height = 9)
  
  path_to_mantel_island <- paste0("outputs/Graph_mantel_island_mantel.png")
  ggsave(filename =  path_to_mantel_island , plot = mm3, width = 10.5, height = 9)
  
  #### conception d'une boucle pour voir comment évolue R2 et coef en fonction du threshold utilisé ####
  

  library(ggbreak)
  
  cc = as.vector(mat.jacc.common)
  dd = as.vector(mat.jacc.rare)
  tt = as.vector(matrix.dist)/1000
  
  #new data frame with vectorized distance matrices
  mat_jacc <- data.frame(
    value = c(cc, dd),
    distance = rep(tt, 2),
    group = rep(c("common", "rare"), each = length(tt))
  )
  
  p <- ggplot(data = mat_jacc, aes(x = distance, y = value, color = group)) +
    geom_point(size = 3, alpha = 0.5) +
    # geom_smooth(method = "lm", se = FALSE, aes(fill = group), alpha = 0.2) +
    scale_color_manual(values = c("common" = "darkblue", "rare" = "darkorange")) +
    labs(x = "Distance géographique(km)",
         y = "Dissimilarité de Jaccard",
         color = "Group",
         fill = "Group") + 
    theme(legend.position = "none")
  
  # Ajout d’une cassure sur l’axe Y (exemple pour sauter les valeurs entre 0.2 et 0.6)
  p_broken <- p + scale_x_break(c(53, 839), scales = 0.5)  # tu peux ajuster les seuils
  
  p_broken
  
  correl = vegan::mantel(matrix.dist, as.dist(mat.jacc.common), method = "pearson")
 
  
  z <- correl$signif
  if (z < 0.001) {
    significativite <- "***"
  } else if (z >= 0.001 && z < 0.01) {
    significativite <- "**"
  } else if (z >= 0.01 && z < 0.05) {
    significativite <- "*"
  } else {
    significativite <- " NS"
  }
  
  p_broken <- p_broken + annotate(geom = "text",  x = min(mat_jacc$distance), y = 0.90, label = paste0("Mantel R = ",  round(correl$statistic, 3), "; p = ", correl$signif,significativite),
                 color = "darkblue", size = 3.5, hjust = 0, vjust = 1)
  
  p_broken
  

  path_to_mantel_island <- paste0("outputs/mantel_bc_broken_jacc.pdf")
  ggsave(filename =  path_to_mantel_island , plot = p_broken, width = 6, height =4)
  
  cc = as.vector(mat.bc.common)
  dd = as.vector(mat.bc.rare)
  tt = as.vector(matrix.dist)/1000
  
  #new data frame with vectorized distance matrices
  mat_bc <- data.frame(
    value = c(cc, dd),
    distance = rep(tt, 2),
    group = rep(c("common", "rare"), each = length(tt))
  )
  
  t <- ggplot(data = mat_bc, aes(x = distance, y = value, color = group)) +
    geom_point(size = 3, alpha = 0.5) +
    # geom_smooth(method = "lm", se = FALSE, aes(fill = group), alpha = 0.2) +
    scale_color_manual(values = c("common" = "darkblue", "rare" = "darkorange")) +
    labs(x = "Distance géographique(km)",
         y = "Dissimilarité de Bray-Curtis",
         color = "Group",
         fill = "Group") + 
    theme(legend.position = "none")
    # theme_minimal()
  
  # Ajout d’une cassure sur l’axe Y (exemple pour sauter les valeurs entre 0.2 et 0.6)
  t_broken <- t + scale_x_break(c(53, 839), scales = 0.5)  # tu peux ajuster les seuils
  
  t_broken
  
  correl = vegan::mantel(matrix.dist, as.dist(mat.bc.common), method = "pearson")
  
  mat.bc.rare[is.na(mat.bc.rare)] <- 0
  correl_rare = vegan::mantel(matrix.dist, as.dist(mat.bc.rare), method = "pearson")
  
  z <- correl$signif
  if (z < 0.001) {
    significativite <- "***"
  } else if (z >= 0.001 && z < 0.01) {
    significativite <- "**"
  } else if (z >= 0.01 && z < 0.05) {
    significativite <- "*"
  } else {
    significativite <- " NS"
  }
  
  t_broken <- t_broken + annotate(geom = "text",  x = min(mat_bc$distance), y = 0.90, label = paste0("Mantel R = ",  round(correl$statistic, 3), "; p = ", correl$signif,significativite),
                                  color = "darkblue", size = 3.5, hjust = 0, vjust = 1)
  
  t_broken

  path_to_mantel_island <- paste0("outputs/mantel_bc_broken_bc.pdf")
  ggsave(filename =  path_to_mantel_island , plot = t_broken, width = 6, height =4)
  ######test####
  
  # Vecteur de seuils de rareté
  thresholds <- seq(0, 1, by = 0.01)
  
  # Liste pour stocker les résultats
  results <- data.frame(threshold = numeric(),
                        intercept = numeric(),
                        slope = numeric(),
                        r_squared = numeric())
  
  # Boucle sur les seuils
  for (rarity_thresh in thresholds) {
    
    # Calcul de l'occurrence moyenne des espèces
    species_occurence <- colSums(data_mean) / nrow(data_mean)
    species_occurence <- data.frame(f = species_occurence,
                                    Species = names(species_occurence))
    
    # Seuil basé sur le quantile
    threshold <- quantile(species_occurence$f, rarity_thresh)
    
    # Définition des espèces rares et communes
    rare_species <- species_occurence[species_occurence$f <= threshold, ]
    common_species <- species_occurence[species_occurence$f >= threshold, ]
    
    # Sous-ensembles de la matrice
    data_mean_rare <- data_mean[, rare_species$Species, drop = FALSE]
    data_mean_common <- data_mean[, common_species$Species, drop = FALSE]
    
    # Calcul des distances de Bray-Curtis
    mat.bc.rare <- vegan::vegdist(data_mean_rare, dist = "bray")
    mat.bc.common <- vegan::vegdist(data_mean_common, dist = "bray")
    mat.bc.total <- vegan::vegdist(data_mean, dist = "bray")
    
    # Vecteurs log10 des distances
    aa <- log10(as.vector(mat.bc.common))
    bb <- log10(as.vector(mat.bc.rare))
    cc <- log10(as.vector(mat.bc.total))
    tt <- log10(as.vector(matrix.dist))
    
    # Régression pour les distances des espèces communes
    mod <- lm(aa ~ tt)
    coef_mod <- coef(mod)
    r2 <- summary(mod)$r.squared
    
    # Stockage des résultats
    results <- rbind(results, data.frame(
      threshold = rarity_thresh,
      intercept = coef_mod[1],
      slope = coef_mod[2],
      r_squared = r2
    ))
  }
  
  # Affichage des résultats
  print(results)
  
  plot(results$threshold, results$slope)
  plot(results$threshold, results$r_squared)
  
  #### test avec jaccard ####
  
  # Vecteur de seuils de rareté
  thresholds <- seq(0, 1, by = 0.01)
  
  # Liste pour stocker les résultats
  results <- data.frame(threshold = numeric(),
                        intercept = numeric(),
                        slope = numeric(),
                        r_squared = numeric())
  
  # Boucle sur les seuils
  for (rarity_thresh in thresholds) {
    
    # Calcul de l'occurrence moyenne des espèces
    species_occurence <- colSums(data_mean) / nrow(data_mean)
    species_occurence <- data.frame(f = species_occurence,
                                    Species = names(species_occurence))
    
    # Seuil basé sur le quantile
    threshold <- quantile(species_occurence$f, rarity_thresh)
    
    # Définition des espèces rares et communes
    rare_species <- species_occurence[species_occurence$f <= threshold, ]
    common_species <- species_occurence[species_occurence$f >= threshold, ]
    
    # Sous-ensembles de la matrice
    data_mean_rare <- data_mean[, rare_species$Species, drop = FALSE]
    data_mean_common <- data_mean[, common_species$Species, drop = FALSE]
    
    # Calcul des distances de Bray-Curtis
    data_mean_rare <- vegan::decostand(data_mean_rare, "pa")
    data_mean_common <- vegan::decostand(data_mean_common, "pa")
    mat.jacc.rare <- vegan::vegdist(data_mean_rare, method = "jaccard")
    mat.jacc.common <- vegan::vegdist(data_mean_common, method = "jaccard")
    mat.jacc.total <- vegan::vegdist(data_mean, method = "jaccard")
    
    # Vecteurs log10 des distances
    aa <- as.vector(mat.jacc.common)
    bb <- as.vector(mat.jacc.rare)
    cc <- as.vector(mat.jacc.total)
    tt <- log10(as.vector(matrix.dist))
    
    # Régression pour les distances des espèces communes
    mod <- lm(aa ~ tt)
    coef_mod <- coef(mod)
    r2 <- summary(mod)$r.squared
    
    # Stockage des résultats
    results <- rbind(results, data.frame(
      threshold = rarity_thresh,
      intercept = coef_mod[1],
      slope = coef_mod[2],
      r_squared = r2
    ))
  }
  
  # Affichage des résultats
  print(results)
  
  plot(results$threshold, results$slope)
  plot(results$threshold, results$r_squared)
  
  
  
  
  return(NULL)
  
}

  
  
  