#' Testing relation btw temperature and community
#' 
#'
#' @return the path to the subseted raw data file
#' @export
#'
fun_temp_comm_RDA <- function(temp_extraction, temp_in_situ, data_and_meta_clean_fullsites){
# 
# temp_extraction = targets::tar_read("temp_extraction")
# temp_in_situ = targets::tar_read("temp_in_situ")
# data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites")
  library(dplyr)
  library(tidyr)
  library(ggplot2) 
  library(vegan)
  library(lubridate)
  library(ggrepel) 

  data_copernicus <- read.csv(temp_extraction[1])
  data_copernicus$date <- as.Date(data_copernicus$date)
  
  
  data_in_situ_shallow <- read.csv(temp_in_situ[grepl("shallow",temp_in_situ)])
  data_in_situ_deep <- read.csv(temp_in_situ[grepl("deep",temp_in_situ)])
  colnames(data_in_situ_deep) <- c("date", "T_1A", "T_1B", "T_2A", "T_2B", "T_3A", "T_3B")
  data_in_situ_deep$date <- as.character(data_in_situ_deep$date)
  data_copernicus$date <- as.character(data_copernicus$date)
  
  # Visualization  ####
  
  combined_data <- left_join(data_copernicus, data_in_situ_shallow, by = "date", relationship = "many-to-many",  suffix = c(".coperni", ".insit"))
  data_copernicus$date <- as.Date(data_copernicus$date)
  
  # Convert date columns to Date type if not already
  data_in_situ_deep$date <- as.Date(data_in_situ_deep$date)
  combined_data$date <- as.Date(combined_data$date)
  
  # Extract the seasons for each dataset
  # For "data_in_situ_deep"
  data_in_situ_deep$season <- ifelse(format(data_in_situ_deep$date, "%m") %in% c("01", "02", "03", "04"),
                                     "cold", 
                                     ifelse(format(data_in_situ_deep$date, "%m") %in% c("07", "08", "09", "10"), "warm", NA))
  
  # For "combined_data"
  combined_data$season <- ifelse(format(combined_data$date, "%m") %in% c("01", "02", "03", "04"),
                                 "cold", 
                                 ifelse(format(combined_data$date, "%m") %in% c("07", "08", "09", "10"), "warm", NA))
  
  # Calculate the seasonal means for both datasets (for all sites mixed)
  seasonal_means_in_situ <- data_in_situ_deep %>%
    group_by(season) %>%
    summarise(mean_temp = mean(c(T_1A, T_2A, T_3A), na.rm = TRUE))
  seasonal_means_in_situ <- seasonal_means_in_situ[-3,]
  
  
  seasonal_means_combined <- combined_data %>%
    group_by(season) %>%
    summarise(mean_temp = mean(c(RUNA1, RUNA5, RUNA9), na.rm = TRUE))
  seasonal_means_combined <- seasonal_means_combined[-3,]
  
  
  cold_season_start <- as.Date("2021-01-01")
  cold_season_end <- as.Date("2021-04-01")
  
  warm_season_start <- as.Date("2021-07-01")
  warm_season_end <- as.Date("2021-10-01")
  y_range = c(22,31)
  
  # Create the plots with the mean lines for each season
  p1_shallow_bis_2 <- ggplot(combined_data, aes(x = date)) +
    geom_line(aes(y = RUNA1, col = "RUNA1"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = RUNA5, col = "RUNA5"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = RUNA9, col = "RUNA9"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = mean.RUNA1, col = "RUNA1_in_situ"), linewidth = 1.1) +
    geom_line(aes(y = mean.RUNA9, col = "RUNA9_in_situ"), linewidth = 1.1) +
    geom_hline(data = seasonal_means_combined, aes(yintercept = mean_temp, color = season), linetype = "dashed", size = 1) +
    geom_text(data = seasonal_means_combined, aes(x = as.Date("2018-11-01"), y = mean_temp, label = paste("Mean:", round(mean_temp, 2))), 
              color = "black", size = 3, vjust = -1) +
    ylab("Mean SST(°C)") +
    xlab("Time") +
    scale_y_continuous(limits = y_range) +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") + 
    theme_minimal() +
    scale_color_manual(values = c("RUNA1" = "blue4", "RUNA5" = "green4", "RUNA9" = "orange2", "RUNA1_in_situ" = "blue4", "RUNA9_in_situ" = "orange2")) +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Create the deep plot with the mean lines for each season
  p1_deep_bis_2 <- ggplot(data_in_situ_deep, aes(x = date)) +
    geom_line(aes(y = T_1A, col = "P50A1"), linewidth = 1.1) +
    geom_line(aes(y = T_2A, col = "P50A2"), linewidth = 1.1) +
    geom_line(aes(y = T_3A, col = "P50A3"), linewidth = 1.1) +
    geom_hline(data = seasonal_means_in_situ, aes(yintercept = mean_temp, color = season), linetype = "dashed", size = 1) +
    geom_text(data = seasonal_means_in_situ, aes(x = as.Date("2021-11-01"), y = mean_temp, label = paste("Mean:", round(mean_temp, 2))), 
              color = "black", size = 3, vjust = -1) +
    ylab("Mean SST(°C)") +
    xlab("Time") +
    scale_y_continuous(limits = y_range) +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") + 
    theme_minimal() +
    scale_color_manual(values = c("P50A1" = "blue4", "P50A2" = "green4", "P50A3" = "orange2")) +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Save the combined plot
  Temp_dev_path <- here::here("outputs/Temperature_comparison_graph(4).pdf")
  cow <- cowplot::plot_grid(p1_shallow_bis_2, p1_deep_bis_2, labels = c("a", "b"), ncol = 1, nrow = 2)

  data_in_situ_deep <- data_in_situ_deep[,c(1,2,4,6)]
  data_ex_situ_shallow <- combined_data[,c(2,3,4,5)]
  
  # Import ecological data
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
  # Définir les triplicats à garder
  triplicats_to_keep <- c("P50ARMS1", "P50ARMS2", "P50ARMS3", "RUNARMS1", "RUNARMS5", "RUNARMS9")
  
  # Filtrer le tableau
  meta_filtered <- meta_mean[meta_mean$triplicat %in% triplicats_to_keep, ]
  data_filtered <- data_mean[meta_mean$triplicat %in% triplicats_to_keep, ]
  
  
  # Calculate the reponse variables
  
  
  data_in_situ_deep$season <- with(
    data_in_situ_deep,
    ifelse(
      as.integer(format(as.Date(date), "%m")) %in% 1:4,
      "warm",
      ifelse(
        as.integer(format(as.Date(date), "%m")) %in% 7:10,
        "cold",
        "inter-season"
      )
    )
  )
  
  
  data_ex_situ_shallow$season <- with(
    data_ex_situ_shallow,
    ifelse(
      as.integer(format(as.Date(date), "%m")) %in% 1:4,
      "warm",
      ifelse(
        as.integer(format(as.Date(date), "%m")) %in% 7:10,
        "cold",
        "inter-season"
      )
    )
  )
  

  temp_season_summary <- function(data, season_name) {
    
    # keep only numeric temperature columns
    temp_cols <- sapply(data, is.numeric)
    
    # subset by season
    data_season <- data[data$season == season_name, temp_cols, drop = FALSE]
    
    # compute metrics
    out <- data.frame(
      site = colnames(data_season),
      mean = sapply(data_season, mean, na.rm = TRUE),
      range = sapply(data_season, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
      max = sapply(data_season, max, na.rm = TRUE),
      sd = sapply(data_season, sd, na.rm = TRUE),
      season = season_name,
      row.names = NULL
    )
    
    return(out)
  }
  
  deep_cold <- temp_season_summary(
    data = data_in_situ_deep,
    season_name = "cold"
  )
  
  deep_warm <- temp_season_summary(
    data = data_in_situ_deep,
    season_name = "warm"
  )
  
  deep_season_stats <- rbind(deep_cold, deep_warm)
  
  shallow_cold <- temp_season_summary(
    data = data_ex_situ_shallow,
    season_name = "cold"
  )
  
  shallow_warm <- temp_season_summary(
    data = data_ex_situ_shallow,
    season_name = "warm"
  )
  
  shallow_season_stats <- rbind(shallow_cold, shallow_warm)
  
  deep_season_stats$depth <- "deep"
  shallow_season_stats$depth <- "shallow"
  
  all_temp_season_stats <- rbind(deep_season_stats, shallow_season_stats)
  
  # For deep sites (P50ARMS)
  temp_deep <- all_temp_season_stats %>%
    filter(depth == "deep") %>%
    select(site, mean, range, max, sd, season) %>%
    pivot_wider(names_from = season, values_from = c(mean, range, max, sd))  # create separate columns for warm/cold
  
  temp_deep <- temp_deep %>%
    slice(rep(1:n(), each = 3))
  
  temp_shallow <- all_temp_season_stats %>%
    filter(depth == "shallow") %>%
    select(site, mean, range, max, sd, season) %>%
    pivot_wider(names_from = season, values_from = c(mean, range, max, sd))  # create separate columns for warm/cold
  
  temp_shallow <- temp_shallow %>%
    slice(rep(1:n(), each = 3))
  
  temp_shallow$site <- gsub("RUNA1", "RUNARMS1", temp_shallow$site)
  temp_shallow$site <- gsub("RUNA5", "RUNARMS5", temp_shallow$site)
  temp_shallow$site <- gsub("RUNA9", "RUNARMS9", temp_shallow$site)
  
  
  temp_deep <- temp_deep %>%
    mutate(site = recode(site,
                         "T_1A" = "P50ARMS1",
                         "T_2A" = "P50ARMS2",
                         "T_3A" = "P50ARMS3"))
  
  temp <- as.data.frame(rbind(temp_deep, temp_shallow))
  
  rownames(temp) <- meta_filtered$arms
  temp$arms <- rownames(temp)
  
  # ACP ####
  temp_vars <- c("mean_cold",  "range_cold", "max_cold", "sd_cold","range_warm","mean_warm", "max_warm", "sd_warm")
  temp_vars_cold <- c("mean_cold",  "range_cold", "max_cold", "sd_cold")
  temp_vars_warm <- c("range_warm","mean_warm", "max_warm", "sd_warm")
  
  

  
  
  # Extraire uniquement les variables numériques
  temp_numeric_warm <-temp[, temp_vars_warm]  # supprime la première colonne "site"
  sites <- temp$site
  temp_scaled_warm <- scale(temp_numeric_warm)
  pca <- prcomp(temp_scaled_warm, center = TRUE, scale. = TRUE)
  summary(pca)
  # Biplot classique
  biplot(pca, scale = 0)
  
  # Extraire uniquement les variables numériques
  temp_numeric_cold <-temp[, temp_vars_cold]  # supprime la première colonne "site"
  sites <- temp$site
  temp_scaled_cold <- scale(temp_numeric_cold)
  pca <- prcomp(temp_scaled_cold, center = TRUE, scale. = TRUE)
  summary(pca)
  # Biplot classique
  biplot(pca, scale = 0)
  
  temp_vars_select <- c("range_warm","mean_warm")
  
  # Extraire uniquement les variables numériques
  temp_numeric_select <-temp[, temp_vars_select]  # supprime la première colonne "site"
  sites <- temp$site
  temp_scaled_select <- scale(temp_numeric_select)
  pca <- prcomp(temp_scaled_select, center = TRUE, scale. = TRUE)
  summary(pca)
  # Biplot classique
  biplot(pca, scale = 0)

  # RDA ####

  meta_temp <- merge(meta_filtered, temp, by = "arms", all.x = TRUE)
  
  # cor_matrix <- cor(meta_temp[, temp_vars], use = "pairwise.complete.obs")
  # print(cor_matrix)
  # library(corrplot)
  # corrplot(cor_matrix, method = "color", tl.cex = 0.8)
  
  # Sélectionner uniquement les colonnes de température
  
  
  # RDA
  rda_temp <- rda(data_filtered ~ ., data = meta_temp[, temp_vars_select])
  
  # Résumé
  summary(rda_temp)
  
  # Test global
  anova(rda_temp)
  
  # Test par axes
  anova(rda_temp, by = "axis")
  
  # Test par variable
  anova(rda_temp, by = "term")
  
  plot(rda_temp, scaling = 2)  # scaling 2 : distances entre échantillons
  
  diss <- vegdist(data_filtered, dist="bray")
  
  # pco
  pco1<-ade4::dudi.pco(diss, scannf = FALSE, nf = 2)
  ade4::s.label(pco1$co)
  ade4::s.label(pco1$li)
  # visualiser
  ade4::scatter(pco1)
  ade4::scatter(pco1,clab=0)
  # faire la db-rda
  dbrda <- ade4::pcaiv(dudi = pco1, df = meta_temp[, temp_vars_select], scannf = FALSE, nf = 2)

  # vif.cca(rda_temp)
  
  plot(dbrda)
  
  ## Selection automatique des variabes ####
  
  # library(vegan)
  # rda_null <- dbrda(diss ~ 1, data = meta_temp[, temp_vars])
  # rda_full <- dbrda(diss ~ mean_cold + mean_warm + range_cold + range_warm +
  #                   max_cold + max_warm + sd_cold + sd_warm,
  #                 data = meta_temp)
  # step_rda <- ordiR2step(rda_null, scope = formula(rda_full), direction = "both")
  # # summary(step_rda)
  # # 
  # ## Final formula ####
  # 
  # temp_vars_final <- meta_temp[, c("mean_warm", "max_warm", "max_cold", "sd_cold")]
  # 
  # dbrda <- ade4::pcaiv(pco1, temp_vars_final)
  # plot(dbrda)
  # 
  # ade4::s.label(dbrda$l1,1,2,clab=0)
  # ade4::s.arrow(dbrda$cor*2,1,2,clab=0.9,add.plot=TRUE)
  # 
  # dbrda2 <- capscale(data_filtered ~ mean_warm + Condition(site.x), meta_temp, dist="bray")
  # anova(dbrda2,by="terms", permu=500)
  # 
  # vif.cca(dbrda2)

  ## Utilisation des axe dela PCA pour expliquer RDA ####
  
  # PCA sur les variables de température
  temp_scaled <- scale(meta_temp[, temp_vars])
  pca_temp <- prcomp(temp_scaled)
  
  # Garde les 2-3 premiers axes qui expliquent >80% de la variance
  meta_temp$PC1 <- pca_temp$x[,1]
  meta_temp$PC2 <- pca_temp$x[,2]
  
  # RDA avec les axes de PCA comme variables explicatives
  
  dbrda3 <- capscale(data_filtered ~ PC1 + PC2, meta_temp, dist="bray")
  plot(dbrda3)
  
  # Extraire les scores des sites (coordonnées des échantillons)
  scores_sites <- scores(dbrda3, display = "sites", scaling = 2)  # scaling 2 = distances entre sites
  scores_sites <- as.data.frame(scores_sites)  # convertir en data.frame
  scores_sites$site <- rownames(scores_sites)  # ajouter le nom des sites
  

 # pour éviter le chevauchement des labels
  
  ggplot(scores_sites, aes(x = CAP1, y = CAP2, label = site)) +
    geom_point(size = 3) +
    geom_text_repel() +
    theme_minimal() +
    labs(x = "dbRDA1 (CAP1)", y = "dbRDA2 (CAP2)", title = "db-RDA - Sites")
  
  scores_vars <- scores(dbrda3, display = "bp", scaling = 2)
  arrows(0, 0, scores_vars[,1], scores_vars[,2], length = 0.1, col = "red")
  text(scores_vars[,1], scores_vars[,2], labels = rownames(scores_vars), col = "red", pos = 4)
  
  
  
  vif.cca(dbrda3)
  
  anova(dbrda3,by="terms", permu=500)
  
  vif.cca(dbrda3)
  
  # Extraire les scores des sites
  scores_sites <- scores(dbrda3, display = "sites", scaling = 2)
  scores_sites <- as.data.frame(scores_sites)
  scores_sites$site <- rownames(scores_sites)
  
  # Ajouter une colonne pour la couleur en fonction du triplicat
  scores_sites$color <- ifelse(grepl("P50ARMS", scores_sites$site), "blue3", "aquamarine3")
  
  # Extraire scores des variables explicatives
  scores_vars <- scores(dbrda3, display = "bp", scaling = 2)
  scores_vars <- as.data.frame(scores_vars)
  scores_vars$var <- rownames(scores_vars)
  
  # Tracer les points avec les bonnes couleurs
  plot(scores_sites$CAP1, scores_sites$CAP2, 
       xlab = "db-RDA1", ylab = "db-RDA2",
       pch = 19, col = scores_sites$color,
       xlim = c(min(c(scores_sites$CAP1, scores_vars$CAP1)) - 0.1,
                max(c(scores_sites$CAP1, scores_vars$CAP1)) + 0.1),
       ylim = c(min(c(scores_sites$CAP2, scores_vars$CAP2)) - 0.1,
                max(c(scores_sites$CAP2, scores_vars$CAP2)) + 0.1))
  
  # Ajouter les labels des sites
  text(scores_sites$CAP1, scores_sites$CAP2, labels = scores_sites$site, pos = 3, cex = 0.8)
  
  # Ajouter les flèches des variables explicatives
  arrows(0, 0, scores_vars$CAP1, scores_vars$CAP2, length = 0.1, col = "red")
  text(scores_vars$CAP1, scores_vars$CAP2, labels = scores_vars$var, col = "red", pos = 4)
  
  res <- anova(dbrda3,by="terms", permu=500)
  as.data.frame(res)
  
  return(res)
  
  }
  