#' Temperature 
#'
#' @param temp_data_1A path to the T° 
#' @param temp_data_1B
#' @param temp_data_2A
#' @param temp_data_2B
#' @param temp_data_3A
#' @param temp_data_3B
#'
#' @return ?
#' @export
#'
fun_temp <- function(temp_1A,
                     temp_1B,
                     temp_2A,
                     temp_2B,
                     temp_3A,
                     temp_3B){
  # temp_1A = targets::tar_read("temp_data_1A")
  # temp_1B = targets::tar_read("temp_data_1B")
  # temp_2A = targets::tar_read("temp_data_2A")
  # temp_2B = targets::tar_read("temp_data_2B")
  # temp_3A = targets::tar_read("temp_data_3A")
  # temp_3B = targets::tar_read("temp_data_3B")
  
  library(zoo)
  
  data_temp_1A <- read.csv(temp_1A, row.names = 1, sep = ";", dec = ",")
  data_temp_1B <- read.csv(temp_1B, row.names = 1, sep = ";", dec = ",")
  data_temp_2A <- read.csv(temp_2A, row.names = 1, sep = ";", dec = ",")
  data_temp_2B <- read.csv(temp_2B, row.names = 1, sep = ";", dec = ",")
  data_temp_3A <- read.csv(temp_3A, row.names = 1, sep = ";", dec = ",")
  data_temp_3B <- read.csv(temp_3B, row.names = 1, sep = ";", dec = ",")

  
  # data_temp_1A$Date=strptime(as.Date(data_temp_1A$Date, tz = "UTC"), format="%Y-%m-%d")
  
  P50ARMS=list(data_temp_1A,data_temp_1B,data_temp_2A,data_temp_2B,data_temp_3A,data_temp_3B)
  library(dplyr)
  P50ARMS <- lapply(P50ARMS, function(new) {
    new <- new %>%
      filter(Date >= as.Date("2021-11-26") & Date <= as.Date("2023-11-25"))
    return(new)
  })
  
  
  for (i in 1:length(P50ARMS)) {
    tableau <- P50ARMS[[i]]
    premiere_date <- min(tableau$Date)
    derniere_date <- max(tableau$Date)
    nrow <- nrow(tableau)
    cat("Tableau", i, ": Première date =", premiere_date, ", Dernière date =", derniere_date, ", nrow = ", nrow, "\n")
  }


  data_temp <- data.frame(Date = P50ARMS[[1]]$Date,
                          T_1A = P50ARMS[[1]]$Temp,
                          T_1B = P50ARMS[[2]]$Temp,
                          T_2A = P50ARMS[[3]]$Temp,
                          T_2B = P50ARMS[[4]]$Temp,
                          T_3A = P50ARMS[[5]]$Temp,
                          T_3B = P50ARMS[[6]]$Temp)
  
  # Version longue ####
  # 
  # T_1A_mean <- aggregate(T_1A ~ Date , FUN = mean, data = data_temp)
  # T_1B_mean <- aggregate(T_1B ~ Date , FUN = mean, data = data_temp)
  # T_2A_mean <- aggregate(T_2A ~ Date , FUN = mean, data = data_temp)
  # T_2B_mean <- aggregate(T_2B ~ Date , FUN = mean, data = data_temp)
  # T_3A_mean <- aggregate(T_3A ~ Date , FUN = mean, data = data_temp)
  # T_3B_mean <- aggregate(T_3B ~ Date , FUN = mean, data = data_temp)
  # 
  # 
  # data_mean_temp <- data.frame(T_1A = T_1A_mean$T_1A,
  #                              T_1B = T_1B_mean$T_1B,
  #                              T_2A = T_2A_mean$T_2A,
  #                              T_2B = T_2B_mean$T_2B,
  #                              T_3A = T_3A_mean$T_3A,
  #                              T_3B = T_3B_mean$T_3B)
  # 
  # rownames(data_mean_temp) <- T_1A_mean$Date
  # 
  # data_mean_temp$T_1A <- zoo::rollmean(as.numeric(data_mean_temp$T_1A), k = 168, align = "left", fill = NA)
  # data_mean_temp$T_1B <- zoo::rollmean(as.numeric(data_mean_temp$T_1B), k = 168, align = "left", fill = NA)
  # data_mean_temp$T_2A <- zoo::rollmean(as.numeric(data_mean_temp$T_2A), k = 168, align = "left", fill = NA)
  # data_mean_temp$T_2B <- zoo::rollmean(as.numeric(data_mean_temp$T_2B), k = 168, align = "left", fill = NA)
  # data_mean_temp$T_3A <- zoo::rollmean(as.numeric(data_mean_temp$T_3A), k = 168, align = "left", fill = NA)
  # data_mean_temp$T_3B <- zoo::rollmean(as.numeric(data_mean_temp$T_3B), k = 168, align = "left", fill = NA)
  
  # Version courte ####
  # Liste des colonnes à traiter
  columns_to_process <- c("T_1A", "T_1B", "T_2A", "T_2B", "T_3A", "T_3B")
  
  # Initialisation du dataframe résultant
  data_mean_temp <- data.frame(Date = unique(data_temp$Date))
  
  # Boucle pour calculer la moyenne par colonne
  for (col in columns_to_process) {
    col_mean <- aggregate(data_temp[[col]] ~ data_temp$Date, FUN = mean)
    col_mean <- zoo::rollmean(as.numeric(col_mean[, 2]), k = 7, align = "center", fill = NA)
    data_mean_temp[[col]] <- col_mean[match(data_mean_temp$Date, unique(data_temp$Date))]
  }
  
  # Affectation des noms de lignes
  rownames(data_mean_temp) <- data_mean_temp$Date
  
  # Suppression des lignes avec des valeurs manquantes (NA)
  data_mean_temp <- na.omit(data_mean_temp)
  
  # Suppression de la colonne "Date" redondante
  data_mean_temp <- data_mean_temp[, -1]
  
  library(ggplot2)
  
  library(tidyr)
  
  
  # Convert Date column to Date class
  data_mean_temp$Date <- as.Date(rownames(data_mean_temp))
  
  # Your ggplot code with legend
  p1 <- ggplot(data_mean_temp, aes(x = Date)) +
    geom_line(aes(y = T_1A, col = "T_1A"), linewidth = 1.1) +
    geom_line(aes(y = T_1B, col = "T_1B"), linewidth = 1.1) +
    geom_line(aes(y = T_2A, col = "T_2A"), linewidth = 1.1) +
    geom_line(aes(y = T_2B, col = "T_2B"), linewidth = 1.1) +
    geom_line(aes(y = T_3A, col = "T_3A"), linewidth = 1.1) +
    geom_line(aes(y = T_3B, col = "T_3B"), linewidth = 1.1) +
    ylab("Temperature mean (°C)") +
    xlab("Time") +
    scale_color_manual(values = c(T_1A = "blue4", T_1B = "blue4", T_2A = "green4", T_2B = "green4", T_3A = "coral", T_3B = "coral")) +
    theme(legend.position = "top", legend.title = element_blank())
  
  Temp_roolmean_path <- here::here("outputs/Temperature_rollmean.pdf")
  
  ggsave(filename =  Temp_roolmean_path, plot = p1 , width = 12, height = 8)
 
  data_mean_temp <- data_mean_temp[,-ncol(data_mean_temp)]
  
  # Ecart à la moyenne ####
  
  # Faire une courbe d'écart à la moyenne
  
  
  data_mean_temp$mean_all_site <- rowMeans(data_mean_temp)
  
  dev_1A <- data_mean_temp$T_1A-data_mean_temp$mean_all_site
  dev_1B <- data_mean_temp$T_1B-data_mean_temp$mean_all_site
  dev_2A <- data_mean_temp$T_2A-data_mean_temp$mean_all_site
  dev_2B <- data_mean_temp$T_2B-data_mean_temp$mean_all_site
  dev_3A <- data_mean_temp$T_3A-data_mean_temp$mean_all_site
  dev_3B <- data_mean_temp$T_3B-data_mean_temp$mean_all_site
  
  temp_deviation <- data.frame(Date = as.Date(rownames(data_mean_temp)),
                               D_1A = dev_1A,
                               D_1B = dev_1B,
                               D_2A = dev_2A,
                               D_2B = dev_2B,
                               D_3A = dev_3A,
                               D_3B = dev_3B)
  
  p2 <- ggplot(temp_deviation, aes(x = Date)) +
    geom_line(aes(y = D_1A, col = "D_1A"), linewidth = 1.1) +
    geom_line(aes(y = D_1B, col = "D_1B"), linewidth = 1.1) +
    geom_line(aes(y = D_2A, col = "D_2A"), linewidth = 1.1) +
    geom_line(aes(y = D_2B, col = "D_2B"), linewidth = 1.1) +
    geom_line(aes(y = D_3A, col = "D_3A"), linewidth = 1.1) +
    geom_line(aes(y = D_3B, col = "D_3B"), linewidth = 1.1) +
    ylab("Temperature deviation (°C)") +
    xlab("Time") +
    scale_color_manual(values = c(D_1A = "blue4", D_1B = "blue4", D_2A = "green4", D_2B = "green4", D_3A = "coral", D_3B = "coral")) +
    theme(legend.position = "top", legend.title = element_blank())
  
  Temp_dev_path <- here::here("outputs/Temperature_deviation.pdf")
  
  ggsave(filename =  Temp_dev_path, plot = p2 , width = 12, height = 8)
  
  
  # install.packages("e1071")
  library(e1071)
  kurtosis_table <- data.frame(Site = character(0), Kurtosis = numeric(0))
  for (col in names(data_mean_temp)) {
    kurt <- kurtosis(data_mean_temp[[col]])
    kurtosis_table <- rbind(kurtosis_table, data.frame(Site = col, Kurtosis = kurt))
  }
  
  ggplot(data_mean_temp, aes(x=T_3A)) + 
    geom_density()
  
  return(Temp_roolmean_path)
}
