#' Temperature comparison
#' 
#'
#' @return the path to the subseted raw data file
#' @export
#'
fun_temperature_comparison <- function(temp_extraction, temp_in_situ){
  # temp_extraction = targets::tar_read("temp_extraction")
  # temp_in_situ = targets::tar_read("temp_in_situ")

  # Comparaison shallow 2019-2020 / Deep 2022-2023 #### 
  
  library(ggplot2) 
  data_copernicus <- read.csv(temp_extraction[1])
  data_copernicus$date <- as.Date(data_copernicus$date)
  
  
  data_in_situ_shallow <- read.csv(temp_in_situ[grepl("shallow",temp_in_situ)])
  data_in_situ_deep <- read.csv(temp_in_situ[grepl("deep",temp_in_situ)])
  colnames(data_in_situ_deep) <- c("date", "T_1A", "T_1B", "T_2A", "T_2B", "T_3A", "T_3B")
  data_in_situ_deep$date <- as.character(data_in_situ_deep$date)
  data_copernicus$date <- as.character(data_copernicus$date)
  
  library(dplyr)
  
  combined_data <- left_join(data_copernicus, data_in_situ_shallow, by = "date", relationship = "many-to-many",  suffix = c(".coperni", ".insit"))
  data_copernicus$date <- as.Date(data_copernicus$date)
  p1_shallow <- ggplot(data_copernicus, aes(x = date)) +
    geom_line(aes(y = RUNA1, col = "RUNA1"), linewidth = 1.1) +
    geom_line(aes(y = RUNA5, col = "RUNA5"), linewidth = 1.1) +
    geom_line(aes(y = RUNA9, col = "RUNA9"), linewidth = 1.1) +
    ylab("Temperature mean (°C)") +
    xlab("Time") +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") + # Affiche mois/année
    theme_minimal() +
    scale_color_manual(values = c(RUNA1 = "blue",  RUNA5 = "green", RUNA9 = "coral")) +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1) # Écrit les dates de biais
    )
  
  
  y_range <- range(c(
    combined_data$RUNA1, combined_data$RUNA5, combined_data$RUNA9, 
    combined_data$mean.RUNA1, combined_data$mean.RUNA2,
    data_in_situ_deep$T_1A, data_in_situ_deep$T_2A, data_in_situ_deep$T_3A
  ), na.rm = TRUE)
  
  combined_data$date <- as.Date(combined_data$date)
  y_range = c(22,31)
  
  p1_shallow_bis <- ggplot(combined_data, aes(x = date)) +
    geom_line(aes(y = RUNA1, col = "RUNA1"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = RUNA5, col = "RUNA5"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = RUNA9, col = "RUNA9"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = mean.RUNA1, col = "RUNA9_in_situ"), linewidth = 1.1) +
    geom_line(aes(y = mean.RUNA2, col = "RUNA1_in_situ"), linewidth = 1.1) +
    ylab("Temperature mean (°C)") +
    xlab("Time") +
    scale_y_continuous(limits = y_range) +  # Set the Y-axis limits
    scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") + # Affiche mois/année
    theme_minimal() +
    scale_color_manual(values = c(RUNA1 = "blue",  RUNA5 = "green", RUNA9 = "coral",RUNA1_in_situ = "blue",RUNA9_in_situ = "coral")) +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1) # Écrit les dates de biais
    )
  
  data_in_situ_deep$date <- as.Date(data_in_situ_deep$date)
  
  p1_deep <- ggplot(data_in_situ_deep, aes(x = date)) +
    geom_line(aes(y = T_1A, col = "P50A1"), linewidth = 1.1) +
    geom_line(aes(y = T_2A, col = "P50A2"), linewidth = 1.1) +
    geom_line(aes(y = T_3A, col = "P50A3"), linewidth = 1.1) +
    ylab("Temperature mean (°C)") +
    xlab("Time") +
    scale_y_continuous(limits = y_range) +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") + # Affiche mois/année
    theme_minimal() +
    scale_color_manual(values = c(P50A1 = "blue4",  P50A2 = "green4", P50A3 = "orange")) +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1) # Écrit les dates de biais
    )
  
  Temp_dev_path <- here::here("outputs/Temperature_comparison_graph(3).pdf")
  
  cow <- cowplot::plot_grid(p1_shallow_bis,
                            p1_deep,
                            labels = c("a", "b"),
                            ncol = 1, 
                            nrow = 2)
  
  ggsave(filename =  Temp_dev_path, plot = cow , width = 13, height = 7)
  
  ##############################################################
  library(dplyr)
  library(ggplot2)
  
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
    ylab("Temperature mean (°C)") +
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
    ylab("Temperature mean (°C)") +
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
  ggsave(Temp_dev_path, cow, width = 13, height = 7)
  
  
  ################################################################
  
  
  data_in_situ_deep <- data_in_situ_deep[,c(1,2,4,6)]
  data_ex_situ_shallow <- combined_data[,c(2,3,4,5)]
  
  
  filter_season <- function(data, date_column, season) {
    data <- data %>%
      mutate(date = as.Date(!!sym(date_column))) %>%
      filter(if (season == "cold") {
        format(date, "%m-%d") >= "07-01" & format(date, "%m-%d") < "10-01"
      } else if (season == "hot") {
        format(date, "%m-%d") >= "01-01" & format(date, "%m-%d") < "04-01"
      } else {
        FALSE
      })
    return(data)
  }
  
  # Filter datasets by season
  hot_deep <- filter_season(data_in_situ_deep, "date", "hot")
  cold_deep <- filter_season(data_in_situ_deep, "date", "cold")
  
  hot_shallow <- filter_season(data_ex_situ_shallow, "date", "hot")
  cold_shallow <- filter_season(data_ex_situ_shallow, "date", "cold")
  
  deep_T <- data.frame(hot_deep_1 = c(mean(hot_deep$T_1A), sd(hot_deep$T_1A)),
                       hot_deep_2 = c(mean(hot_deep$T_2A), sd(hot_deep$T_2A)),
                       hot_deep_3 = c(mean(hot_deep$T_3A), sd(hot_deep$T_3A)),
                       cold_deep_1 = c(mean(cold_deep$T_1A), sd(cold_deep$T_1A)),
                       cold_deep_2 = c(mean(cold_deep$T_2A), sd(cold_deep$T_2A)),
                       cold_deep_3 = c(mean(cold_deep$T_3A), sd(cold_deep$T_3A))
  )
  
  rownames(deep_T) <- c("mean", "sd")
  
  Temp_summary_Hot <- data.frame(Cap_La_Houssaye = c(paste0(round(mean(hot_shallow$RUNA1), 2)," ± ", round(sd(hot_shallow$RUNA1), 2)),
                                                 paste0(round(mean(hot_deep$T_1A), 2)," ± ", round(sd(hot_deep$T_1A), 2))),
                             Saint_Leu = c(paste0(round(mean(hot_shallow$RUNA5), 2)," ± ", round(sd(hot_shallow$RUNA5), 2)),
                                           paste0(round(mean(hot_deep$T_2A), 2)," ± ", round(sd(hot_deep$T_2A), 2))),
                             Grand_Bois = c(paste0(round(mean(hot_shallow$RUNA9), 2)," ± ", round(sd(hot_shallow$RUNA9), 2)),
                                            paste0(round(mean(hot_deep$T_3A), 2)," ± ", round(sd(hot_deep$T_3A), 2))))
                          
  
  rownames(Temp_summary_Hot) <- c("Shallow", "Deep")
  
  Temp_summary_Cold <- data.frame(Cap_La_Houssaye = c(paste0(round(mean(cold_shallow$RUNA1), 2)," ± ", round(sd(cold_shallow$RUNA1), 2)),
                                                     paste0(round(mean(cold_deep$T_1A), 2)," ± ", round(sd(cold_deep$T_1A), 2))),
                                 Saint_Leu = c(paste0(round(mean(cold_shallow$RUNA5), 2)," ± ", round(sd(cold_shallow$RUNA5), 2)),
                                               paste0(round(mean(cold_deep$T_2A), 2)," ± ", round(sd(cold_deep$T_2A), 2))),
                                 Grand_Bois = c(paste0(round(mean(cold_shallow$RUNA9), 2)," ± ", round(sd(cold_shallow$RUNA9), 2)),
                                                paste0(round(mean(cold_deep$T_3A), 2)," ± ", round(sd(cold_deep$T_3A), 2))))
  
  
  rownames(Temp_summary_Cold) <- c("Shallow", "Deep")
  
  hist(hot_shallow$RUNA5)
    
  shapiro.test(hot_deep$T_1A)
  
  
  # Replace t-test with Wilcoxon test
  wilcox_test_hot_RUNA1 <- wilcox.test(hot_deep$T_1A, hot_shallow$RUNA1, na.rm = TRUE)
  wilcox_test_cold_RUNA1 <- wilcox.test(cold_deep$T_1A, cold_shallow$RUNA1, na.rm = TRUE)
  wilcox_test_hot_RUNA5 <- wilcox.test(hot_deep$T_2A, hot_shallow$RUNA5, na.rm = TRUE)
  wilcox_test_cold_RUNA5 <- wilcox.test(cold_deep$T_2A, cold_shallow$RUNA5, na.rm = TRUE)
  wilcox_test_hot_RUNA9 <- wilcox.test(hot_deep$T_3A, hot_shallow$RUNA9, na.rm = TRUE)
  wilcox_test_cold_RUNA9 <- wilcox.test(cold_deep$T_3A, cold_shallow$RUNA9, na.rm = TRUE)
  
  # Comparaison shallow 2022-2023 / Deep 2022-2023 #### 
  
  data_copernicus <- read.csv(temp_extraction[1])
  data_copernicus$date <- as.Date(data_copernicus$date)
  
  
  library(ggplot2)
  data_copernicus <- read.csv(temp_extraction[2])
  data_copernicus$date <- as.Date(data_copernicus$date)
  
  data_in_situ_shallow <- read.csv(temp_in_situ[grepl("shallow",temp_in_situ)])
  data_in_situ_deep <- read.csv(temp_in_situ[grepl("deep",temp_in_situ)])
  colnames(data_in_situ_deep) <- c("date", "T_1A", "T_1B", "T_2A", "T_2B", "T_3A", "T_3B")
  data_in_situ_deep$date <- as.character(data_in_situ_deep$date)
  data_copernicus$date <- as.character(data_copernicus$date)
  
  library(dplyr)
  
  combined_data <- left_join(data_copernicus, data_in_situ_deep, by = "date", relationship = "one-to-one",  suffix = c(".coperni", ".insit"))
  data_copernicus$date <- as.Date(data_copernicus$date)
  
  combined_data$date <- as.Date(combined_data$date)
  
  
  # Extract the seasons for each dataset
  
  # For "combined_data"
  combined_data$season <- ifelse(format(combined_data$date, "%m") %in% c("01", "02", "03", "04"),
                                 "cold", 
                                 ifelse(format(combined_data$date, "%m") %in% c("07", "08", "09", "10"), "warm", NA))
  
  # Calculate the seasonal means for both datasets (for all sites mixed)
 
  # seasonal_means_combined <- combined_data %>%
  #   group_by(season) %>%
  #   summarise(mean_temp = mean(c(RUNA1, RUNA5, RUNA9), na.rm = TRUE))
  # seasonal_means_combined <- seasonal_means_combined[-3,]
  # 
  # seasonal_sd_combined <- combined_data %>%
  #   group_by(season) %>%
  #   summarise(sd_temp = sd(c(RUNA1, RUNA5, RUNA9), na.rm = TRUE))
  # seasonal_sd_combined <- seasonal_sd_combined[-3,]
  # 
  # seasonal_means_combined <- left_join(seasonal_means_combined,seasonal_sd_combined, by = "season" )
  
  p1_deep_bis_3 <- ggplot(combined_data, aes(x = date)) +
    geom_line(aes(y = RUNA1, col = "RUNA1"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = RUNA5, col = "RUNA5"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = RUNA9, col = "RUNA9"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = T_1A, col = "P50A1"), linewidth = 1.1) +
    geom_line(aes(y = T_2A, col = "P50A2"), linewidth = 1.1) +
    geom_line(aes(y = T_3A, col = "P50A3"), linewidth = 1.1) +
    # geom_hline(data = seasonal_means_combined, aes(yintercept = mean_temp, color = season), linetype = "dashed", size = 1) +
    # geom_text(data = seasonal_means_combined, aes(x = as.Date("2021-11-01"), y = mean_temp, label = paste("Mean:", round(mean_temp, 2),"±",round(sd_temp, 2))), 
    #           color = "black", size = 3, vjust = -1) +
    ylab("Temperature mean (°C)") +
    xlab("Time") +
    scale_y_continuous(limits = y_range) +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") + # Affiche mois/année
    theme_minimal() +
    scale_color_manual(values = c(RUNA1 = "blue",  RUNA5 = "green", RUNA9 = "coral", "P50A1" = "blue4", "P50A2" = "green4", "P50A3" = "orange2")) +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1) # Écrit les dates de biais
    )
  
  combined_data$date <- as.Date(combined_data$date)
  
  data_in_situ_deep$date <- as.Date(data_in_situ_deep$date)

  Temp_dev_path <- here::here("outputs/Temperature_comparison_graph(5).pdf")
  
  cow <- cowplot::plot_grid(p1_shallow_bis,
                            p1_deep_bis_3,
                            labels = c("a", "b"),
                            ncol = 1, 
                            nrow = 2)
  
  ggsave(Temp_dev_path, cow, width = 13, height = 7)
  
  combined_data$season <- ifelse(format(combined_data$date, "%m") %in% c("01", "02", "03", "04"),
                                     "hot", 
                                     ifelse(format(combined_data$date, "%m") %in% c("07", "08", "09", "10"), "cool", NA))
 
  
  filter_season <- function(data, date_column, season) {
    # Convert the specified column to Date
    data <- data %>%
      mutate(date = as.Date(.data[[date_column]]))
    
    # Define season-specific filtering logic
    if (season == "cold") {
      filtered_data <- data %>%
        filter(format(date, "%m-%d") >= "07-01" & format(date, "%m-%d") < "10-01")
    } else if (season == "hot") {
      filtered_data <- data %>%
        filter(format(date, "%m-%d") >= "01-01" & format(date, "%m-%d") < "04-01")
    } else {
      stop("Invalid season. Please choose 'cold' or 'hot'.")
    }
    
    return(filtered_data)
  }
  
  # combined_data$date <- as.Date(combined_data$date)
  # Filter datasets by season
  hot <- filter_season(combined_data, "date", "hot")
  cold <- filter_season(combined_data, "date", "cold")
  

  rownames(deep_T) <- c("mean", "sd")
  
  Temp_summary_Hot <- data.frame(Cap_La_Houssaye = c(paste0(round(mean(hot$RUNA1), 2)," ± ", round(sd(hot$RUNA1), 2)),
                                                     paste0(round(mean(hot$T_1A), 2)," ± ", round(sd(hot$T_1A), 2))),
                                 Saint_Leu = c(paste0(round(mean(hot$RUNA5), 2)," ± ", round(sd(hot$RUNA5), 2)),
                                               paste0(round(mean(hot$T_2A), 2)," ± ", round(sd(hot$T_2A), 2))),
                                 Grand_Bois = c(paste0(round(mean(hot$RUNA9), 2)," ± ", round(sd(hot$RUNA9), 2)),
                                                paste0(round(mean(hot$T_3A), 2)," ± ", round(sd(hot$T_3A), 2))))
  
  
  rownames(Temp_summary_Hot) <- c("Shallow", "Deep")
  
  Temp_summary_Cold <- data.frame(Cap_La_Houssaye = c(paste0(round(mean(cold$RUNA1), 2)," ± ", round(sd(cold$RUNA1), 2)),
                                                      paste0(round(mean(cold$T_1A), 2)," ± ", round(sd(cold$T_1A), 2))),
                                  Saint_Leu = c(paste0(round(mean(cold$RUNA5), 2)," ± ", round(sd(cold$RUNA5), 2)),
                                                paste0(round(mean(cold$T_2A), 2)," ± ", round(sd(cold$T_2A), 2))),
                                  Grand_Bois = c(paste0(round(mean(cold$RUNA9), 2)," ± ", round(sd(cold$RUNA9), 2)),
                                                 paste0(round(mean(cold$T_3A), 2)," ± ", round(sd(cold$T_3A), 2))))
  
  
  rownames(Temp_summary_Cold) <- c("Shallow", "Deep")
  
  hist(hot$RUNA5)
  
  shapiro.test(hot$T_1A)
  
  # Replace t-test with Wilcoxon test
  wilcox_test_hot_RUNA1 <- wilcox.test(hot$T_1A, hot$RUNA1, na.rm = TRUE)
  wilcox_test_cold_RUNA1 <- wilcox.test(cold$T_1A, cold$RUNA1, na.rm = TRUE)
  wilcox_test_hot_RUNA5 <- wilcox.test(hot$T_2A, hot$RUNA5, na.rm = TRUE)
  wilcox_test_cold_RUNA5 <- wilcox.test(cold$T_2A, cold$RUNA5, na.rm = TRUE)
  wilcox_test_hot_RUNA9 <- wilcox.test(hot$T_3A, hot$RUNA9, na.rm = TRUE)
  wilcox_test_cold_RUNA9 <- wilcox.test(cold$T_3A, cold$RUNA9, na.rm = TRUE)
  
  
  ## Graph pour montrer toutes les infos dispo ####
  data_rodrigues_22_23 <- read.csv(temp_extraction["data_copernicus_rodrigues_2022_2023"])
  data_rodrigues_22_23$date <- as.Date(data_rodrigues_22_23$date)
  colnames(data_rodrigues_22_23) <- c("X", "date", "RODA1", "RODA2","RODA3")
  
  data_rodrigues_16_17 <- read.csv(temp_extraction["data_copernicus_rodrigues_2016_2017"])
  data_rodrigues_16_17$date <- as.Date(data_rodrigues_16_17$date)
  colnames(data_rodrigues_16_17) <- c("X", "date", "RODA1", "RODA2","RODA3")
  
  combined_data <- left_join(data_copernicus, data_in_situ_deep,by = "date", relationship = "one-to-one",  suffix = c(".coperni", ".insit"))
  data_copernicus$date <- as.Date(data_copernicus$date)
  
  combined_data_roda <- left_join(combined_data, data_rodrigues_22_23,by = "date", relationship = "one-to-one",  suffix = c(".reunion", ".rodri"))
  combined_data_roda$date <- as.Date(combined_data_roda$date)
  
  library(dplyr)
  library(ggplot2)
  
  # Convert date columns to Date type if not already
  data_in_situ_deep$date <- as.Date(data_in_situ_deep$date)
  combined_data_roda$date <- as.Date(combined_data_roda$date)
  
  y_range <- range(c(
    combined_data_roda$RUNA1, 
    combined_data_roda$RUNA5, 
    combined_data_roda$RUNA9,  
    combined_data_roda$T_1A, 
    combined_data_roda$T_2A, 
    combined_data_roda$T_3A, 
    combined_data_roda$RODA1, 
    combined_data_roda$RODA2, 
    combined_data_roda$RODA3,
    data_rodrigues_16_17$RODA1,
    data_rodrigues_16_17$RODA2,
    data_rodrigues_16_17$RODA3
  ), na.rm = TRUE)
  
  y_range = c(22,31)
  
  # Create the plots with the mean lines for each season
    p2 <- ggplot(combined_data_roda, aes(x = date)) +
    geom_line(aes(y = RUNA1, col = "RUNA1"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = RUNA5, col = "RUNA5"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = RUNA9, col = "RUNA9"), linewidth = 0.8, linetype = 2) +
    geom_line(aes(y = T_1A, col = "P50A1"), linewidth = 1.1) +
    geom_line(aes(y = T_2A, col = "P50A2"), linewidth = 1.1) +
    geom_line(aes(y = T_3A, col = "P50A3"), linewidth = 1.1) +
    geom_line(aes(y = RODA1, col = "RODA1"), linewidth = 1.1, linetype = 2) +
    geom_line(aes(y = RODA2, col = "RODA2"), linewidth = 1.1, linetype = 2) +
    geom_line(aes(y = RODA3, col = "RODA3"), linewidth = 1.1, linetype = 2) +
    # geom_hline(data = seasonal_means_combined, aes(yintercept = mean_temp, color = season), linetype = "dashed", size = 1) +
    # geom_text(data = seasonal_means_combined, aes(x = as.Date("2021-11-01"), y = mean_temp, label = paste("Mean:", round(mean_temp, 2),"±",round(sd_temp, 2))), 
    #           color = "black", size = 3, vjust = -1) +
    ylab("Temperature mean (°C)") +
    xlab("Time") +
    scale_y_continuous(limits = y_range) +
    scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") + # Affiche mois/année
    theme_minimal() +
    scale_color_manual(values = c(RUNA1 = "blue",  RUNA5 = "green", RUNA9 = "coral", "P50A1" = "blue4", "P50A2" = "green4", "P50A3" = "orange2", "RODA1" =  "#BC7AF9", "RODA2" = "#8B48BF", "RODA3" = "#5A189A")) +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1) # Écrit les dates de biais
    )
    
    p3 <- ggplot(data_rodrigues_16_17, aes(x = date)) +
      geom_line(aes(y = RODA1, col = "RODA1"), linewidth = 1.1, linetype = 2) +
      geom_line(aes(y = RODA2, col = "RODA2"), linewidth = 1.1, linetype = 2) +
      geom_line(aes(y = RODA3, col = "RODA3"), linewidth = 1.1, linetype = 2) +
      # geom_hline(data = seasonal_means_combined, aes(yintercept = mean_temp, color = season), linetype = "dashed", size = 1) +
      # geom_text(data = seasonal_means_combined, aes(x = as.Date("2021-11-01"), y = mean_temp, label = paste("Mean:", round(mean_temp, 2),"±",round(sd_temp, 2))), 
      #           color = "black", size = 3, vjust = -1) +
      ylab("Temperature mean (°C)") +
      xlab("Time") +
      scale_y_continuous(limits = y_range) +
      scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month") + # Affiche mois/année
      theme_minimal() +
      scale_color_manual(values = c("RODA1" = "#BC7AF9", "RODA2" = "#8B48BF", "RODA3" = "#5A189A")) +
      theme(
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1) # Écrit les dates de biais
      )
  
    cow <- cowplot::plot_grid(p3,
                              p1_shallow_bis,
                              p2,
                              labels = c("a", "b","c"),
                              ncol = 1, 
                              nrow = 3)
    
    Temp_dev_path <- here::here("outputs/Temperature_comparison_graph(6).pdf")
    
    ggsave(Temp_dev_path, cow, width = 14, height = 11)
  

  
  return(NULL) 
}