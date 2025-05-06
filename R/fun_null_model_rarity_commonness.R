#' Rarity and commonness of morphospecies
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
null_model_rarity_commonness <- function(data_and_meta_clean){
  
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)
  # data <- subset(data, meta$island == "Reunion")
  # meta <- subset(meta, meta$island == "Reunion")
  # 
  
  data <- data[, colSums(data) != 0]
  
  msp_list <- names(data) 
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
                                               "Erect_Rhodophyta_algae",
                                               "Cni_Plumulariidae")]

  
  data <- data[, msp_list_filter]
  data_pa <- vegan::decostand(data, "pa")
  
  
  # Quelle forme a ma distribution des frequence d'occurence à large echelle ####
  f <- colSums(data_pa)/nrow(data_pa)
  hist(f, breaks = 50)
  
  f_log <- colSums(data_pa)/nrow(data_pa)
  
  f_log <- f_log[f_log != 0]
  
  df <- data.frame(Frequency = f_log)
  
  # Plot histogram
  aa <- ggplot(df, aes(x = Frequency)) +
    geom_histogram(binwidth = 0.002, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(title = "Histogram of Species Occurrence Frequencies",
         x = "Occurrence Frequency",
         y = "Number of Species") +
    theme_minimal()
  
  # Quelle forme a ma distribution des frequence d'occurence à l'échelle d'un site ? ####
  
  data_RUNA8 <- subset(data_pa, meta$triplicat == "RUNARMS8")
  
  f_RUNA8 <- colSums(data_RUNA8)/nrow(data_RUNA8)
  
  f_RUNA8 <- f_RUNA8[f_RUNA8 != 0]
  
  df <- data.frame(Frequency = f_RUNA8  )
  
  # Plot histogram
  bb <- ggplot(df, aes(x = Frequency)) +
    geom_histogram(binwidth = 0.002, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(title = "Histogram of Species Occurrence Frequencies",
         x = "Occurrence Frequency",
         y = "Number of Species") +
    theme_minimal()
  
  # Compute a null model ####
  
  # Number of simulations
  n_simulations <- 1000  # Adjust as needed
  community_size <- 16 * 3  # Each simulated community has 48 sites
  
  # Define frequency intervals
  intervals <- seq(0, 1, by = 0.002)
  interval_labels <- paste0("[", head(intervals, -1), "-", tail(intervals, -1), "]")
  
  # Initialize a table to count species in each interval
  interval_counts <- setNames(rep(0, length(interval_labels)), interval_labels)
  
  # Loop over simulations
  for (i in 1:n_simulations) {
    # Generate a simulated presence-absence matrix for one community
    sim_matrix <- t(replicate(community_size, rbinom(length(f), 1, f)))
    
    # Associate species names
    colnames(sim_matrix) <- names(f)
    
    # community composition simulated
    community_composition_sim <- colSums(sim_matrix)
    
  
    # P/a transform
    community_composition_sim <- ifelse(community_composition_sim >= 1, 1, 0)
    
    f_sim <- ifelse(community_composition_sim >= 1, f, 0)
    
    # Remove species with f_sim = 0
    f_sim <- f_sim[f_sim > 0]
    
    # Bin frequencies into intervals of 0.02
    bin_counts <- table(cut(f_sim, breaks = intervals, labels = interval_labels, include.lowest = TRUE))
    
    # Accumulate counts
    interval_counts <- interval_counts + bin_counts
  }
  
  # Convert results to a data frame
  null_distribution_df <- data.frame(
    Interval = names(interval_counts),
    Count = (as.numeric(interval_counts)/n_simulations)
  )
  
  # Print results
  print(null_distribution_df)
  
  # Optional: Plot histogram of species counts per interval

  cc <- ggplot(null_distribution_df, aes(x = Interval, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Distribution of Species Frequencies", x = "Frequency Interval", y = "Number of Species")
  
  (sum(null_distribution_df$Count*1000)/4)/1000
  
  null_distribution_df
  # plot a grid of SAD for each unique site ####
  
  
  # Extract unique sites
  sites <- unique(meta$triplicat)
  
  # Initialize an empty list to store data frames
  df_list <- list()
  
  site <- sites[1]
  
  # Loop through each site
  for (site in sites) {
    # Subset data for the current site
    data_site <- subset(data_pa, meta$triplicat == site)
    
    # Compute species occurrence frequency
    f_site <- colSums(data_site) / nrow(data_site)
    
    f_site <- ifelse(f_site > 0, 1, 0)
    
    f_global_dans_site <- ifelse(f_site >= 1, f, 0)
    
    
    # Create a data frame and add a column for the site name
    df_site <- data.frame(Frequency = f_global_dans_site, Site = site)
    
    # Store in the list
    df_list[[site]] <- df_site
  }
  
  # Combine all site data into one data frame
  df_all <- bind_rows(df_list)
  
  df_all <- df_all[df_all$Frequency != 0,]
  
  # Plot histograms for all sites using facet_wrap()
  # ggplot(df_all, aes(x = Frequency)) +
  #   geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black", alpha = 0.7) +
  #   facet_wrap(~ Site, scales = "free_y") +  # Create a grid of plots
  #   labs(title = "Histogram of Species Occurrence Frequencies per Site",
  #        x = "Occurrence Frequency",
  #        y = "Number of Species") +
  #   theme_minimal() +
  #   ylim(0, 20)
  # 
  
  # Define bin width
  bin_width <- 0.002
  bins <- seq(0, 1, by = bin_width)  # Define bin edges
  
  # Create custom interval labels
  interval_labels <- paste0("[", head(bins, -1), "-", tail(bins, -1), "]")
  
  # Compute null_distribution_df for each site
  distribution_df <- df_all %>%
    mutate(Interval = cut(Frequency, breaks = bins, include.lowest = TRUE, right = FALSE, labels = interval_labels)) %>%
    group_by(Site, Interval) %>%
    summarise(Count = n(), .groups = "drop") %>%
    complete(Site, Interval = interval_labels, fill = list(Count = 0)) %>%
    mutate(Interval = factor(Interval, levels = interval_labels))  # Ensure correct ordering
  
  # Plot the histogram for each site with standardized Y-axis (0 to 20)
  dd <- ggplot(distribution_df, aes(x = Interval, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    facet_wrap(~ Site) +  # Create separate plots per site
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, max(distribution_df$Count)) +  # Standardize Y-axis from 0 to 20
    labs(title = "Distribution of Species Frequencies per Site",
         x = "Frequency Interval",
         y = "Number of Species")
  
  
  # Essai sur la base des ARMS plutot que face de plaques ####
  
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[, colSums(data_mean) != 0]
  
  data_mean <- data_mean[, msp_list_filter]
  data_mean_pa <- vegan::decostand(data_mean, "pa")
  
  f_ARMS <- colSums(data_mean_pa)/nrow(data_mean_pa)
  hist(f_ARMS, breaks = 50)
  
  
  df <- data.frame(Frequency = f_ARMS)
  
  # Plot histogram
  ee <- ggplot(df, aes(x = Frequency)) +
    geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(title = "Histogram of Species Occurrence Frequencies",
         x = "Occurrence Frequency",
         y = "Number of Species") +
    theme_minimal()
  
  
  # Essai sur la base des sites ####
  
  data_mean_site <- data %>% 
    group_by(meta$triplicat) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_mean_site <- data_mean_site[,-1]
  
  
  data_mean_site <- data_mean_site[, colSums(data_mean_site) != 0]
  
  data_mean_site <- data_mean_site[, msp_list_filter]
  data_mean_site_pa <- vegan::decostand(data_mean_site, "pa")
  
  f_site <- colSums(data_mean_site_pa)/nrow(data_mean_site_pa)
  hist(f_ARMS, breaks = 50)
  
  
  df <- data.frame(Frequency = f_ARMS)
  
  # Plot histogram
  ff <- ggplot(df, aes(x = Frequency)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(title = "Histogram of Species Occurrence Frequencies",
         x = "Occurrence Frequency",
         y = "Number of Species") +
    theme_minimal()
  
  # Faire un graphique frequence moyenne d'occurence des espèce en fonction de l'emprise #####
  
  tab_f <- data.frame(f_plate = f,
                         f_ARMS = f_ARMS,
                         f_site = f_site)
  
  
  tab_f$species <- rownames(tab_f)
  # Convert to long format
  data_long <- tab_f %>%
    pivot_longer(cols = c(f_plate, f_ARMS, f_site), names_to = "Modality", values_to = "Value")
  
  # # Plot with ggplot2
  # ggplot(data_long, aes(x = factor(Modality, levels = c("f_plate", "f_ARMS", "f_site")), y = Value, color = Modality)) +
  #   geom_point(aes(shape = Modality), size = 3) +
  #   geom_line(aes(group = species), alpha = 0.5) + # Optional: add lines to connect points per species
  #   theme_minimal() +
  #   labs(title = "Species Values Across Different Modalities",
  #        x = "Modality",
  #        y = "Value") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   guides(color = "none")  # Hide redundant color legend
  # 
  # ggplot(data_long, aes(x = factor(Modality, levels = c("f_plate", "f_ARMS", "f_site")), y = Value, color = species)) +
  #   geom_point(aes(shape = Modality), size = 3) +
  #   geom_line(aes(group = species), alpha = 0.5) + # Optional: add lines to connect points per species
  #   theme_minimal() +
  #   labs(title = "Species Values Across Different Modalities",
  #        x = "Modality",
  #        y = "Value") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   guides(color = guide_legend(title = "Species"))  
  
  
  # Step 1: Calculate the mean values for each modality
  mean_values <- data_long %>%
    group_by(Modality) %>%
    summarise(Mean_Value = mean(Value, na.rm = TRUE))
  
  # Step 2: Perform ANOVA to compare mean values across modalities (optional)
  anova_result <- aov(Value ~ Modality, data = data_long)
  summary(anova_result)
  
  # Step 3: Pairwise comparisons using Student's t-test with Bonferroni correction
  pairwise_result <- pairwise.t.test(data_long$Value, data_long$Modality, p.adjust.method = "bonferroni")
  print(pairwise_result)
  
  # Step 4: Set the factor levels for the Modality column to ensure correct order in the plot
  data_long$Modality <- factor(data_long$Modality, levels = c("f_plate", "f_ARMS", "f_site"))
  
  # Step 5: Create a violin plot to compare the distribution of values across modalities
  ggplot(data_long, aes(x = Modality, y = Value, fill = Modality)) +
    geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot
    geom_point(data = mean_values, aes(x = Modality, y = Mean_Value), 
               shape = 18, color = "purple", size = 4) +  # Mean as grey diamonds
    labs(title = "Distribution of Species Occurrence Frequencies by Modality", 
         x = "Modality", 
         y = "Occurrence Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")  # Choose a color palette
  
  # Null Model based on ARMS ####
  
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[, colSums(data_mean) != 0]
  
  data_mean <- data_mean[, msp_list_filter]
  data_mean_pa <- vegan::decostand(data_mean, "pa")
  
  f <- colSums(data_mean_pa)/nrow(data_mean_pa)
  # Compute a null model ####
  
  # Number of simulations
  n_simulations <- 1000  # Adjust as needed
  community_size <- 3  # Each simulated community has 48 sites
  
  # Define frequency intervals
  intervals <- seq(0, 1, by = 0.02)
  interval_labels <- paste0("[", head(intervals, -1), "-", tail(intervals, -1), "]")
  
  # Initialize a table to count species in each interval
  interval_counts <- setNames(rep(0, length(interval_labels)), interval_labels)
  
  # i = 3
  df_list_sim <- list()
  # Loop over simulations
  for (i in 1:n_simulations) {
    # Generate a simulated presence-absence matrix for one community
    sim_matrix <- t(replicate(community_size, rbinom(length(f), 1, f)))
    
    # Associate species names
    colnames(sim_matrix) <- names(f)
    
    # community composition simulated
    community_composition_sim <- colSums(sim_matrix)
    
    
    # P/a transform
    community_composition_sim <- ifelse(community_composition_sim >= 1, 1, 0)
    
    f_sim <- ifelse(community_composition_sim >= 1, f, NA)
    
    
    # Create a data frame and add a column for the site name
    df_sim <- data.frame(Frequency = f_sim, Species = names(f), n_sim = rep(i, length(f_sim)))
    
    # Store in the list
    df_list_sim[[i]] <- df_sim 
    
    
    
    # Remove species with f_sim = 0
    f_sim <- f_sim[f_sim > 0]
 
 
       
       
    # Bin frequencies into intervals of 0.02
    bin_counts <- table(cut(f_sim, breaks = intervals, labels = interval_labels, include.lowest = TRUE))
    
    # Accumulate counts
    interval_counts <- interval_counts + bin_counts
  }
  
  frequency_list <- lapply(df_list_sim, function(df) df$Frequency)
  
  # Combiner toutes les colonnes extraites en un seul DataFrame
  result_df <- bind_cols(frequency_list)
  
  # Optionnel : Nommer les colonnes du résultat en fonction de l'index ou d'un autre critère
  colnames(result_df) <- paste0("Frequency_", seq_along(frequency_list))
  
  # Voir le tableau final
  head(result_df)
  
  mean_values <- rowMeans(result_df, na.rm = TRUE)

  
  # Résultat dans un vecteur

  
  mean_frequency_sim <- data.frame(Frequency = mean_values,
                                   Species = names(f))
  
  vv <- ggplot(mean_frequency_sim, aes(x = Frequency)) +
    geom_density(alpha = 0.3, adjust = 0.4) +
    theme_minimal()
  
  # Convert results to a data frame
  null_distribution_df <- data.frame(
    Interval = names(interval_counts),
    Count = (as.numeric(interval_counts)/n_simulations)
  )
  
  # Print results
  print(null_distribution_df)
  
  # Optional: Plot histogram of species counts per interval
  
  cc <- ggplot(null_distribution_df, aes(x = Interval, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Distribution of Species Frequencies", x = "Frequency Interval", y = "Number of Species")
  
  
  # plot a grid of SAD for each unique site ####
  
  
  # Extract unique sites
  sites <- unique(meta_mean$triplicat)
  
  # Initialize an empty list to store data frames
  df_list <- list()
  
  # site <- sites[1]
  
  # Loop through each site
  for (site in sites) {
    # Subset data for the current site
    data_site <- subset(data_mean_pa, meta_mean$triplicat == site)
    
    # Compute species occurrence frequency
    f_site <- colSums(data_site) / nrow(data_site)
    
    f_site <- ifelse(f_site > 0, 1, 0)
    
    f_global_dans_site <- ifelse(f_site >= 1, f, 0)
    
    
    # Create a data frame and add a column for the site name
    df_site <- data.frame(Frequency = f_global_dans_site, Site = site)
    
    # Store in the list
    df_list[[site]] <- df_site
  }
  
  # Combine all site data into one data frame
  df_all <- bind_rows(df_list)
  
  df_all <- df_all[df_all$Frequency != 0,]
  
  
  color_map <- c(
    "RODARMS1" = "#DD8D29",
    "RODARMS2" = "#DD8D29",
    "RODARMS3" = "#DD8D29",
    "P50ARMS1" = "darkblue",
    "P50ARMS2" = "darkblue",
    "P50ARMS3" = "darkblue",
    "RUNARMS1" = "#46ACC8",
    "RUNARMS2" = "#46ACC8",
    "RUNARMS3" = "#46ACC8",
    "RUNARMS4" = "#46ACC8",
    "RUNARMS5" = "#46ACC8",
    "RUNARMS6" = "#46ACC8",
    "RUNARMS7" = "#46ACC8",
    "RUNARMS8" = "#46ACC8",
    "RUNARMS9" = "#46ACC8"
  )
  
  xx <- ggplot(df_all, aes(x = Frequency, color = Site, fill = Site)) +
    geom_density(alpha = 0.3, adjust = 0.4) +
    facet_wrap(~ Site) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map) +
    labs(title = "Courbes de densité par site",
         x = "Frequency",
         y = "Densité") +
    theme_minimal()
  
  
  
  
  
  # Plot histograms for all sites using facet_wrap()
  # ggplot(df_all, aes(x = Frequency)) +
  #   geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black", alpha = 0.7) +
  #   facet_wrap(~ Site, scales = "free_y") +  # Create a grid of plots
  #   labs(title = "Histogram of Species Occurrence Frequencies per Site",
  #        x = "Occurrence Frequency",
  #        y = "Number of Species") +
  #   theme_minimal() +
  #   ylim(0, 20)
  # 
  
  # Define bin width
  bin_width <- 0.02
  bins <- seq(0, 1, by = bin_width)  # Define bin edges
  
  # Create custom interval labels
  interval_labels <- paste0("[", head(bins, -1), "-", tail(bins, -1), "]")
  
  # Compute null_distribution_df for each site
  distribution_df <- df_all %>%
    mutate(Interval = cut(Frequency, breaks = bins, include.lowest = TRUE, right = FALSE, labels = interval_labels)) %>%
    group_by(Site, Interval) %>%
    summarise(Count = n(), .groups = "drop") %>%
    complete(Site, Interval = interval_labels, fill = list(Count = 0)) %>%
    mutate(Interval = factor(Interval, levels = interval_labels))  # Ensure correct ordering
  
  # Plot the histogram for each site with standardized Y-axis (0 to 20)
  dd <- ggplot(distribution_df, aes(x = Interval, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    facet_wrap(~ Site) +  # Create separate plots per site
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, max(distribution_df$Count)) +  # Standardize Y-axis from 0 to 20
    labs(title = "Distribution of Species Frequencies per Site",
         x = "Frequency Interval",
         y = "Number of Species")
  
  
  #### Kolmogorov smirnov ####
  
  # Reconstruct simulated frequency values from null_distribution_df
  intervals <- seq(0, 1, by = 0.02)
  interval_labels <- paste0("[", head(intervals, -1), "-", tail(intervals, -1), "]")
  
  # Midpoint of each interval
  midpoints <- (head(intervals, -1) + tail(intervals, -1)) / 2
  
  # Repeat midpoints according to count (rounded to nearest integer just in case)
  null_freq_values <- rep(midpoints, times = round(null_distribution_df$Count))
  
  # Fréquences observées pour un site donné
  site_name <- "RUNARMS8"
  observed_freq <- df_all %>%
    filter(Site == site_name) %>%
    pull(Frequency)
  
  ks_result <- ks.test(observed_freq, null_freq_values)
  print(ks_result)
  
  ks_results <- lapply(unique(df_all$Site), function(site_name) {
    observed_freq <- df_all %>%
      filter(Site == site_name) %>%
      pull(Frequency)
    
    res <- ks.test(observed_freq, null_freq_values)
    data.frame(Site = site_name, D = res$statistic, p.value = res$p.value)
  })
  
  ks_results_df <- do.call(rbind, ks_results)
  print(ks_results_df)
  
  library(ggplot2)
  
  ggplot(ks_results_df, aes(x = Site, y = D)) +
    geom_bar(stat = "identity", fill = "tomato", color = "black") +
    labs(title = "Statistique D du test de Kolmogorov-Smirnov par site",
         x = "Site",
         y = "D (écart maximal entre distributions)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  mean(ks_results_df$D)
  
  return(ks_results_df)
}
  