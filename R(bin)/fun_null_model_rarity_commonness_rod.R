# data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 

# Null Model based on ARMS ####

data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)

data <- subset(data, meta$campain == "RODARMS")
meta <- subset(meta, meta$campain == "RODARMS")


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

data <- data[, colSums(data) != 0]

data_pa <- vegan::decostand(data, "pa")

f <- colSums(data_pa)/nrow(data_pa)

df <- data.frame(Frequency = f)

# Plot histogram
aa <- ggplot(df, aes(x = Frequency)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Species Occurrence Frequencies",
       x = "Occurrence Frequency",
       y = "Number of Species") +
  theme_minimal()


# Compute a null model ####

# Number of simulations
n_simulations <- 1000  # Adjust as needed
community_size <- 16*3  # Each simulated community has 48 sites

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

long_freq <- pivot_longer(result_df, cols = everything(), names_to = "Iteration", values_to = "Frequency")

vv <- ggplot(long_freq, aes(x = Frequency)) +
  geom_density(alpha = 0.4, adjust = 2, na.rm = TRUE) +
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
sites <- unique(meta$triplicat)

# Initialize an empty list to store data frames
df_list <- list()

# site <- sites[1]

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
  geom_density(alpha = 0.3, adjust = 1) +
  facet_wrap(~ Site) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  labs(title = "Courbes de densité par site",
       x = "Frequency",
       y = "Densité") +
  theme_minimal()


# 1. Calculer les valeurs de densité globale
global_density <- density(na.omit(long_freq$Frequency), adjust = 5)
df_global_density <- data.frame(
  Frequency = global_density$x,
  Density = global_density$y
)

# 2. Ajouter une colonne fictive "Site" pour répliquer la densité sur tous les sites
unique_sites <- unique(df_all$Site)
df_global_density_all <- df_global_density %>%
  tidyr::crossing(Site = unique_sites)

# 3. Tracer les densités par site avec densité globale superposée
xx <- ggplot(df_all, aes(x = Frequency, color = Site, fill = Site)) +
  geom_density(alpha = 0.3, adjust = 1) +
  geom_line(data = df_global_density_all,
            aes(x = Frequency, y = Density),
            inherit.aes = FALSE,
            color = "black",  linewidth = 0.8) +
  facet_wrap(~ Site) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  labs(title = "Courbes de densité par site (avec densité globale en noir)",
       x = "Frequency",
       y = "Densité") +
  theme_minimal()

# Afficher le plot
print(xx)
