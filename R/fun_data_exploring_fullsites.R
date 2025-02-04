#' Explore the cleaned data
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to the subseted raw data file
#' @export
#'
fun_data_exploring_fullsites <- function(data_and_meta_clean_fullsites){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)

  
  library(iNEXT)
  library(ggpubr)
  library(forcats)
  library(ggplot2)
  library(wesanderson)
  library(vegan)
  library(dplyr)

  
  meta_P50 <- subset(meta, campain == "P50ARMS") 
  data_P50 <- subset(data, meta$campain == "P50ARMS") 
  meta_RODA <- subset(meta, campain == "RODARMS") 
  data_RODA <- subset(data, meta$campain == "RODARMS") 
  meta_RUNA <- subset(meta, campain == "RUNARMS") 
  data_RUNA <- subset(data, meta$campain == "RUNARMS") 
  
  data_RUNA_pa <- t(decostand(data_RUNA, "pa"))
  data_RODA_pa <- t(decostand(data_RODA, "pa"))
  data_P50_pa <- t(decostand(data_P50, "pa"))
  
  # Assuming dat contains your three datasets
  dat <- list(
    RODA = data_RODA_pa, 
    P50 = data_P50_pa, 
    RUNA = data_RUNA_pa
  )
  
  # Process each dataset in the list
  dat_processed <- lapply(dat, function(df) {
    # Convert the data frame to a matrix
    mat <- as.matrix(df)
    
    # Verify dimensions match row and column names
    if (length(rownames(mat)) != nrow(mat)) {
      stop("Row names do not match the number of rows in the data.")
    }
    if (length(colnames(mat)) != ncol(mat)) {
      stop("Column names do not match the number of columns in the data.")
    }
    
    return(mat)
  })
  
  # Rename the list elements with appropriate names (if needed)
  names(dat_processed) <- c("RODA", "P50", "RUNA")
  
  # Testing alpha diversity using jacknife ####
  
  # Check the structure of the processed data
  str(dat_processed)
  t <- seq(1, 1000, by=10)
  obj <- iNEXT::iNEXT(dat_processed, datatype="incidence_raw", size=t)
  
  v1 <- ggiNEXT(obj, type=1, facet.var="Assemblage", color.var = "Assemblage") +
        theme_bw(base_size = 18) +
        theme(legend.position="right")

  
  # Extract species richness estimates from $AsyEst
  extrapolated_richness <- obj$AsyEst %>%
    filter(Diversity == "Species richness") %>%
    select(Assemblage, Estimator, s.e.)
  
  print(extrapolated_richness)

  
  # Add labels for extrapolated species richness
  v1 <- ggiNEXT(obj, type = 1, facet.var = "Assemblage", color.var = "Assemblage") +
    theme_bw(base_size = 18) +
    theme(legend.position = "right") +
    geom_text(
      data = extrapolated_richness,
      aes(
        x = Inf,  # Position at the far right of the plot
        y = Inf,  # Position at the top of the plot
        label = paste("Extrapolated Richness: ", round(Estimator, 1))
      ),
      hjust = 1.1, vjust = 1.1,  # Adjust position relative to plot edges
      inherit.aes = FALSE
    )
  
  
  # Create the combined plot

  v1_combined <- ggiNEXT(obj, type = 1, color.var = "Assemblage") +
    theme_bw(base_size = 18) +
    theme(legend.position = "right") +
    annotate("text", label = paste0("Extrapolated richness = ", 
                                    round(extrapolated_richness[2, 2], 1),
                                    " ± ",
                                    round(extrapolated_richness[2, 3], 1)), 
             x = 575, y = 130, size = 4, hjust = 0) +
    annotate("text", label = paste0("Extrapolated richness = ", 
                                    round(extrapolated_richness[1, 2], 1),
                                    " ± ",
                                    round(extrapolated_richness[1, 3], 1)), 
             x = 575, y = 113, size = 4, hjust = 0) +
    annotate("text", label = paste0("Extrapolated richness = ", 
                                          round(extrapolated_richness[3, 2], 1),
                                          " ± ",
                                          round(extrapolated_richness[3, 3], 1)), 
             x = 575, y = 93, size = 4, hjust = 0) +
    labs(x = "Number of plate faces", y = "Morpho-species diversity")   +
    scale_color_manual(values = c("P50" = "darkblue", "RODA" = "#DD8D29", "RUNA" = "#46ACC8")) # Change the y-axis label here
  
  
  # Show the plot
  print(v1_combined)
  
  
  # Extract the endpoint positions for each assemblage from iNEXT object

  
  alpha_div_fullsites_path <- here::here("outputs/Species_acc_curves_fullsites.pdf")
  
  ggsave(alpha_div_fullsites_path, v1_combined, width = 9.5, height = 5.5)
  
  # Testing alpha diversity differences using bootstrap ####
  meta <- meta %>%
    mutate(
      clust = case_when(
        triplicat == "RUNARMS1" ~ 1,
        triplicat %in% c("RUNARMS2", "RUNARMS3") ~ 2,
        triplicat %in% c("RUNARMS4", "RUNARMS5", "RUNARMS6") ~ 3,
        triplicat %in% c("RUNARMS7", "RUNARMS8", "RUNARMS9") ~ 4,
        TRUE ~ NA_integer_
      )
    )
  

  
  # Function to check if all clusters are represented
  is_valid_sample <- function(sampled_data) {
    all(1:4 %in% sampled_data$clust)
  }
  
  # Bootstrap alpha diversity function
  bootstrap_alpha_div <- function(data, metadata, num_bootstraps = 1000) {
    # Filter shallow metadata
    shallow_meta <- metadata %>%
      filter(campain == "RUNARMS")
    deep_meta <- metadata %>%
      filter(campain == "P50ARMS")
    rodri_meta <- metadata %>%
      filter(campain == "RODARMS")
    
    
    # Bootstrap loop
    results <- replicate(num_bootstraps, {
      valid_sample <- NULL
      
      # Generate a valid sample containing all clusters
      while (is.null(valid_sample)) {
        seq <- sample(unique(shallow_meta$arms), 9) # Randomly sample arms
        sampled <- shallow_meta %>%
          filter(arms %in% seq)
        
        # Check if all clusters are represented
        if (is_valid_sample(sampled)) {
          valid_sample <- sampled
        }
      }
      
      # Get indices for the valid sample
      shallow_names <- valid_sample$name
      deep_names <- deep_meta$name
      rodri_names <- rodri_meta$name
      
      full_names <- c(shallow_names, deep_names, rodri_names)
      
      sub_data <- data[full_names,]
      rownames(meta) <- meta$name
      sub_meta <- meta[full_names,]
  
      
      # Calculate species richness
      richness <- specnumber(sub_data, groups = sub_meta$campain)
      
      # Return results
      list(
        Richness_P50A = richness["P50ARMS"],
        Richness_RODA = richness["RODARMS"],
        Richness_RUNA = richness["RUNARMS"]
      )
    }, simplify = FALSE)
    
    # Combine results into a data frame
    result_df <- do.call(rbind, lapply(results, function(res) {
      data.frame(
        Richness_P50A = res$Richness_P50A,
        Richness_RODA = res$Richness_RODA,
        Richness_RUNA = res$Richness_RUNA
      )
    }))
    
    return(result_df)
  }
  
  # Example usage
  bootstrap_alpha_results <- bootstrap_alpha_div(data, metadata)
  

  
  
  hist(bootstrap_alpha_results$Richness_RUNA, main = "Bootstrap Distribution of Total Richness", xlab = "Richness")
  mean_richness <- mean(bootstrap_alpha_results$Richness_RUNA)
  sd_richness <- sd(bootstrap_alpha_results$Richness_RUNA)

  # Plot the histogram with ggplot2
  alph <- ggplot(bootstrap_alpha_results, aes(x = Richness_RUNA)) +
    geom_histogram(binwidth = 1, fill = "coral", color = "black", alpha = 0.7) +
    geom_vline(xintercept = mean_richness, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = mean_richness, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 85, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 86, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = mean_richness - 10, y = 40, 
             label = paste0("Mean: ", round(mean_richness, 2)), 
             color = "black", hjust = 0) +
    annotate("text", x = mean_richness - 10, y = 35, 
             label = paste0("SD: ", round(sd_richness, 2)), 
             color = "black", hjust = 0) +
    labs(
      title = "Bootstrap Distribution of Total Richness",
      x = "Richness",
      y = "Frequency"
    ) +
    theme_minimal()
  
  
  alpha_div_bootstrap_path <- here::here("outputs/Species_richness_bootstrap_fullsites.pdf")
  
  ggsave(alpha_div_bootstrap_path, alph, width = 9, height = 9)
  
  
  # NMDS ####
  
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  fit <-  vegan::metaMDS(data_mean, distance = "bray")
  col1 <- c(wesanderson::wes_palette("FantasticFox1", 5), wesanderson::wes_palette("Royal2", 1))
  col1 <- c("darkblue","darkblue","darkblue", "#DD8D29","#DD8D29","#DD8D29", "#46ACC8", "#46ACC8", "#46ACC8", "#46ACC8", "#46ACC8", "#46ACC8", "#46ACC8", "#46ACC8", "#46ACC8")
  mds_name <- paste0("NMDS_MESO_bay_2.pdf")
  mds_path <- here::here("outputs/", mds_name)
  pdf(file =  mds_path)
  plot(fit)
  vegan::ordihull(fit, meta_mean$triplicat, col=col1)
  vegan::ordiellipse(fit, meta_mean$triplicat, col=col1, kind = "ehull", lwd=2)
  vegan::ordispider(fit, meta_mean$triplicat, col=col1, cex = 0.5)
  vegan::ordiellipse(fit, meta_mean$triplicat, col=col1, draw="polygon", label = FALSE)
  # vegan::ordilabel(fit, display = "species", choices = c(1, 2), cex = 0.3, border = NA)
  a <- paste0("stress = ",round(fit$stress, 3))
  text(0.8,0.95, a)
  dev.off()
  #Jac
  data_mean_pa <- vegan::decostand(data_mean, "pa")
  fit2 <-  vegan::metaMDS(data_mean_pa, distance = "jaccard")
  ?metaMDS
  mds_name <- paste0("NMDS_MESO_jac_2.pdf")
  mds_path <- here::here("outputs/", mds_name)
  pdf(file =  mds_path)
  plot(fit2)
  vegan::ordihull(fit2, meta_mean$triplicat, col=col1)
  vegan::ordiellipse(fit2, meta_mean$triplicat, col=col1, kind = "ehull", lwd=2)
  vegan::ordispider(fit2, meta_mean$triplicat, col=col1, cex = 0.5)
  vegan::ordiellipse(fit2, meta_mean$triplicat, col=col1, draw="polygon", label = FALSE)
  # vegan::ordilabel(fit2, display = "species", choices = c(1, 2), cex = 0.3, border = NA)
  a <- paste0("stress = ",round(fit2$stress, 3))
  text(0.6,1.3, a)
  dev.off()
  
  return(NULL)
}

  
