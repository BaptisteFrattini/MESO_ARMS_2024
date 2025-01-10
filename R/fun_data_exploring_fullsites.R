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
    select(Assemblage, Estimator)
  
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
                                    round(extrapolated_richness[2, 2], 1)), 
             x = 575, y = 130, size = 4, hjust = 0) +
    annotate("text", label = paste0("Extrapolated richness = ", 
                                    round(extrapolated_richness[1, 2], 1)), 
             x = 575, y = 113, size = 4, hjust = 0) +
    annotate("text", label = paste0("Extrapolated richness = ", 
                                    round(extrapolated_richness[3, 2], 1)), 
             x = 575, y = 93, size = 4, hjust = 0) +
    labs(x = "Number of plate faces", y = "Morpho-species diversity") # Change the y-axis label here
  
  
  # Show the plot
  print(v1_combined)
  
  
  # Extract the endpoint positions for each assemblage from iNEXT object

  
  alpha_div_fullsites_path <- here::here("outputs/Species_acc_curves_fullsites.pdf")
  
  ggsave(alpha_div_fullsites_path, v1_combined, width = 7, height = 4.5)
  
  return(NULL)
}

  
