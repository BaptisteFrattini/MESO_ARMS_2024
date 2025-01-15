#' PERMANOVA and SIMPER
#'
#' @param data_and_meta_clean_fullsites the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_permanova_fullsites <- function(data_and_meta_clean_fullsites){
  
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  # PERMANOVA ####
  ## Effet de la profondeur ####
  
  data_mean_RUN <-  subset(data_mean, meta_mean$island == "Reunion")
  meta_mean_RUN <- subset(meta_mean, meta_mean$island == "Reunion")
  
  
  meta_mean_RUN <- data.frame(
    Sample = rownames(data_mean_RUN),
    Depth = ifelse(grepl("P50", rownames(data_mean_RUN)), "Deep", "Shallow")
  )
  

  bootstrap_permanova_jacc <- function(data, metadata, num_bootstraps = 1000) {
    # Identify indices of shallow and deep samples
    shallow_indices <- which(metadata$Depth == "Shallow")
    deep_indices <- which(metadata$Depth == "Deep")
    
    # Bootstrap loop
    results <- replicate(num_bootstraps, {
      # Randomly select 9 shallow samples with replacement
      selected_shallow <- sample(shallow_indices, size = 9, replace = TRUE)
      selected_indices <- c(selected_shallow, deep_indices)
      
      # Subset data and metadata
      sub_data <- data[selected_indices, ]
      sub_meta <- metadata[selected_indices, ]
      
      #Check dispersion
      dist <- vegan::vegdist(sub_data, "jaccard")
      a <- vegan::betadisper(
        dist,
        sub_meta$Depth,
        type = "median",
        bias.adjust = FALSE,
        sqrt.dist = FALSE,
        add = TRUE
      )
      mod <- anova(a)
      
      # Perform PERMANOVA
      adonis_result <- vegan::adonis2(sub_data ~ Depth, data = sub_meta, method = "jaccard")
      # anosim_results <- vegan::anosim(sub_data, sub_meta$Depth, distance = "jaccard")
      
      # Store results
      list(
        R2 = adonis_result$R2[1], 
        F = adonis_result$F[1],
        p_value = adonis_result$`Pr(>F)`[1],
        mean_permdisp_RUNA = a$group.distances[grepl("Shallow",names(a$group.distances))],
        mean_permdisp_P50A = a$group.distances[grepl("Deep",names(a$group.distances))],
        p_disp = mod$`Pr(>F)`,
        # anosim_R = anosim_results$statistic,
        # anosim_p = anosim_results$signif,
        samples = sub_meta$Sample  # Capture sample names
      )
    }, simplify = FALSE)
    
    # Format the results
    result_df <- do.call(rbind, lapply(results, function(res) {
      c(
        R2 = res$R2, 
        F = res$F, 
        p_value = res$p_value,
        mean_permdisp_RUNA = res$mean_permdisp_RUNA,
        mean_permdisp_P50A = res$mean_permdisp_P50A,
        p_disp = res$p_disp,
        # anosim_R = res$anosim_R,
        # anosim_p = res$anosim_p,
        samples = paste(res$samples, collapse = ", ")
      )
    }))
    
    # Convert to a data frame
    result_df <- as.data.frame(result_df, stringsAsFactors = FALSE)
    return(result_df)
  }
  
  
  
  bootstrap_permanova_bray <- function(data, metadata, num_bootstraps = 1000) {
    # Identify indices of shallow and deep samples
    shallow_indices <- which(metadata$Depth == "Shallow")
    deep_indices <- which(metadata$Depth == "Deep")
    
    # Bootstrap loop
    results <- replicate(num_bootstraps, {
      # Randomly select 9 shallow samples with replacement
      selected_shallow <- sample(shallow_indices, size = 9, replace = TRUE)
      selected_indices <- c(selected_shallow, deep_indices)
      
      # Subset data and metadata
      sub_data <- data[selected_indices, ]
      sub_meta <- metadata[selected_indices, ]
      
      #Check dispersion
      dist <- vegan::vegdist(sub_data, "bray")
      a <- vegan::betadisper(
        dist,
        sub_meta$Depth,
        type = "median",
        bias.adjust = FALSE,
        sqrt.dist = FALSE,
        add = TRUE
      )
      mod <- anova(a)
      
      # Perform PERMANOVA
      adonis_result <- vegan::adonis2(sub_data ~ Depth, data = sub_meta, method = "bray")
      # anosim_results <- vegan::anosim(sub_data, sub_meta$Depth, distance = "bray")
      
      # Store results
      list(
        R2 = adonis_result$R2[1], 
        F = adonis_result$F[1],
        p_value = adonis_result$`Pr(>F)`[1],
        mean_permdisp_RUNA = a$group.distances[grepl("Shallow",names(a$group.distances))],
        mean_permdisp_P50A = a$group.distances[grepl("Deep",names(a$group.distances))],
        p_disp = mod$`Pr(>F)`,
        # anosim_R = anosim_results$statistic,
        # anosim_p = anosim_results$signif,
        samples = sub_meta$Sample  # Capture sample names
      )
    }, simplify = FALSE)
    
    # Format the results
    result_df <- do.call(rbind, lapply(results, function(res) {
      c(
        R2 = res$R2, 
        F = res$F, 
        p_value = res$p_value,
        mean_permdisp_RUNA = res$mean_permdisp_RUNA,
        mean_permdisp_P50A = res$mean_permdisp_P50A,
        p_disp = res$p_disp,
        # anosim_R = res$anosim_R,
        # anosim_p = res$anosim_p,
        samples = paste(res$samples, collapse = ", ")
      )
    }))
    
    # Convert to a data frame
    result_df <- as.data.frame(result_df, stringsAsFactors = FALSE)
    return(result_df)
  }
  
  
  
  # Run the bootstrap analysis
  data_mean_RUN_pa <- vegan::decostand(data_mean_RUN, "pa")
  bootstrap_jacc_results <- bootstrap_permanova_jacc(data_mean_RUN_pa, meta_mean_RUN)
  bootstrap_bray_results <- bootstrap_permanova_bray(data_mean_RUN, meta_mean_RUN)
  
  write.csv(bootstrap_jacc_results, file = "outputs/table_bootstrap_permanova_jaccard.csv", row.names = TRUE)
  write.csv(bootstrap_bray_results, file = "outputs/table_bootstrap_permanova_bray.csv", row.names = TRUE)
  #Dispersion
  
  dist <- vegan::vegdist(data_mean_RUN, "bray")
  a <- vegan::betadisper(
    dist,
    meta_mean_RUN$Depth,
    type = "median",
    bias.adjust = TRUE,
    sqrt.dist = FALSE,
    add = TRUE
  )
  
  
  mod <- anova(a)
  
  anosim_results <- vegan::anosim(data_mean_RUN, meta_mean_RUN$Depth, distance = "bray")

  anosim_results
  
  # View results
  # t_test <- t.test(as.numeric(bootstrap_bray_results$mean_permdisp_RUNA.Shallow), mu = a$group.distances[1])
  # mean(as.numeric(bootstrap_bray_results$mean_permdisp_RUNA))

  #Dissimilarité des communautés
  
  # range(bootstrap_bray_results$R2)
  # range(bootstrap_jacc_results$R2)
  
  
  return(NULL)
}