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
  
  # Shallow VS Deep ####
  ## PERMANOVA ####

  
  data_mean_RUN <-  subset(data_mean, meta_mean$island == "Reunion")
  meta_mean_RUN <- subset(meta_mean, meta_mean$island == "Reunion")
  
  
  meta_mean_RUN <- data.frame(
    Sample = rownames(data_mean_RUN),
    Depth = ifelse(grepl("P50", rownames(data_mean_RUN)), "Deep", "Shallow")
  )
  
  meta_mean_RUN$clust <- c(rep("NA",9),rep(1,3),rep(2,6),rep(3,9), rep(4,9))
 
  is_valid_sample <- function(sampled_data) {
    all(1:4 %in% sampled_data$clust)
  }
  
  
  bootstrap_permanova_jacc <- function(data, metadata, num_bootstraps = 1000) {
    # Identify indices of shallow and deep samples
    shallow_meta <- subset(metadata, Depth == "Shallow")
    deep_indices <- which(metadata$Depth == "Deep")
    
    # Bootstrap loop
    results <- replicate(num_bootstraps, {
      # Randomly select 9 shallow samples with replacement
      
      #insert here the loop for 
      
      valid_sample <- NULL
      
      while (is.null(valid_sample)) {
        # Randomly sample 9 rows
        sampled <- shallow_meta[sample(nrow(shallow_meta), 9), ]
        
        # Check if the sampled subset contains all clusters
        if (is_valid_sample(sampled)) {
          valid_sample <- sampled
        }
      }
      
      shallow_indices <- as.numeric(rownames(valid_sample))
      
      selected_indices <- c(shallow_indices, deep_indices)
      
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
    shallow_meta <- subset(metadata, Depth == "Shallow")
    deep_indices <- which(metadata$Depth == "Deep")
    
    # Bootstrap loop
    results <- replicate(num_bootstraps, {
      # Randomly select 9 shallow samples with replacement
      
      #insert here the loop for 
      
      valid_sample <- NULL
      
      while (is.null(valid_sample)) {
        # Randomly sample 9 rows
        sampled <- shallow_meta[sample(nrow(shallow_meta), 9), ]
        
        # Check if the sampled subset contains all clusters
        if (is_valid_sample(sampled)) {
          valid_sample <- sampled
        }
      }
      
      shallow_indices <- as.numeric(rownames(valid_sample))
      
      selected_indices <- c(shallow_indices, deep_indices)
      
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
  mean(as.numeric(bootstrap_bray_results$R2))
  sd(as.numeric(bootstrap_bray_results$R2))
  hist(as.numeric(bootstrap_bray_results$R2), breaks = 100)
  
  mean(as.numeric(bootstrap_jacc_results$R2))
  sd(as.numeric(bootstrap_jacc_results$R2))
  hist(as.numeric(bootstrap_jacc_results$R2), breaks = 100)
  
  range(as.numeric(bootstrap_bray_results$p_value))
  hist(as.numeric(bootstrap_bray_results$p_value))
  
  range(as.numeric(bootstrap_jacc_results$p_value))
  hist(as.numeric(bootstrap_jacc_results$p_value))
  
  
  library(ggplot2)
  library(dplyr)
  library(hrbrthemes)
  library(data.table)
  
  tab_bray <- data.frame(shallow = as.numeric(bootstrap_bray_results$mean_permdisp_RUNA.Shallow),
                         deep = as.numeric(bootstrap_bray_results$mean_permdisp_P50A.Deep))
  tab_bray <- tab_bray %>%
    tidyr::pivot_longer(cols = c(shallow, deep), 
                 names_to = "depth", 
                 values_to = "permdisp_value")
  
  
  p <- tab_bray %>%
    ggplot( aes(x=permdisp_value, fill=depth)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill="")
  
  tab_jacc <- data.frame(shallow = as.numeric(bootstrap_jacc_results$mean_permdisp_RUNA.Shallow),
                         deep = as.numeric(bootstrap_jacc_results$mean_permdisp_P50A.Deep))
  tab_jacc <- tab_jacc %>%
    tidyr::pivot_longer(cols = c(shallow, deep), 
                        names_to = "depth", 
                        values_to = "permdisp_value")
  

  t <- tab_jacc %>%
    ggplot( aes(x=permdisp_value, fill=depth)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
    theme_ipsum() +
    labs(fill="")
  
  mean(as.numeric(bootstrap_bray_results$p_disp1))
  sd(as.numeric(bootstrap_bray_results$p_disp1))
  hist(as.numeric(bootstrap_bray_results$p_disp1), breaks = 100)
  
  mean(as.numeric(bootstrap_jacc_results$p_disp1))
  sd(as.numeric(bootstrap_jacc_results$p_disp1))
  hist(as.numeric(bootstrap_jacc_results$p_disp1), breaks = 100)
  
  # View results
  # t_test <- t.test(as.numeric(bootstrap_bray_results$mean_permdisp_RUNA.Shallow), mu = a$group.distances[1])
  # mean(as.numeric(bootstrap_bray_results$mean_permdisp_RUNA))

  #Dissimilarité des communautés
  
  # range(bootstrap_bray_results$R2)
  # range(bootstrap_jacc_results$R2)
  
  # Reunion VS Rodrigues ####
  ## PERMANOVA ####
  
  meta_mean_Shallow <- subset(meta_mean, campain %in% c("RODARMS", "RUNARMS"))
  data_mean_Shallow <- subset(data_mean, meta_mean$campain %in% c("RODARMS", "RUNARMS"))
  
  meta_mean_Shallow <- data.frame(
    Sample = rownames(data_mean_Shallow),
    Island = ifelse(grepl("ROD", rownames(data_mean_Shallow)), "Rodrigues", "Reunion")
  )
  
  meta_mean_Shallow$clust <- c(rep("NA",8),rep(1,3),rep(2,6),rep(3,9), rep(4,9))
  
  is_valid_sample <- function(sampled_data) {
    all(1:4 %in% sampled_data$clust)
  }
  
  
  bootstrap_permanova_jacc <- function(data, metadata, num_bootstraps = 1000) {
    # Identify indices of shallow and deep samples
    Reunion_meta <- subset(metadata, Island == "Reunion")
    Rodrigues_indices <- which(metadata$Island == "Rodrigues")
    
    # Bootstrap loop
    results <- replicate(num_bootstraps, {
      # Randomly select 9 shallow samples with replacement
      
      #insert here the loop for 
      
      valid_sample <- NULL
      
      while (is.null(valid_sample)) {
        # Randomly sample 9 rows
        sampled <- Reunion_meta[sample(nrow(Reunion_meta), 9), ]
        
        # Check if the sampled subset contains all clusters
        if (is_valid_sample(sampled)) {
          valid_sample <- sampled
        }
      }
      
      Reunion_indices <- as.numeric(rownames(valid_sample))
      
      selected_indices <- c(Reunion_indices, Rodrigues_indices)
      
      # Subset data and metadata
      sub_data <- data[selected_indices, ]
      sub_meta <- metadata[selected_indices, ]
      
      #Check dispersion
      dist <- vegan::vegdist(sub_data, "jaccard")
      a <- vegan::betadisper(
        dist,
        sub_meta$Island,
        type = "median",
        bias.adjust = FALSE,
        sqrt.dist = FALSE,
        add = TRUE
      )
      mod <- anova(a)
      
      # Perform PERMANOVA
      adonis_result <- vegan::adonis2(sub_data ~ Island, data = sub_meta, method = "jaccard")
      # anosim_results <- vegan::anosim(sub_data, sub_meta$Depth, distance = "jaccard")
      
      # Store results
      list(
        R2 = adonis_result$R2[1], 
        F = adonis_result$F[1],
        p_value = adonis_result$`Pr(>F)`[1],
        mean_permdisp_Reunion = a$group.distances[grepl("Reunion",names(a$group.distances))],
        mean_permdisp_Rodrigues = a$group.distances[grepl("Rodrigues",names(a$group.distances))],
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
        mean_permdisp_Reunion = res$mean_permdisp_Reunion,
        mean_permdisp_Rodrigues = res$mean_permdisp_Rodrigues,
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
    Reunion_meta <- subset(metadata, Island == "Reunion")
    Rodrigues_indices <- which(metadata$Island == "Rodrigues")
    
    # Bootstrap loop
    results <- replicate(num_bootstraps, {
      # Randomly select 9 shallow samples with replacement
      
      #insert here the loop for 
      
      valid_sample <- NULL
      
      while (is.null(valid_sample)) {
        # Randomly sample 9 rows
        sampled <- Reunion_meta[sample(nrow(Reunion_meta), 9), ]
        
        # Check if the sampled subset contains all clusters
        if (is_valid_sample(sampled)) {
          valid_sample <- sampled
        }
      }
      
      Reunion_indices <- as.numeric(rownames(valid_sample))
      
      selected_indices <- c(Reunion_indices, Rodrigues_indices)
      
      # Subset data and metadata
      sub_data <- data[selected_indices, ]
      sub_meta <- metadata[selected_indices, ]
      
      #Check dispersion
      dist <- vegan::vegdist(sub_data, "bray")
      a <- vegan::betadisper(
        dist,
        sub_meta$Island,
        type = "median",
        bias.adjust = FALSE,
        sqrt.dist = FALSE,
        add = TRUE
      )
      mod <- anova(a)
      
      # Perform PERMANOVA
      adonis_result <- vegan::adonis2(sub_data ~ Island, data = sub_meta, method = "bray")
      # anosim_results <- vegan::anosim(sub_data, sub_meta$Depth, distance = "bray")
      
      # Store results
      list(
        R2 = adonis_result$R2[1], 
        F = adonis_result$F[1],
        p_value = adonis_result$`Pr(>F)`[1],
        mean_permdisp_Reunion = a$group.distances[grepl("Reunion",names(a$group.distances))],
        mean_permdisp_Rodrigues = a$group.distances[grepl("Rodrigues",names(a$group.distances))],
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
        mean_permdisp_Reunion = res$mean_permdisp_Reunion,
        mean_permdisp_Rodrigues = res$mean_permdisp_Rodrigues,
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
  data_mean_Shallow_pa <- vegan::decostand(data_mean_Shallow, "pa")
  bootstrap_jacc_results <- bootstrap_permanova_jacc(data_mean_Shallow_pa, meta_mean_Shallow)
  bootstrap_bray_results <- bootstrap_permanova_bray(data_mean_Shallow, meta_mean_Shallow)
  
  write.csv(bootstrap_jacc_results, file = "outputs/table_bootstrap_permanova_jaccard_island.csv", row.names = TRUE)
  write.csv(bootstrap_bray_results, file = "outputs/table_bootstrap_permanova_bray_island.csv", row.names = TRUE)
  #Dispersion
  mean(as.numeric(bootstrap_bray_results$R2))
  sd(as.numeric(bootstrap_bray_results$R2))
  hist(as.numeric(bootstrap_bray_results$R2), breaks = 100)
  
  mean(as.numeric(bootstrap_jacc_results$R2))
  sd(as.numeric(bootstrap_jacc_results$R2))
  hist(as.numeric(bootstrap_jacc_results$R2), breaks = 100)
  
  range(as.numeric(bootstrap_bray_results$p_value))
  hist(as.numeric(bootstrap_bray_results$p_value))
  
  range(as.numeric(bootstrap_jacc_results$p_value))
  hist(as.numeric(bootstrap_jacc_results$p_value))
  
  
  library(ggplot2)
  library(dplyr)
  library(hrbrthemes)
  library(data.table)
  
  tab_bray <- data.frame(Rodrigues = as.numeric(bootstrap_bray_results$mean_permdisp_Rodrigues.Rodrigues),
                         Reunion = as.numeric(bootstrap_bray_results$mean_permdisp_Reunion.Reunion))
  tab_bray <- tab_bray %>%
    tidyr::pivot_longer(cols = c(Reunion, Rodrigues), 
                        names_to = "island", 
                        values_to = "permdisp_value")
  
  
  p <- tab_bray %>%
    ggplot( aes(x=permdisp_value, fill=island)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#DD8D29")) +
    theme_ipsum() +
    labs(fill="")
  
  tab_jacc <- data.frame(Rodrigues = as.numeric(bootstrap_jacc_results$mean_permdisp_Rodrigues.Rodrigues),
                         Reunion = as.numeric(bootstrap_jacc_results$mean_permdisp_Reunion.Reunion))
  tab_jacc <- tab_jacc %>%
    tidyr::pivot_longer(cols = c(Reunion, Rodrigues), 
                        names_to = "island", 
                        values_to = "permdisp_value")
  
  
  t <- tab_jacc %>%
    ggplot( aes(x=permdisp_value, fill=island)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#DD8D29")) +
    theme_ipsum() +
    labs(fill="")
  
  mean(as.numeric(bootstrap_bray_results$p_disp1))
  sd(as.numeric(bootstrap_bray_results$p_disp1))
  hist(as.numeric(bootstrap_bray_results$p_disp1), breaks = 100)
  
  mean(as.numeric(bootstrap_jacc_results$p_disp1))
  sd(as.numeric(bootstrap_jacc_results$p_disp1))
  hist(as.numeric(bootstrap_jacc_results$p_disp1), breaks = 100)
  
  return(NULL)
}