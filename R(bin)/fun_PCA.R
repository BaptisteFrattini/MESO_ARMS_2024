#' 
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_PCoA <- function(data_and_meta_clean){
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  library(ggpubr)
  library(forcats)
  library(reshape2)
  library(ggplot2)
  library(plotly)
  library(ggfortify)
  library(factoextra)
  library(dplyr)
  library(BiodiversityR)
  data <- read.csv(data_and_meta_clean["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean["path_meta"], row.names = 1)
  
  # PERMANOVA Test + SIMPER ####
  data_P50 <- subset(data, meta$campain == "P50ARMS")
  meta_P50 <- subset(meta, meta$campain == "P50ARMS")
  
  perm <- vegan::adonis2(data_P50 ~ orientation*open_close, 
                         data = meta_P50, 
                         method = "bray", 
                         strata = meta_P50$triplicat, 
                         permutations = 999)
  perm
  summary(perm)
  
  perm <- vegan::adonis2(data ~ orientation*campain, 
                         data = meta, 
                         method = "bray",
                         permutations = 999)
  
  # pairwise_ado <- pairwise.adonis2(df_mean ~ imm_time, data = meta_mean, method = "bray")
  
  
  # PCA ####
  
    ## orientation + campain ####
  
  bc <- vegan::vegdist(data, method = "bray")
  
  pca.res.bc <- prcomp(bc,
                       center = TRUE,
                       scale. = TRUE)
  
  plot_pca <- ggplot2::autoplot(pca.res.bc,
                                data = meta,
                                colour = "campain",
                                shape = "orientation", frame = TRUE)
  
    ## orientation-arms ####
  
  meta$orientation_arms <- paste0(meta$arms,"_", meta$orientation)
  
  data_mean_orient <- data %>% 
    group_by(meta$orientation_arms) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_mean_orient <- as.data.frame(data_mean_orient)
  row.names(data_mean_orient) <- data_mean_orient$`meta$orientation_arms`
  data_mean_orient <- data_mean_orient[,-1]
 
  name <- rownames(data_mean_orient)
  arms <- substr(rownames(data_mean_orient), 1,9)
  campain <- substr(rownames(data_mean_orient), 1,7)
  triplicat <- substr(rownames(data_mean_orient), 1,8) 
  orientation <- substr(rownames(data_mean_orient), 11, 11)
  meta_mean_orient <- data.frame(name, arms, campain, triplicat,orientation)
  meta_mean_orient$arms_orient <- paste0(meta_mean_orient$triplicat,"_", meta_mean_orient$orientation)
  
  bc <- vegan::vegdist(data_mean_orient, method = "bray")
  
  pca.res.bc <- prcomp(bc,
                       center = TRUE,
                       scale. = TRUE)
  
  plot_pca <- ggplot2::autoplot(pca.res.bc,
                                data = meta_mean_orient,
                                colour = "arms_orient", frame = TRUE)
  
      ###PCoA ####
  df_bray <- vegan::vegdist(data_mean_orient, method = "bray")
  pcoa_res <- ape::pcoa(df_bray)
  
  sites.scores <- pcoa_res$vectors[,c(1,2)]
  
  pcoa_res2 <- cmdscale(df_bray)
  species.scores <- BiodiversityR::add.spec.scores(pcoa_res2 , data_mean_orient, method="cor.scores", multi=1, Rscale=F, scaling="1")
  species.scores <- species.scores$cproj
  nrow(species.scores)

  
  sim_imm_tim <- summary(vegan::simper(data_mean_orient, meta_mean_orient$arms_orient, permutations = 999))
 
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$cumsum < 0.95),]
  
  sim_imm_tim_2 <- sim_imm_tim[[2]]
  contrib_imm_tim_2 <- sim_imm_tim_2[(sim_imm_tim_2$cumsum < 0.95),]
  
  sim_imm_tim_3 <- sim_imm_tim[[3]]
  contrib_imm_tim_3 <- sim_imm_tim_3[(sim_imm_tim_3$cumsum < 0.95),]
  
  sel <- levels(as.factor(c(rownames(contrib_imm_tim_1), rownames(contrib_imm_tim_2), rownames(contrib_imm_tim_3))))
  
  species.scores <- species.scores[sel,]
  species.scores <- species.scores*0.5
  
  
  
  colnames(sites.scores) <- c("Dim1", "Dim2")
  
  biplot1 <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_mean_orient$campain, shape = meta_mean_orient$orientation, size = 3)) +
    scale_color_manual(values = c("P50ARMS" = "darkblue", "RUNARMS" = "blue")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = NA, color = meta_mean_orient$campain, fontface = "bold"), vjust = -1.5) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  
  
  
  
  # PCoA data_mean ####
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  return(NULL)
  
}