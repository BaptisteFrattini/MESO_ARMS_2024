#' Beta diversity decomposing
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_decomp_betadiv <- function(data_and_meta_clean){
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  library(betapart)
  library(reshape2)
  library(stringr)
  library(dplyr)
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)

  data_mean_pa <- vegan::decostand(data_mean, "pa") 
  
  B.pair.pa <- betapart::beta.pair(data_mean_pa, index.family = "jaccard")
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- B.pair.pa$beta.jac
  
  ## Comparison between ARMS of the same depth and ARMS from =/= depth ####
    ### turn ####
  
  df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  df.turn <- subset(df.turn, row != col)
  
  # Convert factors to characters
  df.turn$row <- as.character(df.turn$row)
  df.turn$col <- as.character(df.turn$col)
  
  # Create a new column with sorted combinations
  df.turn$sorted_comparison <- apply(df.turn[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  
  # Identify and remove redundant comparisons
  df_unique <- df.turn %>% distinct(sorted_comparison, .keep_all = TRUE)
  
  # Remove the temporary column
  df.turn <- df_unique[, -ncol(df_unique)]
  
  merged_df <- merge(df.turn, meta_mean, by.x = c("row"), by.y = c("arms"))
  merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  
  df.turn <- data.frame(merged_df2$col, merged_df2$site.y, merged_df2$row, merged_df2$site.x, merged_df2$value) 
  df.turn$same_site <- ifelse(df.turn$merged_df2.site.y == df.turn$merged_df2.site.x, "Yes", "No")
  subset_df.turn <- subset(df.turn, same_site == "Yes")
  subset_df.turn$row2 <- substr(subset_df.turn$merged_df2.row, 1, 7)
  subset_df.turn$col2 <- substr(subset_df.turn$merged_df2.col, 1, 7)
  subset_df.turn$same_depth <- ifelse(subset_df.turn$row2 == subset_df.turn$col2, "Yes", "No")
  subset_df.turn <- data.frame(col = subset_df.turn$merged_df2.col, row = subset_df.turn$merged_df2.row, site = subset_df.turn$merged_df2.site.y, same_depth = subset_df.turn$same_depth, value = subset_df.turn$merged_df2.value)
  
  turn.mean = tapply(subset_df.turn$value, list(subset_df.turn$same_depth, subset_df.turn$site), mean)
  turn.sd = tapply(subset_df.turn$value, list(subset_df.turn$same_depth, subset_df.turn$site), sd)
  
  library(ggplot2)
  library(ggpubr)
  library(forcats)
  
  
  intrainter = c("between ARMS \n of the same depth", "between mesophotic \n ARMS and shallow ARMS")
  
  df_saint_leu <- subset(subset_df.turn, subset_df.turn$site == "Saint_Leu")
  df_saint_leu <- dplyr::distinct(df_saint_leu, col, row, .keep_all = TRUE)
  
  df_grand_bois <- subset(subset_df.turn, subset_df.turn$site == "Grand_Bois")
  df_grand_bois <- dplyr::distinct(df_grand_bois, col, row, .keep_all = TRUE)
  
  df_cap_houssaye <- subset(subset_df.turn, subset_df.turn$site == "Cap_La_Houssaye")
  df_cap_houssaye <- dplyr::distinct(df_cap_houssaye, col, row, .keep_all = TRUE)
  
  ggplot(df_saint_leu, aes(x=value)) + 
    geom_density() #Données pas normales normales
  
  p.sed <- rstatix::wilcox_test(df_saint_leu, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  a1 <- ggplot(df_saint_leu, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("coral","coral") ) +
    labs(title = NULL,
         x = "Comparisons",
         y = "Turnover component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = NULL) +
    theme_classic() + 
    scale_y_continuous(limits = c(0, 0.9)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  a1
  
  ggplot(df_grand_bois, aes(x=value)) + 
    geom_density() #Données pas normales normales
  shapiro.test(df_grand_bois$value) #Shapiro not OK 
  
  p.sed <- rstatix::t_test(df_grand_bois, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  b1 <- ggplot(df_grand_bois, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("chartreuse4","chartreuse4") ) +
    labs(title = NULL,
         x = "Comparisons",
         y = "Turnover component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = NULL) +
    theme_classic() + 
    scale_y_continuous(limits = c(0, 0.9)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  b1
  
  ggplot(df_cap_houssaye, aes(x=value)) + 
    geom_density() #Données pas normales normales
  shapiro.test(df_cap_houssaye$value) #Shapiro not OK 
  
  p.sed <- rstatix::wilcox_test(df_cap_houssaye, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  c1 <- ggplot(df_cap_houssaye, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("blue3","blue3") ) +
    labs(title = NULL,
         x = "Comparisons",
         y = "Turnover component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = NULL) +
    theme_classic() + 
    scale_y_continuous(limits = c(0, 0.9)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  c1
  

    ### nest ####
  
  df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  df.nest <- subset(df.nest, row != col)
  
  # Convert factors to characters
  df.nest$row <- as.character(df.nest$row)
  df.nest$col <- as.character(df.nest$col)
  
  # Create a new column with sorted combinations
  df.nest$sorted_comparison <- apply(df.nest[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  
  # Identify and remove redundant comparisons
  df_unique <- df.nest %>% distinct(sorted_comparison, .keep_all = TRUE)
  
  # Remove the temporary column
  df.nest <- df_unique[, -ncol(df_unique)]
  
  
  merged_df <- merge(df.nest, meta_mean, by.x = c("row"), by.y = c("arms"))
  merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  
  df.nest <- data.frame(merged_df2$col, merged_df2$site.y, merged_df2$row, merged_df2$site.x, merged_df2$value) 
  df.nest$same_site <- ifelse(df.nest$merged_df2.site.y == df.nest$merged_df2.site.x, "Yes", "No")
  subset_df.nest <- subset(df.nest, same_site == "Yes")
  subset_df.nest$row2 <- substr(subset_df.nest$merged_df2.row, 1, 7)
  subset_df.nest$col2 <- substr(subset_df.nest$merged_df2.col, 1, 7)
  subset_df.nest$same_depth <- ifelse(subset_df.nest$row2 == subset_df.nest$col2, "Yes", "No")
  subset_df.nest <- data.frame(col = subset_df.nest$merged_df2.col, row = subset_df.nest$merged_df2.row, site = subset_df.nest$merged_df2.site.y, same_depth = subset_df.nest$same_depth, value = subset_df.nest$merged_df2.value)
  
  nest.mean = tapply(subset_df.nest$value, list(subset_df.nest$same_depth, subset_df.nest$site), mean)
  nest.sd = tapply(subset_df.nest$value, list(subset_df.nest$same_depth, subset_df.nest$site), sd)
  
  intrainter = c("between ARMS \n of the same depth", "between mesophotic \n ARMS and shallow ARMS")
  
  df_saint_leu <- subset(subset_df.nest, subset_df.nest$site == "Saint_Leu")
  df_saint_leu <- dplyr::distinct(df_saint_leu, col, row, .keep_all = TRUE)
  
  df_grand_bois <- subset(subset_df.nest, subset_df.nest$site == "Grand_Bois")
  df_grand_bois <- dplyr::distinct(df_grand_bois, col, row, .keep_all = TRUE)
  
  df_cap_houssaye <- subset(subset_df.nest, subset_df.nest$site == "Cap_La_Houssaye")
  df_cap_houssaye <- dplyr::distinct(df_cap_houssaye, col, row, .keep_all = TRUE)
  
  ggplot(df_saint_leu, aes(x=value)) + 
    geom_density() #Données pas normales normales
  shapiro.test(df_saint_leu$value) #Shapiro not OK 
  
  p.sed <- rstatix::wilcox_test(df_saint_leu, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  a2 <- ggplot(df_saint_leu, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("coral","coral") ) +
    labs(title = NULL,
         x = "Comparisons",
         y = "Nestedness component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intrainter) +
    theme_classic() + 
    scale_y_continuous(limits = c(0, 0.25)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  a2
  
  ggplot(df_grand_bois, aes(x=value)) + 
    geom_density() #Données pas normales normales
  shapiro.test(df_grand_bois$value) #Shapiro not OK 
  
  p.sed <- rstatix::wilcox_test(df_grand_bois, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  b2 <- ggplot(df_grand_bois, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("chartreuse4","chartreuse4") ) +
    labs(title = NULL,
         x = "Comparisons",
         y = "Nestedness component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intrainter) +
    theme_classic() + 
    scale_y_continuous(limits = c(0, 0.25)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  b2
  
  ggplot(df_cap_houssaye, aes(x=value)) + 
    geom_density() #Données pas normales normales
  shapiro.test(df_cap_houssaye$value) #Shapiro not OK 
  
  p.sed <- rstatix::wilcox_test(df_cap_houssaye, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  c2 <- ggplot(df_cap_houssaye, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("blue3","blue3") ) +
    labs(title = NULL,
         x = "Comparisons",
         y = "Nestedness component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intrainter) +
    theme_classic() + 
    scale_y_continuous(limits = c(0, 0.25)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  c2
  
    ### jacc ####
  
  df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  df.jacc <- subset(df.jacc, row != col)
  
  # Convert factors to characters
  df.jacc$row <- as.character(df.jacc$row)
  df.jacc$col <- as.character(df.jacc$col)
  
  # Create a new column with sorted combinations
  df.jacc$sorted_comparison <- apply(df.jacc[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  
  # Identify and remove redundant comparisons
  df_unique <- df.jacc %>% distinct(sorted_comparison, .keep_all = TRUE)
  
  # Remove the temporary column
  df.jacc <- df_unique[, -ncol(df_unique)]
  
  
  
  merged_df <- merge(df.jacc, meta_mean, by.x = c("row"), by.y = c("arms"))
  merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  
  df.jacc <- data.frame(merged_df2$col, merged_df2$site.y, merged_df2$row, merged_df2$site.x, merged_df2$value) 
  df.jacc$same_site <- ifelse(df.jacc$merged_df2.site.y == df.jacc$merged_df2.site.x, "Yes", "No")
  subset_df.jacc <- subset(df.jacc, same_site == "Yes")
  subset_df.jacc$row2 <- substr(subset_df.jacc$merged_df2.row, 1, 7)
  subset_df.jacc$col2 <- substr(subset_df.jacc$merged_df2.col, 1, 7)
  subset_df.jacc$same_depth <- ifelse(subset_df.jacc$row2 == subset_df.jacc$col2, "Yes", "No")
  subset_df.jacc <- data.frame(col = subset_df.jacc$merged_df2.col, row = subset_df.jacc$merged_df2.row, site = subset_df.jacc$merged_df2.site.y, same_depth = subset_df.jacc$same_depth, value = subset_df.jacc$merged_df2.value)
  
  jacc.mean = tapply(subset_df.jacc$value, list(subset_df.jacc$same_depth, subset_df.jacc$site), mean)
  jacc.sd = tapply(subset_df.jacc$value, list(subset_df.jacc$same_depth, subset_df.jacc$site), sd)
  
  intrainter = c("between ARMS \n of the same depth", "between mesophotic \n ARMS and shallow ARMS")
  
  df_saint_leu <- subset(subset_df.jacc, subset_df.jacc$site == "Saint_Leu")
  df_saint_leu <- dplyr::distinct(df_saint_leu, col, row, .keep_all = TRUE)
  
  df_grand_bois <- subset(subset_df.jacc, subset_df.jacc$site == "Grand_Bois")
  df_grand_bois <- dplyr::distinct(df_grand_bois, col, row, .keep_all = TRUE)
  
  df_cap_houssaye <- subset(subset_df.jacc, subset_df.jacc$site == "Cap_La_Houssaye")
  df_cap_houssaye <- dplyr::distinct(df_cap_houssaye, col, row, .keep_all = TRUE)
  
  ggplot(df_saint_leu, aes(x=value)) + 
    geom_density() #Données pas normales normales
  shapiro.test(df_saint_leu$value) #Shapiro not OK 
  
  p.sed <- rstatix::wilcox_test(df_saint_leu, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  a3 <- ggplot(df_saint_leu, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("coral","coral") ) +
    labs(title = "Saint Leu",
         x = "Comparisons",
         y = "Jaccard component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = NULL) +
    theme_classic() + 
    scale_y_continuous(limits = c(0, 0.9)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  a3
  
  ggplot(df_grand_bois, aes(x=value)) + 
    geom_density() #Données pas normales normales
  shapiro.test(df_grand_bois$value) #Shapiro not OK 
  
  p.sed <- rstatix::wilcox_test(df_grand_bois, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  b3 <- ggplot(df_grand_bois, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("chartreuse4","chartreuse4") ) +
    labs(title = "Grand Bois",
         x = "Comparisons",
         y = "Jaccard component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = NULL) +
    theme_classic() + 
    scale_y_continuous(limits = c(0,0.9)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  b3
  
  ggplot(df_cap_houssaye, aes(x=value)) + 
    geom_density() #Données pas normales normales
  shapiro.test(df_cap_houssaye$value) #Shapiro not OK 
  
  p.sed <- rstatix::t_test(df_cap_houssaye, value ~ same_depth) #parametrique
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  c3 <- ggplot(df_cap_houssaye, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
    geom_boxplot(fill =  c("blue3","blue3") ) +
    labs(title = "Cap La Houssaye",
         x = "Comparisons",
         y = "Jaccard component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=NULL) +
    theme_classic() + 
    scale_y_continuous(limits = c(0, 0.9)) +
    stat_pvalue_manual(p.sed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  c3
  
  fin <- cowplot::plot_grid(c3, a3, b3, c1, a1, b1, c2, a2, b2, 
                            ncol = 3,
                            nrow = 3)
  
  path_to_boxplot_betadiv <- paste0("outputs/boxplot_beta_decomp.pdf")
  ggsave(filename =  path_to_boxplot_betadiv, plot = fin, width = 13, height = 15.5)
  
  
  
  ## Comparison between ARMS of differents site and ARMS of different depth ####
    ### turn ####
  df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  df.turn <- subset(df.turn, row != col)
  
  # Convert factors to characters
  df.turn$row <- as.character(df.turn$row)
  df.turn$col <- as.character(df.turn$col)
  
  # Create a new column with sorted combinations
  df.turn$sorted_comparison <- apply(df.turn[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  
  # Identify and remove redundant comparisons
  df_unique <- df.turn %>% distinct(sorted_comparison, .keep_all = TRUE)
  
  # Remove the temporary column
  
  df.turn <- df_unique[, -ncol(df_unique)]
  
  
  merged_df <- merge(df.turn, meta_mean, by.x = c("row"), by.y = c("arms"))
  merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  df.turn <- subset(merged_df2, campain.x == campain.y)
  df.turn$comparisons <- ifelse(df.turn$site.x == df.turn$site.y, "Same_site", "Different_site")

  
  
  ggplot(df.turn, aes(x=value)) + 
    geom_density() #Données normales
  
  
  library(ggsignif)
  library(ggpubr)
  p3 <- ggboxplot(df.turn, x = "campain.x", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "campain.x",
                  palette = c("blue4", "aquamarine3"),
                  facet.by = "comparisons", ylab = "Turnover")+ stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) 
  
  p3 
  
  ### trying to include comparison between deep and shallow
  df.turn.2 <- subset(merged_df2, campain.x != campain.y)
  df.turn.2$comparisons <- ifelse(df.turn.2$site.x == df.turn.2$site.y, "Same_site", "Different_site")
  df.turn.2 <- subset(df.turn.2, df.turn.2$comparisons == "Same_site")
  
  
  Z = data.frame(value = df.turn.2$value, 
                 comp = rep("Z", nrow(df.turn.2)))
  
  df.turn.Y <- subset(df.turn, df.turn$comparisons == "Different_site" & df.turn$campain.x == "RUNARMS")
  
  Y = data.frame(value = df.turn.Y$value, 
                 comp = rep("Y", nrow(df.turn.Y)))
  
  df.turn.X <- subset(df.turn, df.turn$comparisons == "Different_site" & df.turn$campain.x == "P50ARMS")
    
  X = data.frame(value = df.turn.X$value, 
                 comp = rep("X", nrow(df.turn.X)))
  
  df <- rbind.data.frame(Z,Y,X)
  
  my_comparisons <- list( c("Z", "Y"), c("Y", "X"), c("Z", "X"))
  
  v1 <- ggboxplot(df, x = "comp", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "comp",
                  palette = c("black", "aquamarine3", "blue3"),
                  ylab = "Turnover") + 
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "anova")           
  
  v1 
  
  
    ### nest ####
  df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  df.nest <- subset(df.nest, row != col)
  
  # Convert factors to characters
  df.nest$row <- as.character(df.nest$row)
  df.nest$col <- as.character(df.nest$col)
  
  # Create a new column with sorted combinations
  df.nest$sorted_comparison <- apply(df.nest[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  
  # Identify and remove redundant comparisons
  df_unique <- df.nest %>% distinct(sorted_comparison, .keep_all = TRUE)
  
  # Remove the temporary column
  
  df.nest <- df_unique[, -ncol(df_unique)]
  
  
  merged_df <- merge(df.nest, meta_mean, by.x = c("row"), by.y = c("arms"))
  merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  df.nest <- subset(merged_df2, campain.x == campain.y)
  df.nest$comparisons <- ifelse(df.nest$site.x == df.nest$site.y, "Same_site", "Different_site")
  df.nest$same_site_campain <- paste0(df.nest$same_site, "_", df.nest$campain.x)
  
  
  ggplot(df.nest, aes(x=value)) + 
    geom_density() #Données normales

  p4 <- ggboxplot(df.nest, x = "campain.x", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "campain.x",
                  palette = c("blue4", "aquamarine3"),
                  facet.by = "comparisons", ylab = "Nestedness") + stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) 
  
  p4 
  
    ### jacc ####
  df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  df.jacc <- subset(df.jacc, row != col)
  
  # Convert factors to characters
  df.jacc$row <- as.character(df.jacc$row)
  df.jacc$col <- as.character(df.jacc$col)
  
  # Create a new column with sorted combinations
  df.jacc$sorted_comparison <- apply(df.jacc[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  
  # Identify and remove redundant comparisons
  df_unique <- df.jacc %>% distinct(sorted_comparison, .keep_all = TRUE)
  
  # Remove the temporary column
  
  df.jacc <- df_unique[, -ncol(df_unique)]
  
  
  merged_df <- merge(df.jacc, meta_mean, by.x = c("row"), by.y = c("arms"))
  merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  df.jacc <- subset(merged_df2, campain.x == campain.y)
  df.jacc$comparisons <- ifelse(df.jacc$site.x == df.jacc$site.y, "Same_site", "Different_site")
  df.jacc$same_site_campain <- paste0(df.jacc$same_site, "_", df.jacc$campain.x)
  
  
  ggplot(df.jacc, aes(x=value)) + 
    geom_density() #Données normales
  
  p5 <- ggboxplot(df.jacc, x = "campain.x", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "campain.x",
                  palette = c("blue4", "aquamarine3"),
                  facet.by = "comparisons", ylab = "jaccard dissimilarity") + stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) 
  
  
  p5 
  
  fin <- cowplot::plot_grid(p5, p3, p4, 
                            ncol = 1,
                            nrow = 3)
  
  path_to_boxplot_betadiv_site <- paste0("outputs/boxplot_beta_decomp_site.pdf")
  ggsave(filename =  path_to_boxplot_betadiv_site, plot = fin, width = 10, height = 13)
  
  
  return(path_to_boxplot_betadiv)
}