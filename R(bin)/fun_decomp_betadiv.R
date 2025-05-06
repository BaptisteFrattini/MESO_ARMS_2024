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
  library(ggplot2)
  library(forcats)
  library(ggsignif)
  library(ggpubr)
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean <- subset(data_mean, meta_mean$island == "Reunion")
  meta_mean <- subset(meta_mean, meta_mean$island == "Reunion")
  
  data_mean_pa <- vegan::decostand(data_mean, "pa") 
  
  B.pair.pa <- betapart::beta.pair(data_mean_pa, index.family = "jaccard")
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- B.pair.pa$beta.jac
  # #Comparisons at reunion scale ####
  #  ## Comparison between ARMS of the same depth and ARMS from =/= depth ####
  #   ### turn ####
  # 
  # df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  # df.turn <- subset(df.turn, row != col)
  # 
  # # Convert factors to characters
  # df.turn$row <- as.character(df.turn$row)
  # df.turn$col <- as.character(df.turn$col)
  # 
  # # Create a new column with sorted combinations
  # df.turn$sorted_comparison <- apply(df.turn[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  # 
  # # Identify and remove redundant comparisons
  # df_unique <- df.turn %>% distinct(sorted_comparison, .keep_all = TRUE)
  # 
  # # Remove the temporary column
  # df.turn <- df_unique[, -ncol(df_unique)]
  # 
  # merged_df <- merge(df.turn, meta_mean, by.x = c("row"), by.y = c("arms"))
  # merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  # 
  # df.turn <- data.frame(merged_df2$col, merged_df2$site.y, merged_df2$island.y, merged_df2$row, merged_df2$site.x, merged_df2$island.x, merged_df2$value) 
  # df.turn$same_site <- ifelse(df.turn$merged_df2.site.y == df.turn$merged_df2.site.x, "Yes", "No")
  # subset_df.turn <- subset(df.turn, same_site == "Yes")
  # subset_df.turn$row2 <- substr(subset_df.turn$merged_df2.row, 1, 7)
  # subset_df.turn$col2 <- substr(subset_df.turn$merged_df2.col, 1, 7)
  # subset_df.turn$same_depth <- ifelse(subset_df.turn$row2 == subset_df.turn$col2, "Yes", "No")
  # subset_df.turn <- data.frame(col = subset_df.turn$merged_df2.col, row = subset_df.turn$merged_df2.row, site = subset_df.turn$merged_df2.site.y, same_depth = subset_df.turn$same_depth, value = subset_df.turn$merged_df2.value)
  # 
  # turn.mean = tapply(subset_df.turn$value, list(subset_df.turn$same_depth, subset_df.turn$site), mean)
  # turn.sd = tapply(subset_df.turn$value, list(subset_df.turn$same_depth, subset_df.turn$site), sd)
  # 
 
  # 
  # 
  # intrainter = c("between ARMS \n of the same depth", "between mesophotic \n ARMS and shallow ARMS")
  # 
  # df_saint_leu <- subset(subset_df.turn, subset_df.turn$site == "Saint_Leu")
  # df_saint_leu <- dplyr::distinct(df_saint_leu, col, row, .keep_all = TRUE)
  # 
  # df_grand_bois <- subset(subset_df.turn, subset_df.turn$site == "Grand_Bois")
  # df_grand_bois <- dplyr::distinct(df_grand_bois, col, row, .keep_all = TRUE)
  # 
  # df_cap_houssaye <- subset(subset_df.turn, subset_df.turn$site == "Cap_La_Houssaye")
  # df_cap_houssaye <- dplyr::distinct(df_cap_houssaye, col, row, .keep_all = TRUE)
  # 
  # ggplot(df_saint_leu, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # 
  # p.sed <- rstatix::wilcox_test(df_saint_leu, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # a1 <- ggplot(df_saint_leu, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("coral","coral") ) +
  #   labs(title = NULL,
  #        x = "Comparisons",
  #        y = "Turnover component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels = NULL) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0, 0.9)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # a1
  # 
  # ggplot(df_grand_bois, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # shapiro.test(df_grand_bois$value) #Shapiro not OK 
  # 
  # p.sed <- rstatix::t_test(df_grand_bois, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # b1 <- ggplot(df_grand_bois, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("chartreuse4","chartreuse4") ) +
  #   labs(title = NULL,
  #        x = "Comparisons",
  #        y = "Turnover component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels = NULL) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0, 0.9)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # b1
  # 
  # ggplot(df_cap_houssaye, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # shapiro.test(df_cap_houssaye$value) #Shapiro not OK 
  # 
  # p.sed <- rstatix::wilcox_test(df_cap_houssaye, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # c1 <- ggplot(df_cap_houssaye, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("blue3","blue3") ) +
  #   labs(title = NULL,
  #        x = "Comparisons",
  #        y = "Turnover component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels = NULL) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0, 0.9)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # c1
  # 
  # 
  #   ### nest ####
  # 
  # df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  # df.nest <- subset(df.nest, row != col)
  # 
  # # Convert factors to characters
  # df.nest$row <- as.character(df.nest$row)
  # df.nest$col <- as.character(df.nest$col)
  # 
  # # Create a new column with sorted combinations
  # df.nest$sorted_comparison <- apply(df.nest[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  # 
  # # Identify and remove redundant comparisons
  # df_unique <- df.nest %>% distinct(sorted_comparison, .keep_all = TRUE)
  # 
  # # Remove the temporary column
  # df.nest <- df_unique[, -ncol(df_unique)]
  # 
  # 
  # merged_df <- merge(df.nest, meta_mean, by.x = c("row"), by.y = c("arms"))
  # merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  # 
  # df.nest <- data.frame(merged_df2$col, merged_df2$site.y, merged_df2$row, merged_df2$site.x, merged_df2$value) 
  # df.nest$same_site <- ifelse(df.nest$merged_df2.site.y == df.nest$merged_df2.site.x, "Yes", "No")
  # subset_df.nest <- subset(df.nest, same_site == "Yes")
  # subset_df.nest$row2 <- substr(subset_df.nest$merged_df2.row, 1, 7)
  # subset_df.nest$col2 <- substr(subset_df.nest$merged_df2.col, 1, 7)
  # subset_df.nest$same_depth <- ifelse(subset_df.nest$row2 == subset_df.nest$col2, "Yes", "No")
  # subset_df.nest <- data.frame(col = subset_df.nest$merged_df2.col, row = subset_df.nest$merged_df2.row, site = subset_df.nest$merged_df2.site.y, same_depth = subset_df.nest$same_depth, value = subset_df.nest$merged_df2.value)
  # 
  # nest.mean = tapply(subset_df.nest$value, list(subset_df.nest$same_depth, subset_df.nest$site), mean)
  # nest.sd = tapply(subset_df.nest$value, list(subset_df.nest$same_depth, subset_df.nest$site), sd)
  # 
  # intrainter = c("between ARMS \n of the same depth", "between mesophotic \n ARMS and shallow ARMS")
  # 
  # df_saint_leu <- subset(subset_df.nest, subset_df.nest$site == "Saint_Leu")
  # df_saint_leu <- dplyr::distinct(df_saint_leu, col, row, .keep_all = TRUE)
  # 
  # df_grand_bois <- subset(subset_df.nest, subset_df.nest$site == "Grand_Bois")
  # df_grand_bois <- dplyr::distinct(df_grand_bois, col, row, .keep_all = TRUE)
  # 
  # df_cap_houssaye <- subset(subset_df.nest, subset_df.nest$site == "Cap_La_Houssaye")
  # df_cap_houssaye <- dplyr::distinct(df_cap_houssaye, col, row, .keep_all = TRUE)
  # 
  # ggplot(df_saint_leu, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # shapiro.test(df_saint_leu$value) #Shapiro not OK 
  # 
  # p.sed <- rstatix::wilcox_test(df_saint_leu, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # a2 <- ggplot(df_saint_leu, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("coral","coral") ) +
  #   labs(title = NULL,
  #        x = "Comparisons",
  #        y = "Nestedness component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intrainter) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0, 0.25)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # a2
  # 
  # ggplot(df_grand_bois, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # shapiro.test(df_grand_bois$value) #Shapiro not OK 
  # 
  # p.sed <- rstatix::wilcox_test(df_grand_bois, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # b2 <- ggplot(df_grand_bois, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("chartreuse4","chartreuse4") ) +
  #   labs(title = NULL,
  #        x = "Comparisons",
  #        y = "Nestedness component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intrainter) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0, 0.25)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # b2
  # 
  # ggplot(df_cap_houssaye, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # shapiro.test(df_cap_houssaye$value) #Shapiro not OK 
  # 
  # p.sed <- rstatix::wilcox_test(df_cap_houssaye, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # c2 <- ggplot(df_cap_houssaye, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("blue3","blue3") ) +
  #   labs(title = NULL,
  #        x = "Comparisons",
  #        y = "Nestedness component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intrainter) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0, 0.25)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # c2
  # 
  #   ### jacc ####
  # 
  # df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  # df.jacc <- subset(df.jacc, row != col)
  # 
  # # Convert factors to characters
  # df.jacc$row <- as.character(df.jacc$row)
  # df.jacc$col <- as.character(df.jacc$col)
  # 
  # # Create a new column with sorted combinations
  # df.jacc$sorted_comparison <- apply(df.jacc[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  # 
  # # Identify and remove redundant comparisons
  # df_unique <- df.jacc %>% distinct(sorted_comparison, .keep_all = TRUE)
  # 
  # # Remove the temporary column
  # df.jacc <- df_unique[, -ncol(df_unique)]
  # 
  # 
  # 
  # merged_df <- merge(df.jacc, meta_mean, by.x = c("row"), by.y = c("arms"))
  # merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  # 
  # df.jacc <- data.frame(merged_df2$col, merged_df2$site.y, merged_df2$row, merged_df2$site.x, merged_df2$value) 
  # df.jacc$same_site <- ifelse(df.jacc$merged_df2.site.y == df.jacc$merged_df2.site.x, "Yes", "No")
  # subset_df.jacc <- subset(df.jacc, same_site == "Yes")
  # subset_df.jacc$row2 <- substr(subset_df.jacc$merged_df2.row, 1, 7)
  # subset_df.jacc$col2 <- substr(subset_df.jacc$merged_df2.col, 1, 7)
  # subset_df.jacc$same_depth <- ifelse(subset_df.jacc$row2 == subset_df.jacc$col2, "Yes", "No")
  # subset_df.jacc <- data.frame(col = subset_df.jacc$merged_df2.col, row = subset_df.jacc$merged_df2.row, site = subset_df.jacc$merged_df2.site.y, same_depth = subset_df.jacc$same_depth, value = subset_df.jacc$merged_df2.value)
  # 
  # jacc.mean = tapply(subset_df.jacc$value, list(subset_df.jacc$same_depth, subset_df.jacc$site), mean)
  # jacc.sd = tapply(subset_df.jacc$value, list(subset_df.jacc$same_depth, subset_df.jacc$site), sd)
  # 
  # intrainter = c("between ARMS \n of the same depth", "between mesophotic \n ARMS and shallow ARMS")
  # 
  # df_saint_leu <- subset(subset_df.jacc, subset_df.jacc$site == "Saint_Leu")
  # df_saint_leu <- dplyr::distinct(df_saint_leu, col, row, .keep_all = TRUE)
  # 
  # df_grand_bois <- subset(subset_df.jacc, subset_df.jacc$site == "Grand_Bois")
  # df_grand_bois <- dplyr::distinct(df_grand_bois, col, row, .keep_all = TRUE)
  # 
  # df_cap_houssaye <- subset(subset_df.jacc, subset_df.jacc$site == "Cap_La_Houssaye")
  # df_cap_houssaye <- dplyr::distinct(df_cap_houssaye, col, row, .keep_all = TRUE)
  # 
  # ggplot(df_saint_leu, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # shapiro.test(df_saint_leu$value) #Shapiro not OK 
  # 
  # p.sed <- rstatix::wilcox_test(df_saint_leu, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # a3 <- ggplot(df_saint_leu, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("coral","coral") ) +
  #   labs(title = "Saint Leu",
  #        x = "Comparisons",
  #        y = "Jaccard component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels = NULL) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0, 0.9)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # a3
  # 
  # ggplot(df_grand_bois, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # shapiro.test(df_grand_bois$value) #Shapiro not OK 
  # 
  # p.sed <- rstatix::wilcox_test(df_grand_bois, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # b3 <- ggplot(df_grand_bois, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("chartreuse4","chartreuse4") ) +
  #   labs(title = "Grand Bois",
  #        x = "Comparisons",
  #        y = "Jaccard component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels = NULL) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0,0.9)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # b3
  # 
  # ggplot(df_cap_houssaye, aes(x=value)) + 
  #   geom_density() #Données pas normales normales
  # shapiro.test(df_cap_houssaye$value) #Shapiro not OK 
  # 
  # p.sed <- rstatix::t_test(df_cap_houssaye, value ~ same_depth) #parametrique
  # p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.3)
  # c3 <- ggplot(df_cap_houssaye, aes(x = fct_relevel(same_depth, "Yes", "No"), y = value)) +
  #   geom_boxplot(fill =  c("blue3","blue3") ) +
  #   labs(title = "Cap La Houssaye",
  #        x = "Comparisons",
  #        y = "Jaccard component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=NULL) +
  #   theme_classic() + 
  #   scale_y_continuous(limits = c(0, 0.9)) +
  #   stat_pvalue_manual(p.sed) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=9))
  # c3
  # 
  # fin <- cowplot::plot_grid(c3, a3, b3, c1, a1, b1, c2, a2, b2, 
  #                           ncol = 3,
  #                           nrow = 3)
  # 
  # path_to_boxplot_betadiv <- paste0("outputs/boxplot_beta_decomp.pdf")
  # ggsave(filename =  path_to_boxplot_betadiv, plot = fin, width = 13, height = 15.5)
  # 
  # 
  # 
  #   ## Comparison between ARMS of differents site and ARMS of different depth ####
  #     ### turn ####
  # df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  # df.turn <- subset(df.turn, row != col)
  # 
  # # Convert factors to characters
  # df.turn$row <- as.character(df.turn$row)
  # df.turn$col <- as.character(df.turn$col)
  # 
  # # Create a new column with sorted combinations
  # df.turn$sorted_comparison <- apply(df.turn[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  # 
  # # Identify and remove redundant comparisons
  # df_unique <- df.turn %>% distinct(sorted_comparison, .keep_all = TRUE)
  # 
  # # Remove the temporary column
  # 
  # df.turn <- df_unique[, -ncol(df_unique)]
  # 
  # 
  # merged_df <- merge(df.turn, meta_mean, by.x = c("row"), by.y = c("arms"))
  # merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  # df.turn <- subset(merged_df2, campain.x == campain.y)
  # df.turn$comparisons <- ifelse(df.turn$site.x == df.turn$site.y, "Same_site", "Different_site")
  # 
  # 
  # 
  # ggplot(df.turn, aes(x=value)) + 
  #   geom_density() #Données normales
  # 
  # 

  # p3 <- ggboxplot(df.turn, x = "campain.x", 
  #                 y = "value",
  #                 add = "jitter", 
  #                 short.panel.labs = FALSE,
  #                 color = "campain.x",
  #                 palette = c("blue4", "aquamarine3"),
  #                 facet.by = "comparisons", ylab = "Turnover")+ stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) 
  # 
  # p3 
  # 
  # ### trying to include comparison between deep and shallow
  # df.turn.2 <- subset(merged_df2, campain.x != campain.y)
  # df.turn.2$comparisons <- ifelse(df.turn.2$site.x == df.turn.2$site.y, "Same_site", "Different_site")
  # df.turn.2 <- subset(df.turn.2, df.turn.2$comparisons == "Same_site")
  # 
  # 
  # Z = data.frame(value = df.turn.2$value, 
  #                comp = rep("Z", nrow(df.turn.2)))
  # 
  # df.turn.Y <- subset(df.turn, df.turn$comparisons == "Different_site" & df.turn$campain.x == "RUNARMS")
  # 
  # Y = data.frame(value = df.turn.Y$value, 
  #                comp = rep("Y", nrow(df.turn.Y)))
  # 
  # df.turn.X <- subset(df.turn, df.turn$comparisons == "Different_site" & df.turn$campain.x == "P50ARMS")
  #   
  # X = data.frame(value = df.turn.X$value, 
  #                comp = rep("X", nrow(df.turn.X)))
  # 
  # df <- rbind.data.frame(Z,Y,X)
  # 
  # my_comparisons <- list( c("Z", "Y"), c("Y", "X"), c("Z", "X"))
  # 
  # v1 <- ggboxplot(df, x = "comp", 
  #                 y = "value",
  #                 add = "jitter", 
  #                 short.panel.labs = FALSE,
  #                 color = "comp",
  #                 palette = c("black", "aquamarine3", "blue3"),
  #                 ylab = "Turnover") + 
  #   stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  #   stat_compare_means(label.y = 1, method = "anova") + theme_classic()          
  # 
  # v1 
  # 
  # 
  #   ### nest ####
  # df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  # df.nest <- subset(df.nest, row != col)
  # 
  # # Convert factors to characters
  # df.nest$row <- as.character(df.nest$row)
  # df.nest$col <- as.character(df.nest$col)
  # 
  # # Create a new column with sorted combinations
  # df.nest$sorted_comparison <- apply(df.nest[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  # 
  # # Identify and remove redundant comparisons
  # df_unique <- df.nest %>% distinct(sorted_comparison, .keep_all = TRUE)
  # 
  # # Remove the temporary column
  # 
  # df.nest <- df_unique[, -ncol(df_unique)]
  # 
  # 
  # merged_df <- merge(df.nest, meta_mean, by.x = c("row"), by.y = c("arms"))
  # merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  # df.nest <- subset(merged_df2, campain.x == campain.y)
  # df.nest$comparisons <- ifelse(df.nest$site.x == df.nest$site.y, "Same_site", "Different_site")
  # df.nest$same_site_campain <- paste0(df.nest$same_site, "_", df.nest$campain.x)
  # 
  # 
  # ggplot(df.nest, aes(x=value)) + 
  #   geom_density() #Données normales
  # 
  # p4 <- ggboxplot(df.nest, x = "campain.x", 
  #                 y = "value",
  #                 add = "jitter", 
  #                 short.panel.labs = FALSE,
  #                 color = "campain.x",
  #                 palette = c("blue4", "aquamarine3"),
  #                 facet.by = "comparisons", ylab = "Nestedness") + stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) 
  # 
  # p4 
  # 
  # ### trying to include comparison between deep and shallow
  # df.nest.2 <- subset(merged_df2, campain.x != campain.y)
  # df.nest.2$comparisons <- ifelse(df.nest.2$site.x == df.nest.2$site.y, "Same_site", "Different_site")
  # df.nest.2 <- subset(df.nest.2, df.nest.2$comparisons == "Same_site")
  # 
  # 
  # Z = data.frame(value = df.nest.2$value, 
  #                comp = rep("Z", nrow(df.nest.2)))
  # 
  # df.nest.Y <- subset(df.nest, df.nest$comparisons == "Different_site" & df.nest$campain.x == "RUNARMS")
  # 
  # Y = data.frame(value = df.nest.Y$value, 
  #                comp = rep("Y", nrow(df.nest.Y)))
  # 
  # df.nest.X <- subset(df.nest, df.nest$comparisons == "Different_site" & df.nest$campain.x == "P50ARMS")
  # 
  # X = data.frame(value = df.nest.X$value, 
  #                comp = rep("X", nrow(df.nest.X)))
  # 
  # df <- rbind.data.frame(Z,Y,X)
  # 
  # my_comparisons <- list( c("Z", "Y"), c("Y", "X"), c("Z", "X"))
  # 
  # v2 <- ggboxplot(df, x = "comp", 
  #                 y = "value",
  #                 add = "jitter", 
  #                 short.panel.labs = FALSE,
  #                 color = "comp",
  #                 palette = c("black", "aquamarine3", "blue3"),
  #                 ylab = "Nestedness") + 
  #   stat_compare_means(label.y = 1, method = "kruskal.test")  + theme_classic()         
  # 
  # v2 
  # 
  # 
  #   ### jacc ####
  # df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  # df.jacc <- subset(df.jacc, row != col)
  # 
  # # Convert factors to characters
  # df.jacc$row <- as.character(df.jacc$row)
  # df.jacc$col <- as.character(df.jacc$col)
  # 
  # # Create a new column with sorted combinations
  # df.jacc$sorted_comparison <- apply(df.jacc[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  # 
  # # Identify and remove redundant comparisons
  # df_unique <- df.jacc %>% distinct(sorted_comparison, .keep_all = TRUE)
  # 
  # # Remove the temporary column
  # 
  # df.jacc <- df_unique[, -ncol(df_unique)]
  # 
  # 
  # merged_df <- merge(df.jacc, meta_mean, by.x = c("row"), by.y = c("arms"))
  # merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  # df.jacc <- subset(merged_df2, campain.x == campain.y)
  # df.jacc$comparisons <- ifelse(df.jacc$site.x == df.jacc$site.y, "Same_site", "Different_site")
  # df.jacc$same_site_campain <- paste0(df.jacc$same_site, "_", df.jacc$campain.x)
  # 
  # 
  # ggplot(df.jacc, aes(x=value)) + 
  #   geom_density() #Données normales
  # 
  # p5 <- ggboxplot(df.jacc, x = "campain.x", 
  #                 y = "value",
  #                 add = "jitter", 
  #                 short.panel.labs = FALSE,
  #                 color = "campain.x",
  #                 palette = c("blue4", "aquamarine3"),
  #                 facet.by = "comparisons", ylab = "jaccard dissimilarity") + stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) 
  # 
  # 
  # p5 
  # 
  # ### trying to include comparison between deep and shallow
  # df.jacc.2 <- subset(merged_df2, campain.x != campain.y)
  # df.jacc.2$comparisons <- ifelse(df.jacc.2$site.x == df.jacc.2$site.y, "Same_site", "Different_site")
  # df.jacc.2 <- subset(df.jacc.2, df.jacc.2$comparisons == "Same_site")
  # 
  # 
  # Z = data.frame(value = df.jacc.2$value, 
  #                comp = rep("Z", nrow(df.jacc.2)))
  # 
  # df.jacc.Y <- subset(df.jacc, df.jacc$comparisons == "Different_site" & df.jacc$campain.x == "RUNARMS")
  # 
  # Y = data.frame(value = df.jacc.Y$value, 
  #                comp = rep("Y", nrow(df.jacc.Y)))
  # 
  # df.jacc.X <- subset(df.jacc, df.jacc$comparisons == "Different_site" & df.jacc$campain.x == "P50ARMS")
  # 
  # X = data.frame(value = df.jacc.X$value, 
  #                comp = rep("X", nrow(df.jacc.X)))
  # 
  # df <- rbind.data.frame(Z,Y,X)
  # 
  # my_comparisons <- list( c("Z", "Y"), c("Y", "X"), c("Z", "X"))
  # 
  # v3 <- ggboxplot(df, x = "comp", 
  #                 y = "value",
  #                 add = "jitter", 
  #                 short.panel.labs = FALSE,
  #                 color = "comp",
  #                 palette = c("black", "aquamarine3", "blue3"),
  #                 ylab = "jaccard dissimilarity") + 
  #   stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
  #   stat_compare_means(label.y = 1, method = "anova") + theme_classic()          
  # 
  # v3 
  # 
  # fin <- cowplot::plot_grid(v3, v1, v2, 
  #                           ncol = 1,
  #                           nrow = 3)
  # 
  # path_to_boxplot_betadiv_XYZ <- paste0("outputs/boxplot_beta_decomp_XYZ.pdf")
  # ggsave(filename =  path_to_boxplot_betadiv_XYZ, plot = fin, width = 6, height = 11.5)
  # 
  # 
  # 
  # fin <- cowplot::plot_grid(p5, p3, p4, 
  #                           ncol = 1,
  #                           nrow = 3)
  # 
  # path_to_boxplot_betadiv_site <- paste0("outputs/boxplot_beta_decomp_site.pdf")
  # ggsave(filename =  path_to_boxplot_betadiv_site, plot = fin, width = 10, height = 13)
  # 
  # Comparison at mascarene scale ####
  
  data_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  data_mean_pa <- vegan::decostand(data_mean, "pa") 
  
  B.pair.pa <- betapart::beta.pair(data_mean_pa, index.family = "jaccard")
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- B.pair.pa$beta.jac
  mat.bray <- vegan::vegdist(data_mean, "bray")
  
  #### bray ####
  
  df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
  df.bray <- subset(df.bray, row != col)
  
  # Convert factors to characters
  df.bray$row <- as.character(df.bray$row)
  df.bray$col <- as.character(df.bray$col)
  
  # Create a new column with sorted combinations
  df.bray$sorted_comparison <- apply(df.bray[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))
  
  # Identify and remove redundant comparisons
  df_unique <- df.bray %>% distinct(sorted_comparison, .keep_all = TRUE)
  
  # Remove the temporary column
  df.bray <- df_unique[, -ncol(df_unique)]
  
  merged_df <- merge(df.bray, meta_mean, by.x = c("row"), by.y = c("arms"))
  merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))
  
  df.bray.V <- subset(merged_df2,  campain.x == "RODARMS" & campain.y == "RODARMS")
  df.bray.V$comparisons <- ifelse(df.bray.V$site.x == df.bray.V$site.y, "Same_site", "Different_site")
  df.bray.V <- subset(df.bray.V, df.bray.V$comparisons == "Different_site")
  
  V = data.frame(value = df.bray.V$value, 
                 comp = rep("V", nrow(df.bray.V)))
  
  df.bray.W <- subset(merged_df2, campain.x != "P50ARMS" & campain.y != "P50ARMS")
  df.bray.W$comparisons <- ifelse(df.bray.W$campain.x == df.bray.W$campain.y, "Same_Island", "Different_Island")
  df.bray.W <- subset(df.bray.W, df.bray.W$comparisons == "Different_Island")
  
  
  W = data.frame(value = df.bray.W$value, 
                 comp = rep("W", nrow(df.bray.W)))
  
  df.bray.Z <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.bray.Z <- subset(df.bray.Z, campain.x != campain.y)
  df.bray.Z$comparisons <- ifelse(df.bray.Z$site.x == df.bray.Z$site.y, "Same_site", "Different_site")
  df.bray.Z <- subset(df.bray.Z, df.bray.Z$comparisons == "Same_site")
  
  Z = data.frame(value = df.bray.Z$value, 
                 comp = rep("Z", nrow(df.bray.Z)))
  
  df.bray.Y <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.bray.Y <- subset(df.bray.Y, campain.x == campain.y)
  df.bray.Y$comparisons <- ifelse(df.bray.Y$site.x == df.bray.Y$site.y, "Same_site", "Different_site")
  df.bray.Y <- subset(df.bray.Y, df.bray.Y$comparisons == "Different_site" & df.bray.Y$campain.x == "RUNARMS")
  
  Y = data.frame(value = df.bray.Y$value, 
                 comp = rep("Y", nrow(df.bray.Y)))
  
  df.bray.X <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.bray.X <- subset(df.bray.X, campain.x == campain.y)
  df.bray.X$comparisons <- ifelse(df.bray.X$site.x == df.bray.X$site.y, "Same_site", "Different_site")
  df.bray.X <- subset(df.bray.X, df.bray.X$comparisons == "Different_site" & df.bray.X$campain.x == "P50ARMS")
  
  X = data.frame(value = df.bray.X$value, 
                 comp = rep("X", nrow(df.bray.X)))
  
  df <- rbind(V,W,X,Y,Z)
  
  my_comparisons <- list( c("V", "W"), c("V", "X"), c("V", "Y"), c("V","Z"), c("W","X"), c("W", "Y"), c("W","Z"), c("X","Y"), c("X","Z"), c("Y","Z"))
  
  
  library(ggplot2)
  
  df$comp <- factor(df$comp, levels = c("X", "Y", "Z", "W", "V"))
  
  q1 <- ggboxplot(df, 
                  x = "comp", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "comp",
                  palette = c("purple4" , "darkgreen", "blue3","aquamarine3", "black"),
                  ylab = "Bray-Curtis") + 
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "anova") + theme_classic()          
  
  
  q1 
  
  # Calcul des comparaisons post-hoc avec t.test
  comp_results <- compare_means(value ~ comp, data = df, method = "t.test", p.adjust.method = "bonferroni")
  
  # Extraction des p-values ajustées
  pvals <- comp_results$p.adj
  names(pvals) <- apply(comp_results[c("group1", "group2")], 1, paste, collapse = "-")
  
  # Génération des lettres
  
  letters <- multcompView::multcompLetters(pvals)$Letters
  comp_means <- df %>% group_by(comp) %>% summarize(mean = mean(value))
  
  # Ajouter les lettres aux moyennes
  comp_means$letters <- letters[as.character(comp_means$comp)]
  
  
  q1bis <- ggboxplot(df, x = "comp", 
                     y = "value",
                     add = "jitter", 
                     short.panel.labs = FALSE,
                     color = "comp",
                     palette = c("blue3", "aquamarine3", "black",  "darkgreen","purple4"),
                     ylab = "Bray-Curtis") +  # Change x-axis label 
    scale_x_discrete(name = NULL, labels = NULL)+
    theme(legend.position = "none")+
    ylim(0, 0.8)
  
  
  # Ajouter les lettres pour montrer les différences significatives
  q1bis <- q1bis + geom_text(data = comp_means, aes(x = comp, y = 0.75, label = letters),
                             color = "black", vjust = 0)

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

  df.turn.V <- subset(merged_df2,  campain.x == "RODARMS" & campain.y == "RODARMS")
  df.turn.V$comparisons <- ifelse(df.turn.V$site.x == df.turn.V$site.y, "Same_site", "Different_site")
  df.turn.V <- subset(df.turn.V, df.turn.V$comparisons == "Different_site")
  
  V = data.frame(value = df.turn.V$value, 
                 comp = rep("V", nrow(df.turn.V)))
  
  df.turn.W <- subset(merged_df2, campain.x != "P50ARMS" & campain.y != "P50ARMS")
  df.turn.W$comparisons <- ifelse(df.turn.W$campain.x == df.turn.W$campain.y, "Same_Island", "Different_Island")
  df.turn.W <- subset(df.turn.W, df.turn.W$comparisons == "Different_Island")
  
  
  W = data.frame(value = df.turn.W$value, 
                 comp = rep("W", nrow(df.turn.W)))
  
  df.turn.Z <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.turn.Z <- subset(df.turn.Z, campain.x != campain.y)
  df.turn.Z$comparisons <- ifelse(df.turn.Z$site.x == df.turn.Z$site.y, "Same_site", "Different_site")
  df.turn.Z <- subset(df.turn.Z, df.turn.Z$comparisons == "Same_site")
  
  Z = data.frame(value = df.turn.Z$value, 
                 comp = rep("Z", nrow(df.turn.Z)))
  
  df.turn.Y <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.turn.Y <- subset(df.turn.Y, campain.x == campain.y)
  df.turn.Y$comparisons <- ifelse(df.turn.Y$site.x == df.turn.Y$site.y, "Same_site", "Different_site")
  df.turn.Y <- subset(df.turn.Y, df.turn.Y$comparisons == "Different_site" & df.turn.Y$campain.x == "RUNARMS")
  
  Y = data.frame(value = df.turn.Y$value, 
                 comp = rep("Y", nrow(df.turn.Y)))
  
  df.turn.X <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.turn.X <- subset(df.turn.X, campain.x == campain.y)
  df.turn.X$comparisons <- ifelse(df.turn.X$site.x == df.turn.X$site.y, "Same_site", "Different_site")
  df.turn.X <- subset(df.turn.X, df.turn.X$comparisons == "Different_site" & df.turn.X$campain.x == "P50ARMS")
  
  X = data.frame(value = df.turn.X$value, 
                 comp = rep("X", nrow(df.turn.X)))
  
  df <- rbind(V,W,X,Y,Z)
  
  my_comparisons <- list( c("V", "W"), c("V", "X"), c("V", "Y"), c("V","Z"), c("W","X"), c("W", "Y"), c("W","Z"), c("X","Y"), c("X","Z"), c("Y","Z"))
  

  library(ggplot2)
  
  df$comp <- factor(df$comp, levels = c("X", "Y", "Z", "W", "V"))
  
  r1 <- ggboxplot(df, 
                  x = "comp", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "comp",
                  palette = c("purple4" , "darkgreen", "blue3","aquamarine3", "black"),
                  ylab = "Turnover") + 
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "anova") + theme_classic()          
  
  
  r1 
  
  # Calcul des comparaisons post-hoc avec t.test
  comp_results <- compare_means(value ~ comp, data = df, method = "t.test", p.adjust.method = "bonferroni")
  
  # Extraction des p-values ajustées
  pvals <- comp_results$p.adj
  names(pvals) <- apply(comp_results[c("group1", "group2")], 1, paste, collapse = "-")
  
  # Génération des lettres
 
  letters <- multcompView::multcompLetters(pvals)$Letters
  comp_means <- df %>% group_by(comp) %>% summarize(mean = mean(value))
  
  # Ajouter les lettres aux moyennes
  comp_means$letters <- letters[as.character(comp_means$comp)]
  
  new_labels <- c("V" = paste0("Between RODARMS \n sites (Shallow) - N = ", nrow(V)),
                  "W" = paste0("Between Islands \n (Shallow) - N = ", nrow(W)),
                  "X" = paste0("Between RUNARMS \n sites (Mesophotic) - N = ", nrow(X)),
                  "Y" = paste0("Between RUNARMS \n sites (Shallow) - N = ", nrow(Y)), 
                  "Z" = paste0("Between depth within \n sites (RUNARMS) - N = ", nrow(Z)))
  
  r1bis <- ggboxplot(df, x = "comp", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "comp",
                  palette = c("blue3", "aquamarine3", "black",  "darkgreen","purple4"),
                  ylab = "Turnover") + 
    scale_x_discrete(name = "Comparisons", labels = new_labels) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), legend.position = "none")+
    ylim(0, 0.8)
    
  
  # Ajouter les lettres pour montrer les différences significatives
  r1bis <- r1bis + geom_text(data = comp_means, aes(x = comp, y = 0.75, label = letters),
                       color = "black", vjust = 0)
  
  r1bis_2 <- ggboxplot(df, x = "comp", 
                     y = "value",
                     add = "jitter", 
                     short.panel.labs = FALSE,
                     color = "comp",
                     palette = c("blue3", "aquamarine3", "black",  "darkgreen","purple4"),
                     ylab = "Turnover") + 
    scale_x_discrete(name = NULL, labels = NULL) +
    theme(legend.position = "none")+
    ylim(0, 0.8)
  
  
  # Ajouter les lettres pour montrer les différences significatives
  r1bis_2 <- r1bis_2 + geom_text(data = comp_means, aes(x = comp, y = 0.75, label = letters),
                             color = "black", vjust = 0)
  
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
  
  df.jacc.V <- subset(merged_df2,  campain.x == "RODARMS" & campain.y == "RODARMS")
  df.jacc.V$comparisons <- ifelse(df.jacc.V$site.x == df.jacc.V$site.y, "Same_site", "Different_site")
  df.jacc.V <- subset(df.jacc.V, df.jacc.V$comparisons == "Different_site")
  
  V = data.frame(value = df.jacc.V$value, 
                 comp = rep("V", nrow(df.jacc.V)))
  
  df.jacc.W <- subset(merged_df2, campain.x != "P50ARMS" & campain.y != "P50ARMS")
  df.jacc.W$comparisons <- ifelse(df.jacc.W$campain.x == df.jacc.W$campain.y, "Same_Island", "Different_Island")
  df.jacc.W <- subset(df.jacc.W, df.jacc.W$comparisons == "Different_Island")
  
  
  W = data.frame(value = df.jacc.W$value, 
                 comp = rep("W", nrow(df.jacc.W)))
  
  df.jacc.Z <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.jacc.Z <- subset(df.jacc.Z, campain.x != campain.y)
  df.jacc.Z$comparisons <- ifelse(df.jacc.Z$site.x == df.jacc.Z$site.y, "Same_site", "Different_site")
  df.jacc.Z <- subset(df.jacc.Z, df.jacc.Z$comparisons == "Same_site")
  
  Z = data.frame(value = df.jacc.Z$value, 
                 comp = rep("Z", nrow(df.jacc.Z)))
  
  df.jacc.Y <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.jacc.Y <- subset(df.jacc.Y, campain.x == campain.y)
  df.jacc.Y$comparisons <- ifelse(df.jacc.Y$site.x == df.jacc.Y$site.y, "Same_site", "Different_site")
  df.jacc.Y <- subset(df.jacc.Y, df.jacc.Y$comparisons == "Different_site" & df.jacc.Y$campain.x == "RUNARMS")
  
  Y = data.frame(value = df.jacc.Y$value, 
                 comp = rep("Y", nrow(df.jacc.Y)))
  
  df.jacc.X <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.jacc.X <- subset(df.jacc.X, campain.x == campain.y)
  df.jacc.X$comparisons <- ifelse(df.jacc.X$site.x == df.jacc.X$site.y, "Same_site", "Different_site")
  df.jacc.X <- subset(df.jacc.X, df.jacc.X$comparisons == "Different_site" & df.jacc.X$campain.x == "P50ARMS")
  
  X = data.frame(value = df.jacc.X$value, 
                 comp = rep("X", nrow(df.jacc.X)))
  
  df <- rbind(V,W,X,Y,Z)
  
  df$comp <- factor(df$comp, levels = c("X", "Y", "Z", "W", "V"))
  
  my_comparisons <- list( c("V", "W"), c("V", "X"), c("V", "Y"), c("V","Z"), c("W","X"), c("W", "Y"), c("W","Z"), c("X","Y"), c("X","Z"), c("Y","Z"))
  library(ggplot2)
  s1 <- ggboxplot(df, x = "comp", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "comp",
                  palette = c("blue3", "aquamarine3", "black",  "darkgreen","purple4"),
                  ylab = "jaccard") + 
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "anova") + theme_classic()          
  
  
  s1 
  
  # Calcul des comparaisons post-hoc avec t.test
  comp_results <- compare_means(value ~ comp, data = df, method = "t.test", p.adjust.method = "bonferroni")
  
  # Extraction des p-values ajustées
  pvals <- comp_results$p.adj
  names(pvals) <- apply(comp_results[c("group1", "group2")], 1, paste, collapse = "-")
  
  # Génération des lettres
  
  letters <- multcompView::multcompLetters(pvals)$Letters
  comp_means <- df %>% group_by(comp) %>% summarize(mean = mean(value))
  
  # Ajouter les lettres aux moyennes
  comp_means$letters <- letters[as.character(comp_means$comp)]
  
  s1bis <- ggboxplot(df, x = "comp", 
                     y = "value",
                     add = "jitter", 
                     short.panel.labs = FALSE,
                     color = "comp",
                     palette = c("blue3", "aquamarine3", "black",  "darkgreen","purple4"),
                     ylab = "Jaccard") +  # Change x-axis label 
    scale_x_discrete(name = NULL, labels = NULL) +
    theme(legend.position = "none")+
    ylim(0, 0.8)
  
  # Ajouter les lettres pour montrer les différences significatives
  s1bis <- s1bis + geom_text(data = comp_means, aes(x = comp, y = 0.78, label = letters),
                             color = "black", vjust = 0)
  
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
  
  df.nest.V <- subset(merged_df2,  campain.x == "RODARMS" & campain.y == "RODARMS")
  df.nest.V$comparisons <- ifelse(df.nest.V$site.x == df.nest.V$site.y, "Same_site", "Different_site")
  df.nest.V <- subset(df.nest.V, df.nest.V$comparisons == "Different_site")
  
  V = data.frame(value = df.nest.V$value, 
                 comp = rep("V", nrow(df.nest.V)))
  
  df.nest.W <- subset(merged_df2, campain.x != "P50ARMS" & campain.y != "P50ARMS")
  df.nest.W$comparisons <- ifelse(df.nest.W$campain.x == df.nest.W$campain.y, "Same_Island", "Different_Island")
  df.nest.W <- subset(df.nest.W, df.nest.W$comparisons == "Different_Island")
  
  
  W = data.frame(value = df.nest.W$value, 
                 comp = rep("W", nrow(df.nest.W)))
  
  df.nest.Z <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.nest.Z <- subset(df.nest.Z, campain.x != campain.y)
  df.nest.Z$comparisons <- ifelse(df.nest.Z$site.x == df.nest.Z$site.y, "Same_site", "Different_site")
  df.nest.Z <- subset(df.nest.Z, df.nest.Z$comparisons == "Same_site")
  
  Z = data.frame(value = df.nest.Z$value, 
                 comp = rep("Z", nrow(df.nest.Z)))
  
  df.nest.Y <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.nest.Y <- subset(df.nest.Y, campain.x == campain.y)
  df.nest.Y$comparisons <- ifelse(df.nest.Y$site.x == df.nest.Y$site.y, "Same_site", "Different_site")
  df.nest.Y <- subset(df.nest.Y, df.nest.Y$comparisons == "Different_site" & df.nest.Y$campain.x == "RUNARMS")
  
  Y = data.frame(value = df.nest.Y$value, 
                 comp = rep("Y", nrow(df.nest.Y)))
  
  df.nest.X <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
  df.nest.X <- subset(df.nest.X, campain.x == campain.y)
  df.nest.X$comparisons <- ifelse(df.nest.X$site.x == df.nest.X$site.y, "Same_site", "Different_site")
  df.nest.X <- subset(df.nest.X, df.nest.X$comparisons == "Different_site" & df.nest.X$campain.x == "P50ARMS")
  
  X = data.frame(value = df.nest.X$value, 
                 comp = rep("X", nrow(df.nest.X)))
  
  df <- rbind(V,W,X,Y,Z)
  
  new_labels <- c("V" = paste0("Between RODARMS \n sites (Shallow) - N = ", nrow(V)),
                  "W" = paste0("Between Islands \n (Shallow) - N = ", nrow(W)),
                  "X" = paste0("Between RUNARMS \n sites (Mesophotic) - N = ", nrow(X)),
                  "Y" = paste0("Between RUNARMS \n sites (Shallow) - N = ", nrow(Y)), 
                  "Z" = paste0("Between depth within \n sites (RUNARMS) - N = ", nrow(Z)))
  
  
  my_comparisons <- list( c("V", "W"), c("V", "X"), c("V", "Y"), c("V","Z"), c("W","X"), c("W", "Y"), c("W","Z"), c("X","Y"), c("X","Z"), c("Y","Z"))
  library(ggplot2)
  
  df$comp <- factor(df$comp, levels = c("X", "Y", "Z", "W", "V"))
  
  t1 <- ggboxplot(df, x = "comp", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "comp",
                  palette =c("blue3", "aquamarine3", "black",  "darkgreen","purple4"),
                  ylab = "Nestedness") + 
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "anova") + theme_classic()          
  
  
  t1 
  
  # Calcul des comparaisons post-hoc avec t.test
  comp_results <- compare_means(value ~ comp, data = df, method = "t.test", p.adjust.method = "bonferroni")
  
  # Extraction des p-values ajustées
  pvals <- comp_results$p.adj
  names(pvals) <- apply(comp_results[c("group1", "group2")], 1, paste, collapse = "-")
  
  # Génération des lettres
  
  letters <- multcompView::multcompLetters(pvals)$Letters
  comp_means <- df %>% group_by(comp) %>% summarize(mean = mean(value))
  
  # Ajouter les lettres aux moyennes
  comp_means$letters <- letters[as.character(comp_means$comp)]
  
  t1bis <- ggboxplot(df, x = "comp", 
                     y = "value",
                     add = "jitter", 
                     short.panel.labs = FALSE,
                     color = "comp",
                     palette = c("blue3", "aquamarine3", "black",  "darkgreen","purple4"),
                     ylab = "Nestedness") + 
    scale_x_discrete(name = "Comparisons", labels = new_labels) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), legend.position = "none")+
    ylim(0, 0.8)

  
  
  # Ajouter les lettres pour montrer les différences significatives
  t1bis <- t1bis + geom_text(data = comp_means, aes(x = comp, y = 0.25, label = letters),
                             color = "black", vjust = 0)
 
  t1bis_2 <- ggboxplot(df, x = "comp", 
                     y = "value",
                     add = "jitter", 
                     short.panel.labs = FALSE,
                     color = "comp",
                     palette = c("blue3", "aquamarine3", "black",  "darkgreen","purple4"),
                     ylab = "Nestedness") + 
    scale_x_discrete(name = NULL, labels = NULL) +
    theme(legend.position = "none") +
    ylim(0, 0.8)
  
  
  # Ajouter les lettres pour montrer les différences significatives
  t1bis_2 <- t1bis_2 + geom_text(data = comp_means, aes(x = comp, y = 0.25, label = letters),
                             color = "black", vjust = 0)
  
  
  fin_1 <- cowplot::plot_grid(q1bis,  r1bis_2, s1bis, t1bis_2, 
                            ncol = 2,
                            nrow = 2)
  
  fin_2 <- cowplot::plot_grid(q1bis,s1bis, r1bis, t1bis,
                            ncol = 2,
                            nrow = 2)
  
  path_to_boxplot_betadiv_XYZ_bis_2 <- paste0("outputs/boxplot_beta_decomp_XYZ_bis_2.pdf")
  ggsave(filename =  path_to_boxplot_betadiv_XYZ_bis_2, plot = fin_1, width = 10, height = 9)
  
  path_to_boxplot_betadiv_XYZ_bis <- paste0("outputs/boxplot_beta_decomp_XYZ_bis.pdf")
  ggsave(filename =  path_to_boxplot_betadiv_XYZ_bis, plot = fin_2, width = 10, height = 9)
  
  return(path_to_boxplot_betadiv_XYZ_bis_2)
}






# df.turn <- data.frame(merged_df2$col, merged_df2$site.y, merged_df2$island.y, merged_df2$row, merged_df2$site.x, merged_df2$island.x, merged_df2$value) 
# df.turn$same_site <- ifelse(df.turn$merged_df2.site.y == df.turn$merged_df2.site.x, "Yes", "No")
# df.turn$same_island <- ifelse(df.turn$merged_df2.island.y == df.turn$merged_df2.island.x, "Yes", "No")
# subset_df.turn <- subset(df.turn, same_site == "Yes")
# subset_df.turn$row2 <- substr(subset_df.turn$merged_df2.row, 1, 7)
# subset_df.turn$col2 <- substr(subset_df.turn$merged_df2.col, 1, 7)
# subset_df.turn$same_depth <- ifelse(subset_df.turn$row2 == subset_df.turn$col2, "Yes", "No")
# subset_df.turn <- data.frame(col = subset_df.turn$merged_df2.col, row = subset_df.turn$merged_df2.row, site = subset_df.turn$merged_df2.site.y, same_depth = subset_df.turn$same_depth, value = subset_df.turn$merged_df2.value)
