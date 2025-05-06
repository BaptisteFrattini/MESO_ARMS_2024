#' Beta diversity decomposing
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_decomp_betadiv_fullsites <- function(data_and_meta_clean){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  library(betapart)
  library(reshape2)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(ggsignif)
  library(ggpubr)
  library(car)
  library(ggplot2)
  
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  
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
  
  
  # ANOVA
  anova_res <- aov(value ~ comp, data = df)
  summary(anova_res)
  
  # Extraire les résidus
  residus <- residuals(anova_res)
  
  # Test de Shapiro-Wilk
  shapiro.test(residus)
  
  # Graphiques diagnostiques de l'ANOVA
  par(mfrow = c(2, 2))  # Affiche 4 graphiques à la suite
  plot(anova_res)
  
  
  bartlett.test(residus ~ comp, data = df)
  
  df$comp <- factor(df$comp, levels = c("X", "Y", "Z", "W", "V"))
  
  q1 <- ggboxplot(df, 
                  x = "comp", 
                  y = "value",
                  add = "jitter", 
                  short.panel.labs = FALSE,
                  color = "comp",
                  palette = c("purple4" , "darkgreen", "blue3","aquamarine3", "black"),
                  ylab = "Bray-Curtis") + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "kruskal.test") + theme_classic()          
  ?stat_compare_means
  
  
  q1 
  
  # Calcul des comparaisons post-hoc avec t.test
  comp_results <- compare_means(value ~ comp, data = df, method = "wilcox.test", p.adjust.method = "bonferroni")
  
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
    ylim(0, 0.9) +
    stat_summary(fun = mean, 
                 geom = "point", 
                 shape = 23,   # Diamond shape
                 size = 2, 
                 fill = "gray")
  
  
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
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "kruskal.test") + theme_classic()          
  
  
  r1 
  
  # Calcul des comparaisons post-hoc avec t.test
  comp_results <- compare_means(value ~ comp, data = df, method = "wilcox.test", p.adjust.method = "bonferroni")
  
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
    ylim(0, 0.9) +
    stat_summary(fun = mean, 
                 geom = "point", 
                 shape = 23,   # Diamond shape
                 size = 2, 
                 fill = "gray")
  
  
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
    ylim(0, 0.9) +
    stat_summary(fun = mean, 
                 geom = "point", 
                 shape = 23,   # Diamond shape
                 size = 2, 
                 fill = "gray")
  
  
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
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "kruskal.test") + theme_classic()          
  
  
  s1 
  
  # Calcul des comparaisons post-hoc avec t.test
  comp_results <- compare_means(value ~ comp, data = df, method = "wilcox.test", p.adjust.method = "bonferroni")
  
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
    ylim(0, 0.9)  +
    stat_summary(fun = mean, 
                 geom = "point", 
                 shape = 23,   # Diamond shape
                 size = 2, 
                 fill = "gray")
  
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
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + 
    stat_compare_means(label.y = 1, method = "kruskal.test") + theme_classic()          
  
  
  t1 
  
  # Calcul des comparaisons post-hoc avec t.test
  comp_results <- compare_means(value ~ comp, data = df, method = "wilcox.test", p.adjust.method = "bonferroni")
  
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
    ylim(0, 0.9) +
    stat_summary(fun = mean, 
                 geom = "point", 
                 shape = 23,   # Diamond shape
                 size = 2, 
                 fill = "gray")
  
  
  
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
    ylim(0, 0.9) +
    stat_summary(fun = mean, 
                 geom = "point", 
                 shape = 23,   # Diamond shape
                 size = 2, 
                 fill = "gray")
  
  
  # Ajouter les lettres pour montrer les différences significatives
  t1bis_2 <- t1bis_2 + geom_text(data = comp_means, aes(x = comp, y = 0.30, label = letters),
                                 color = "black", vjust = 0)
  
  
  fin_1 <- cowplot::plot_grid(q1bis,  r1bis_2, s1bis, t1bis_2, 
                              ncol = 2,
                              nrow = 2)
  
  fin_2 <- cowplot::plot_grid(q1bis,s1bis, r1bis, t1bis,
                              ncol = 2,
                              nrow = 2)
  
  path_to_boxplot_betadiv_XYZ_bis_2 <- paste0("outputs/boxplot_beta_decomp_XYZ_bis_2_fullsites.pdf")
  ggsave(filename =  path_to_boxplot_betadiv_XYZ_bis_2, plot = fin_1, width = 10, height = 9)
  
  path_to_boxplot_betadiv_XYZ_bis <- paste0("outputs/boxplot_beta_decomp_XYZ_bis_fullsites.pdf")
  ggsave(filename =  path_to_boxplot_betadiv_XYZ_bis, plot = fin_2, width = 10, height = 9)
  
  return(path_to_boxplot_betadiv_XYZ_bis_2)
}
