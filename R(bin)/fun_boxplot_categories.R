#' Boxplot of categories
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to...
#' @export
#'
fun_boxplot_categories <- function(data_and_meta_clean){
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  library(ggpubr)
  library(forcats)
  library(ggplot2)
  library(dplyr)
  
  data <- read.csv(data_and_meta_clean["path_data_pool"], row.names = 1)
  meta <- read.csv(data_and_meta_clean["path_meta"], row.names = 1)
  
  df <- data.frame(data, meta)
  
  
  
  # Exemple de données
  # data <- data.frame(
  #   Site = rep(c("Site1", "Site2", "Site3"), each = 6),
  #   Depth = rep(c("Depth1", "Depth2"), times = 9),
  #   Value = rnorm(18)  # Remplacez ceci par vos valeurs réelles
  # )
  
  # Définir l'ordre souhaité des sites
  site_order <- c("Cap_Houssaye", "Saint_Leu","Grand_Bois")
  
  ggplot(df, aes(x=bivalvia)) + 
    geom_density()
  ggpubr::ggqqplot(df$bivalvia)
  shapiro.test(df$bivalvia)
  
  library(ggsignif)
  colors = c("blue4", "aquamarine3")
  
  # bivalvia ####
  ggplot(df, aes(x=bivalvia)) + 
    geom_density()
  ggpubr::ggqqplot(df$bivalvia)
  shapiro.test(df$bivalvia)
  
  
  p1 <- ggboxplot(df, x = "campain", y = "bivalvia",
                 color = "campain", palette = colors,
                 add = "jitter",
                 facet.by = "site", short.panel.labs = FALSE) 
  
  p1 <- p1 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  
  
  # porifera ####
  ggplot(df, aes(x=porifera)) + 
    geom_density()
  ggpubr::ggqqplot(df$porifera)
  shapiro.test(df$porifera)
  
  
  p2 <- ggboxplot(df, x = "campain", y = "porifera",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p2 <- p2 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  
  # bryozoa ####
  ggplot(df, aes(x=bryozoa)) + 
    geom_density()
  ggpubr::ggqqplot(df$bryozoa)
  shapiro.test(df$bryozoa)
  
  
  p3 <- ggboxplot(df, x = "campain", y = "bryozoa",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p3 <- p3 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  
  # ascidiacea ####
  df$ascidiacea <- df$ascidiacea_c+df$ascidiacea_s
  
  ggplot(df, aes(x=ascidiacea)) + 
    geom_density()
  ggpubr::ggqqplot(df$ascidiacea)
  shapiro.test(df$ascidiacea)
  
  
  p3 <- ggboxplot(df, x = "campain", y = "ascidiacea",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p3 <- p3 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  
  # foraminifera ####
  
  ggplot(df, aes(x=foraminifera)) + 
    geom_density()
  ggpubr::ggqqplot(df$foraminifera)
  shapiro.test(df$foraminifera)
  
  
  p4 <- ggboxplot(df, x = "campain", y = "foraminifera",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p4 <- p4 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  
  # annelida ####
  
  ggplot(df, aes(x=annelida)) + 
    geom_density()
  ggpubr::ggqqplot(df$annelida)
  shapiro.test(df$annelida)
  
  
  p5 <- ggboxplot(df, x = "campain", y = "annelida",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p5 <- p5 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  
  # hydrozoa ####
  
  ggplot(df, aes(x=hydrozoa)) + 
    geom_density()
  ggpubr::ggqqplot(df$hydrozoa)
  shapiro.test(df$hydrozoa)
  
  
  p6 <- ggboxplot(df, x = "campain", y = "hydrozoa",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p6 <- p6 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  
  # CCA ####
  
  ggplot(df, aes(x=CCA)) + 
    geom_density()
  ggpubr::ggqqplot(df$CCA)
  shapiro.test(df$CCA)
  
  
  p7 <- ggboxplot(df, x = "campain", y = "CCA",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p7 <- p7 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  # No_Recruitment ####
  
  ggplot(df, aes(x=No_Recruitment)) + 
    geom_density()
  ggpubr::ggqqplot(df$No_Recruitment)
  shapiro.test(df$No_Recruitment)
  
  
  p8 <- ggboxplot(df, x = "campain", y = "No_Recruitment",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p8 <- p8 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  # Sediments ####
  
  ggplot(df, aes(x=Sediments)) + 
    geom_density()
  ggpubr::ggqqplot(df$Sediments)
  shapiro.test(df$Sediments)
  
  
  p9 <- ggboxplot(df, x = "campain", y = "Sediments",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p9 <- p9 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  # Prokariotic biotas ####
  
  ggplot(df, aes(x=prokariotic_biotas)) + 
    geom_density()
  ggpubr::ggqqplot(df$prokariotic_biotas)
  shapiro.test(df$prokariotic_biotas)
  
  
  p10 <- ggboxplot(df, x = "campain", y = "prokariotic_biotas",
                  color = "campain", palette = colors,
                  add = "jitter",
                  facet.by = "site", short.panel.labs = FALSE) 
  
  p10 <- p10 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  # other_algae ####
  
  ggplot(df, aes(x=other_algae)) + 
    geom_density()
  ggpubr::ggqqplot(df$other_algae)
  shapiro.test(df$other_algae)
  
  
  p11 <- ggboxplot(df, x = "campain", y = "other_algae",
                   color = "campain", palette = colors,
                   add = "jitter",
                   facet.by = "site", short.panel.labs = FALSE) 
  
  p11 <- p11 + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")
  
  fin <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, 
                            ncol = 3,
                            nrow = 3)
  path_to_boxplot_categories <- paste0("outputs/boxplot_categories.pdf")
  ggsave(filename =  path_to_boxplot_categories, plot = fin, width = 19, height = 13.5)
  
  # # Créer un boxplot avec des paires de boîtes pour chaque site et chaque profondeur
  # p <- ggplot(df, aes(x = factor(site, levels = site_order), y = bivalvia, fill = campain)) +
  #      geom_boxplot(position = position_dodge(0.8), width = 0.7) +
  #      labs(title = "Boxplot par Site et Profondeur", x = "Site", y = "mean cover") +
  #      theme_minimal() +
  #      scale_fill_manual(values = colors) 
  # 
  # p + geom_text(aes(x = factor(site, levels = site_order), y = max(df$bivalvia) + 0.2, label = ifelse(p_values < 0.05, "*", "")),
  #               position = position_dodge(0.8), size = 4, vjust = 0)
  # 
  # p + geom_signif(comparisons = list(c("P50ARMS", "RUNARMS")), textsize = 5, vjust = -0.5)
  
  
  return(path_to_boxplot_categories)
}
  