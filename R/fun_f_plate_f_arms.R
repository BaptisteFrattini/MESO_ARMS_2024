#' occ freq of species in plates arms and sites
#'
#' @param data_and_meta_clean_fullsites the path to the clean data
#' 
#'
#' @return NULL
#' @export
#'
fun_f_plate_f_arms <- function(data_and_meta_clean_fullsites){
  # data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(mgcv)
  data <- read.csv(data_and_meta_clean_fullsites["path_data"], row.names = 1)
  meta <- read.csv(data_and_meta_clean_fullsites["path_meta"], row.names = 1)
  # data <- subset(data, meta$island == "Reunion")
  # meta <- subset(meta, meta$island == "Reunion")
  # 
  
  data <- data[, colSums(data) != 0]
  
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
  data_pa <- vegan::decostand(data, "pa")
  
  species_names <- data.frame(Species = colnames(data))
  
  corr_taxa <- species_names %>%
    mutate(
      taxa = case_when(
        grepl("bry", Species, ignore.case = TRUE) ~ "Bryozoa",
        grepl("SPON", Species, ignore.case = TRUE) ~ "Porifera",
        grepl("ASCC|ASCS", Species, ignore.case = TRUE) ~ "Ascidiacea",
        grepl("FOR", Species, ignore.case = TRUE) ~ "Foraminifera",
        grepl("BIV|pina", Species, ignore.case = TRUE) ~ "Bivalvia",
        grepl("algae|CCA", Species, ignore.case = TRUE) ~ "Algae",
        grepl("BIOF|CYAN", Species, ignore.case = TRUE) ~ "Prokaryotic biotas",
        grepl("WORM", Species, ignore.case = TRUE) ~ "Annelida",
        grepl("HYD|Cni", Species, ignore.case = TRUE) ~ "Cnidaria",
        grepl("Cirr", Species, ignore.case = TRUE) ~ "Cirripedia",
        TRUE ~ "Other"  # Pour toutes les autres espèces
      )
    )
  
  
  
  # fréquence d'occurence à l'echele des faces de plaques ####
  f <- colSums(data_pa)/nrow(data_pa)
  hist(f, breaks = 50)
  
  f_log <- colSums(data_pa)/nrow(data_pa)
  
  f_log <- f_log[f_log != 0]
  
  df <- data.frame(Frequency = f_log)
  
  # Fréquence d'occurence à l'echele des ARMS ####
  
  
  data_mean <- read.csv(data_and_meta_clean_fullsites["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean_fullsites["path_meta_mean"], row.names = 1)
  
  data_mean <- data_mean[, colSums(data_mean) != 0]
  
  data_mean <- data_mean[, msp_list_filter]
  data_mean_pa <- vegan::decostand(data_mean, "pa")
  
  f_ARMS <- colSums(data_mean_pa)/nrow(data_mean_pa)
  hist(f_ARMS, breaks = 50)
  
  
  df <- data.frame(Frequency = f_ARMS)
  
  
  # Fréquence d'occurence à l'echele des site ####
  
  
  data_mean_site <- data %>% 
    group_by(meta$triplicat) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_mean_site <- data_mean_site[,-1]
  
  
  data_mean_site <- data_mean_site[, colSums(data_mean_site) != 0]
  
  data_mean_site <- data_mean_site[, msp_list_filter]
  data_mean_site_pa <- vegan::decostand(data_mean_site, "pa")
  
  f_site <- colSums(data_mean_site_pa)/nrow(data_mean_site_pa)
  
  df <- data.frame(Frequency = f_ARMS)
  
  # Faire un graphique frequence moyenne d'occurence des espèce en fonction de l'emprise #####
  
  tab_f <- data.frame(f_plate = f,
                      f_ARMS = f_ARMS,
                      f_site = f_site)
  
  
  tab_f$species <- rownames(tab_f)
  # Convert to long format
  data_long <- tab_f %>%
    pivot_longer(cols = c(f_plate, f_ARMS, f_site), names_to = "Modality", values_to = "Value")
  
  # # Plot with ggplot2
  
  # Step 1: Calculate the mean values for each modality
  mean_values <- data_long %>%
    group_by(Modality) %>%
    summarise(Mean_Value = mean(Value, na.rm = TRUE))
  
  # Step 2: Perform ANOVA to compare mean values across modalities (optional)
  anova_result <- aov(Value ~ Modality, data = data_long)
  summary(anova_result)
  
  # Step 3: Pairwise comparisons using Student's t-test with Bonferroni correction
  pairwise_result <- pairwise.t.test(data_long$Value, data_long$Modality, p.adjust.method = "bonferroni")
  print(pairwise_result)
  
  # Step 4: Set the factor levels for the Modality column to ensure correct order in the plot
  data_long$Modality <- factor(data_long$Modality, levels = c("f_plate", "f_ARMS", "f_site"))
  
  # Step 5: Create a violin plot to compare the distribution of values across modalities
  ggplot(data_long, aes(x = Modality, y = Value, fill = Modality)) +
    geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot
    geom_point(data = mean_values, aes(x = Modality, y = Mean_Value), 
               shape = 18, color = "purple", size = 4) +  # Mean as grey diamonds
    labs(title = "Distribution of Species Occurrence Frequencies by Modality", 
         x = "Modality", 
         y = "Occurrence Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")  # Choose a color palette
  
  plot(tab_f$f_site, tab_f$f_plate)
  plot(tab_f$f_ARMS, tab_f$f_plate)
  
  # 3. GAM avec pénalisation et k modéré
  mod_gam <- gam(f_plate ~ s(f_ARMS, k = 6), data = tab_f, select = TRUE)
  
  # plot(residuals(mod_gam))
  # k.check(mod_gam)
  par(mfrow = c(2, 2))
  gam.check(mod_gam)
  
  ggplot(tab_f, aes(x = log(f_ARMS), y = log(f_plate))) +
    geom_point() +
    geom_smooth(method = "lm")
  
  #AIC
  
  # 1. Régression polynomiale (ordre 2)
  mod_poly2 <- lm(f_plate ~ poly(f_ARMS, 2), data = tab_f)
  
  # 2. Régression exponentielle (f_plate = a * exp(b * f_ARMS))
  mod_exp <- nls(f_plate ~ a * exp(b * f_ARMS),
                 data = tab_f,
                 start = list(a = 0.01, b = 2),
                 control = list(maxiter = 100))
  
  model_exp <- nls(f_plate ~ a * exp(b * f_ARMS), data = tab_f, start = list(a = 0.01, b = 2))
  
  
  tab_f$pred_exp <- predict(model_exp)
  
  coefs <- coef(model_exp)
  a <- signif(coefs["a"], 3)
  b <- signif(coefs["b"], 3)
  formule_texte <- paste0("f = ", a, " * exp(", b, " * x)")
  
  vv <- ggplot(tab_f, aes(x = f_ARMS, y = f_plate)) +
    geom_point() +
    geom_line(aes(y = pred_exp), color = "red") +
    annotate("text", x = min(tab_f$f_ARMS), y = max(tab_f$f_plate),
             label = formule_texte, hjust = 0, vjust = 1, size = 6, fontface = "italic") +
    labs(
      title = "",
      x = "Occurrence frequency of MSPs on ARMS", 
      y = "Occurrence frequency of MSPs on ARMS plate faces"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16)
    )
  
  # 4. Comparaison des AIC
  aic_values <- AIC(mod_poly2, mod_exp, mod_gam)
  print(aic_values)
  
  path_to_stack_c_chart <-  here::here("outputs/plot_f_plate_f_arms.pdf")
  ggsave(path_to_stack_c_chart, plot = vv, width = 12, height = 9)

  return(NULL)
  
}



