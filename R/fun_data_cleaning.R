#' Import data MESO_ARMS
#'
#' @param raw_data the path to the raw data file
#' @param arms_id the ID of the arms to subset for
#'
#' @return the path to the subseted raw data file
#' @export
#'
fun_data_cleaning <- function(raw_data){
  # raw_data = targets::tar_read("raw_data") 
  
  library(dplyr)
  library(vegan)
  
  dat_path <- here::here(raw_data)
  data_raw <- read.table(dat_path, 
                     header = TRUE, 
                     sep = ";", 
                     dec = ",")
  
  #### change the morpho-species names ####
  
  meta_data <- data_raw[-(nrow(data_raw)),c(1:4)]
  data <- data_raw[-(nrow(data_raw)), c(5:ncol(data_raw))]
  data <- data[,!colSums(data) == 0]
  
  data_filtered <- data[grepl("RUNARMS1|RUNARMS5|RUNARMS9|P50", meta_data$Image.name),]
  meta_filtered <- meta_data[grepl("RUNARMS1|RUNARMS5|RUNARMS9|P50", meta_data$Image.name),]
  
  data_filtered2 <- data_filtered[!grepl("9T", meta_filtered$Image.name),]
  meta_filtered2 <- meta_filtered[!grepl("9T", meta_filtered$Image.name),]
  
  rownames(data_filtered2) <- substr(meta_filtered2$Image.name, 1, 12)
  
  data_filtered2$MSP1_BRYO <- NULL
  
  data_filtered2$MSP8_BRYO <- data_filtered2$MSP18_BRYO + data_filtered2$MSP8_BRYO
  data_filtered2$MSP18_BRYO <- NULL
  
  data_filtered2$MSP3_BRYO <- data_filtered2$MSP15_BRYO + data_filtered2$MSP3_BRYO
  data_filtered2$MSP15_BRYO <- NULL
  
  data_filtered2$MSP2_ASCC <- data_filtered2$MSP2_ASCC + data_filtered2$MSP18_ASCC
  data_filtered2$MSP18_ASCC <- NULL
  
  data_filtered2$MSP8_ASCC <- data_filtered2$MSP8_ASCC + data_filtered2$MSP20_ASCC
  data_filtered2$MSP20_ASCC <- NULL
  
  data_filtered2$MSP16_ASCC <- data_filtered2$MSP16_ASCC + data_filtered2$MSP23_ASCC
  data_filtered2$MSP23_ASCC <- NULL
  
  data_filtered2$MSP3_ASCS <- data_filtered2$MSP3_ASCS + data_filtered2$MSP7_ASCS
  data_filtered2$MSP7_ASCS <- NULL
  
  data_filtered2$MSP26_ASCC <- data_filtered2$MSP31_ASCC + data_filtered2$MSP26_ASCC
  data_filtered2$MSP31_ASCC <- NULL
  
  data_filtered2$MSP15_SPON <- data_filtered2$MSP13_SP + data_filtered2$MSP15_SPON + data_filtered2$MSP17_SPON
  data_filtered2$MSP13_SPON <- NULL
  data_filtered2$MSP17_SPON <- NULL
  
  data_filtered2$MSP9_SPON <- data_filtered2$MSP32_SPON + data_filtered2$MSP9_SPON
  data_filtered2$MSP32_SPON <- NULL
  
  data_filtered2$MSP15_ASCC <- data_filtered2$MSP6_SPON + data_filtered2$MSP15_ASCC
  data_filtered2$MSP6_SPON <- NULL
  
  data_filtered2$MSP10_SPON <- data_filtered2$MSP40_SPON + data_filtered2$MSP10_SPON
  data_filtered2$MSP40_SPON <- NULL
  
  data_filtered2$OROS <- data_filtered2$OROS + data_filtered2$X_SED
  data_filtered2$X_SED <- NULL
  
  # rename foram
  data_filtered2 <- data_filtered2 %>%
    rename(For_Miniacina_sp_red = X_BRY)
  data_filtered2 <- data_filtered2 %>%
    rename(For_Miniacina_sp_juv = MSP9_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(For_Miniacina_sp_adu = MSP10_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Other_foraminifera = X_FORM)
  
  # rename Bryozoa
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Parasmittina_margaritata = MSP11_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Smittipora_harmeriana = MSP3_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Schizomavella_sp = MSP16_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Smittoidea_cautela = MSP17_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Watersipora_subtorquata = MSP2_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Disporella_novaehollandiae = MSP4_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Crisia_elongata = MSP6_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Tubulipora_sp = MSP7_BRYO)
  data_filtered2 <- data_filtered2 %>%
    rename(Bry_Scrupocellaria_sp = MSP8_BRYO)
  
  # rename Ascidie
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Aplidium_sp = MSP1_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascs_Polycarpa_sp = MSP1_ASCS)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Lissoclinum_sp1 = MSP10_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Ascidia_archaia = MSP10_ASCS)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Lissoclinum_sp2 = MSP11_ASCS)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Didemnidae_sp1 = MSP11_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Cystodystes_sp = MSP12_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Polysyncraton_milleporae = MSP13_ASCS)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Botryllus_sp1 = MSP14_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Botryllus_gregalis = MSP16_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Botryllus_sp2 = MSP19_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Eusynstyela_hartmeyeri = MSP2_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Trididemnum_sp = MSP21_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Botryllus_sp3 = MSP22_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Botryllus_sp4 = MSP3_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Styela_canopus = MSP3_ASCS)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Didemnidae_sp2 = MSP4_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Ascidia_sydneiensis = MSP4_ASCS)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Botryllus_sp5 = MSP5_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Polysyncraton_rostrum = MSP6_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Didemnidae_sp3 = MSP8_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Ascidia_fictile = MSP8_ASCS)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Botryllus_tuberatus = MSP24_ASCC)
  data_filtered2 <- data_filtered2 %>%
    rename(Ascc_Botryllus_sp6 = MSP25_ASCC)
  

  # Other
  data_filtered2 <- data_filtered2 %>%
    rename(Erect_Phaeophyceae_algae = X_BRUP)
  data_filtered2 <- data_filtered2 %>%
    rename(Erect_Chlorophyta_algae = X_GRUP)
  data_filtered2 <- data_filtered2 %>%
    rename(Erect_Rhodophyta_algae = X_RDUP)
  data_filtered2 <- data_filtered2 %>%
    rename(Encrusting_Phaeophyceae_algae = X_BREN)
  data_filtered2 <- data_filtered2 %>%
    rename(Encrusting_Chlorophyta_algae = X_GREN)
  data_filtered2 <- data_filtered2 %>%
    rename(No_Recruitment = X_NR)
  data_filtered2 <- data_filtered2 %>%
    rename(Sediments = OROS)
  data_filtered2 <- data_filtered2 %>%
    rename(Bivalvia = X_BI)
  data_filtered2 <- data_filtered2 %>%
    rename(Calcareous_worm_tubes = X_CAWT)
  data_filtered2 <- data_filtered2 %>%
    rename(Soft_worm_tubes = X_SOWT)
  data_filtered2 <- data_filtered2 %>%
    rename(Eggs = X_EGG)
  data_filtered2 <- data_filtered2 %>%
    rename(Cni_Plumulariidae = X_HYD)
  data_filtered2 <- data_filtered2 %>%
    rename(CCA = X_CCA)
  data_filtered2 <- data_filtered2 %>%
    rename(Cirripedia = Root.Barn)
  
  #### filtering and transform data ####  
  data_filtered2 <- data_filtered2[, !(colnames(data_filtered2) %in% c("X_UNAV", "X_UNK", "UNK_2", "X_CMOR", "X_OMOL", "X_TUNC", "TUNI","X_CO","X_SP", "X_MOBF", "Eggs"))]
  colnames(data_filtered2)
  data_filtered2 <- vegan::decostand(data_filtered2, method = "total")
  data <- data_filtered2*100
  rowSums(data)
  
  #### build meta_data and data ####
  name <- rownames(data)
  arms <- substr(rownames(data), 1,9)
  campain <- substr(rownames(data), 1,7)
  triplicat <- substr(rownames(data), 1,8) 
  meta_data <- data.frame(name, arms, campain, triplicat)
  
  path <- "data/derived-data/" 
  
  name_data <- "data_clean.csv" 
  name_metadata <- "metadata_clean.csv"

  #dans ce dossier 
  meta_out_path <- paste0(path, name_data) #Nom du fichier de 
  data_out_path <- paste0(path, name_metadata)
  #metadata généré dans ce dossier

  write.csv(data, file = data_out_path, row.names = TRUE)
  write.csv(meta_data, file = meta_out_path, row.names = TRUE)
  
  
  #### Build the mean meta data and data ####
  
  data_mean <- data %>% 
    group_by(meta_data$arms) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_mean <- as.data.frame(data_mean)
  row.names(data_mean) <- data_mean$`meta_data$arms`
  data_mean <- data_mean[,-1]
 
  arms <- rownames(data_mean)
  campain <- substr(rownames(data_mean), 1,7)
  triplicat <- substr(rownames(data_mean), 1,8) 
  meta_data_mean <- data.frame(arms, campain, triplicat)
  
  path <- "data/derived-data/" 
  
  name_data_mean <- "data_mean_clean.csv" 
  name_metadata_mean <- "metadata_mean_clean.csv"
  
  #dans ce dossier 
  meta_mean_out_path <- paste0(path, name_metadata_mean) #Nom du fichier de 
  data_mean_out_path <- paste0(path, name_data_mean)
  
  #metadata généré dans ce dossier
  write.csv(data_mean, file = data_mean_out_path, row.names = TRUE)
  write.csv(meta_data_mean, file = meta_mean_out_path, row.names = TRUE)

  
  
  
  #### Build the pool data ####

  spo_columns <- grep("SPON", names(data), value = TRUE)
  spo_mean <- rowSums(data[, spo_columns])
  
  bryo_columns <- c(grep("Bry", names(data), value = TRUE),grep("BRY", names(data), value = TRUE))
  bryo_mean <- rowSums(data[,bryo_columns])
  
  ascc_columns <- c(grep("Ascc", names(data), value = TRUE),grep("ASCC", names(data), value = TRUE))
  ascc_mean <- rowSums(data[,ascc_columns])
  
  ascs_columns <- c(grep("Ascs", names(data), value = TRUE),grep("ASCS", names(data), value = TRUE))
  ascs_mean <- rowSums(data[,ascs_columns])
  
  hyd_columns <- c(grep("HYD", names(data), value = TRUE), grep("Cni", names(data), value = TRUE))
  hyd_mean <- rowSums(data[, hyd_columns])
  
  biv_columns <- c(grep("Biv", names(data), value = TRUE), "Pina")
  biv_mean <- rowSums(data[, biv_columns])
  
  for_columns <- c(grep("For", names(data), value = TRUE),"Other_foraminifera")
  for_mean <- rowSums(data[, for_columns])
  
  algae_columns <- grep("algae", names(data), value = TRUE)
  algae_mean <- rowSums(data[,algae_columns])
  
  worm_columns <- grep("worm", names(data), value = TRUE)
  worm_mean <- rowSums(data[,worm_columns])
  
  prokariot_columns <- c(grep("BIOF", names(data), value = TRUE), grep("CYAN", names(data), value = TRUE))
  prokariot_mean <- rowSums(data[,prokariot_columns])
  
  rest_columns <- names(data)
  rest_columns <- rest_columns[rest_columns %in% c("Cirripedia", "No_Recruitment", "Sediments", "CCA")]
  rest_data <- data[,rest_columns]
  
  data_pool <- data.frame(porifera = spo_mean,
                          bryozoa = bryo_mean,
                          ascidiacea_c = ascc_mean,
                          ascidiacea_s = ascs_mean,
                          foraminifera = for_mean,
                          other_algae = algae_mean,
                          annelida = worm_mean,
                          bivalvia = biv_mean,
                          hydrozoa = hyd_mean,
                          prokariotic_biotas = prokariot_mean)
  
  data_pool <- data.frame(cbind(data_pool,rest_data))
  
  path <- "data/derived-data/" 
  
  name_data_mean <- "data_pool_clean.csv" 

  data_pool_out_path <- paste0(path, name_data_mean)
  
  write.csv(data_pool, file = data_pool_out_path, row.names = TRUE)

  #### Build the pool & mean data ####
  
  data_pool_mean <- data_pool %>% 
    group_by(meta_data$arms) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  data_pool_mean <- as.data.frame(data_pool_mean)
  row.names(data_pool_mean) <- data_pool_mean$`meta_data$arms`
  data_pool_mean <- data_pool_mean[,-1]
  rowSums(data_pool_mean)
  
  path <- "data/derived-data/" 
  
  name_data_mean_pool <- "data_mean_pool_clean.csv" 
  
  data_mean_pool_out_path <- paste0(path, name_data_mean_pool)
  
  write.csv(data_pool_mean, file = data_mean_pool_out_path, row.names = TRUE)
  
  res <- c(data_out_path, meta_out_path, data_mean_out_path, meta_mean_out_path, data_pool_out_path, data_mean_pool_out_path)
  names(res) <- c("path_data", "path_meta", "path_data_mean", "path_meta_mean", "path_data_pool", "path_data_pool_mean")
  
  return(res)
  
}

