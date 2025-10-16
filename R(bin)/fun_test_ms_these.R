# raw_data = targets::tar_read("raw_data") 

library(dplyr)
library(vegan)
library(ggplot2)

dat_path <- here::here(raw_data)
data_raw <- read.table(dat_path, 
                       header = TRUE, 
                       sep = ";", 
                       dec = ",")

#### change the morpho-species names ####

meta_data <- data_raw[-(nrow(data_raw)),c(1:4)]
data <- data_raw[-(nrow(data_raw)), c(5:ncol(data_raw))]
data <- data[,!colSums(data) == 0]

data_filtered <- data[grepl("RUNARMS1|RUNARMS2|RUNARMS3|RUNARMS4|RUNARMS5|RUNARMS6|RUNARMS7|RUNARMS8|RUNARMS9|P50|RODRARMS|CINARMS", meta_data$Image.name),]
meta_filtered <- meta_data[grepl("RUNARMS1|RUNARMS2|RUNARMS3|RUNARMS4|RUNARMS5|RUNARMS6|RUNARMS7|RUNARMS8|RUNARMS9|P50|RODRARMS|CINARMS", meta_data$Image.name),]

data_filtered2 <- data_filtered[!grepl("Unconfirmed", meta_filtered$Annotation.status),]
meta_filtered2 <- meta_filtered[!grepl("Unconfirmed", meta_filtered$Annotation.status),]

data_filtered2 <- data_filtered2[!grepl("9T", meta_filtered2$Image.name),]
meta_filtered2 <- meta_filtered2[!grepl("9T", meta_filtered2$Image.name),]


rownames(data_filtered2) <- substr(meta_filtered2$Image.name, 1, nchar(meta_filtered2$Image.name) - 4)

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

data_filtered2$MSP15_SPON <- data_filtered2$MSP13_SP + data_filtered2$MSP15_SP + data_filtered2$MSP17_SP
data_filtered2$MSP13_SPON <- NULL
data_filtered2$MSP17_SPON <- NULL

data_filtered2$MSP9_SPON <- data_filtered2$MSP32_SP + data_filtered2$MSP9_SP
data_filtered2$MSP32_SPON <- NULL

data_filtered2$MSP15_ASCC <- data_filtered2$MSP6_SP + data_filtered2$MSP15_ASCC
data_filtered2$MSP6_SPON <- NULL

data_filtered2$MSP10_SPON <- data_filtered2$MSP40_SP + data_filtered2$MSP10_SP
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
  rename(Bry_Schizoporellidae_sp1 = MSP12_BRYO)
data_filtered2 <- data_filtered2 %>%
  rename(Bry_Schizoporellidae_sp2 = MSP14_BRYO)
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
data_filtered2 <- data_filtered2 %>%
  rename(Ascc_Polysyncraton_sp1 = MSP37_ASCC)
data_filtered2 <- data_filtered2 %>%
  rename(Ascc_Symplegma_sp1 = MSP38_ASCC)
data_filtered2 <- data_filtered2 %>%
  rename(Ascc_Polyclinum_sp1 = MSP39_ASCC)
data_filtered2 <- data_filtered2 %>%
  rename(Ascc_Botrylloides_sp1 = MSP40_ASCC)
data_filtered2 <- data_filtered2 %>%
  rename(Ascc_Didemnidae_sp4 = MSP42_ASCC)
data_filtered2 <- data_filtered2 %>%
  rename(Ascc_Symplegma_sp2 = MSP45_ASCC)
data_filtered2 <- data_filtered2 %>%
  rename(Ascc_Botryllus_sp7 = MSP46_ASCC)

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
rownames(data) <- gsub("RODRARMS", "RODARMS", rownames(data))

colnames(data)
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


nrow(corr_taxa[corr_taxa$taxa == "Ascidiacea",])
