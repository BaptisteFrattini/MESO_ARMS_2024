library(targets)
library(tarchetypes)
targets::tar_source()

list(

  tar_target(raw_data, "data/raw-data/data_MESO_ARMS(3).csv", format = "file")
  
  #site location ####
  
  ,tar_target(data_gps_sites, "data/raw-data/GPS/Site_location_long_lat.csv", format = "file")
  
  # Temperature ####
  
  ,tar_target(temp_data_1A, "data/raw-data/Temperature/P50ARMS1A.csv", format = "file")
  
  ,tar_target(temp_data_1B, "data/raw-data/Temperature/P50ARMS1B.csv", format = "file")
  
  ,tar_target(temp_data_2A, "data/raw-data/Temperature/P50ARMS2A.csv", format = "file")
  
  ,tar_target(temp_data_2B, "data/raw-data/Temperature/P50ARMS2B.csv", format = "file")
  
  ,tar_target(temp_data_3A, "data/raw-data/Temperature/P50ARMS3A.csv", format = "file")
  
  ,tar_target(temp_data_3B, "data/raw-data/Temperature/P50ARMS3B.csv", format = "file")
  
  ,tar_target(temp_data_shallow, "data/raw-data/Temperature/temperature_in_situ.csv")
  
  ,tar_target(temp_in_situ, fun_temp(temp_1A = temp_data_1A,
                                     temp_1B = temp_data_1B,
                                     temp_2A = temp_data_2A,
                                     temp_2B = temp_data_2B,
                                     temp_3A = temp_data_3A,
                                     temp_3B = temp_data_3B,
                                     temp_shallow = temp_data_shallow))
  
  ,tar_target(temp_copernicus, "data/raw-data/Temperature/cmems_mod_glo_phy_my_0.083deg_P1D-m_1735902316351.nc", format = "file")
  
  ,tar_target(temp_copernicus_2, "data/raw-data/Temperature/cmems_mod_glo_phy_myint_0.083deg_P1D-m_1737626858552.nc", format = "file")
  
  ,tar_target(temp_copernicus_rodrigues, "data/raw-data/Temperature/cmems_mod_glo_phy_my_0.083deg_P1D-m_1738770214961.nc", format = "file")
  
  ,tar_target(temp_copernicus_rodrigues_2, "data/raw-data/Temperature/cmems_mod_glo_phy_myint_0.083deg_P1D-m_1738676877692.nc", format = "file")
  
  ,tar_target(temp_extraction, fun_temperature_extraction(data_temp_copernicus = temp_copernicus,
                                                          data_temp_copernicus_2 = temp_copernicus_2,
                                                          data_temp_copernicus_rodrigues = temp_copernicus_rodrigues,
                                                          data_temp_copernicus_rodrigues_2 = temp_copernicus_rodrigues_2,
                                                          gps_sites = data_gps_sites))
  
  ,tar_target(temp_comparison, fun_temperature_comparison(temp_extraction = temp_extraction,
                                                          temp_in_situ = temp_in_situ))
  
  # Cleaning data ####
  
  ,tar_target(clean_data_metadata, fun_data_cleaning(raw_data = raw_data))
  
  ,tar_target(clean_data_metadata_fullsites, fun_data_cleaning_fullsites(raw_data = raw_data))
  
  # Data Analysis ####
  
  ,tar_target(explore_data, fun_data_exploring(data_and_meta_clean = clean_data_metadata))

  ,tar_target(stack_c_chart, fun_stack_c_chart(data_and_meta_clean = clean_data_metadata))

  ,tar_target(venn_diag, fun_venn_diag(data_and_meta_clean = clean_data_metadata))
   
  ,tar_target(beta_decomp, fun_decomp_betadiv(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(permanova_simper, fun_permanova(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(mantel_test, fun_mantel(data_and_meta_clean = clean_data_metadata,
                                              gps_sites = data_gps_sites))
  
  ,tar_target(graph, fun_graph(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(taxo_overlap, fun_taxo_overlap(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(taxo_overlap_perc, fun_taxo_overlap_perc(data_and_meta_clean = clean_data_metadata))
 
  ,tar_target(taxo_overlap_rarity, fun_taxo_overlap_rarity(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(taxo_overlap_quartile, fun_taxo_overlap_quart(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(LCBD, fun_LCBD(data_and_meta_clean = clean_data_metadata))
  
  
  # 
  # ,tar_target(PCoA, fun_PCoA(data_and_meta_clean = clean_data_metadata))
  # 
  # ,tar_target(boxplot_categories, fun_boxplot_categories(data_and_meta_clean = clean_data_metadata))
  # 
  # ,tarchetypes::tar_quarto(report, "quarto_file.qmd")
  
  # Data analysis - fullsites ####
  
  ,tar_target(permanova_simper_fullsites, fun_permanova_fullsites(data_and_meta_clean_fullsites = clean_data_metadata_fullsites))
  
  ,tar_target(taxo_overlap_rarity_fullsites, fun_taxo_overlap_rarity_fullsites(data_and_meta_clean_fullsites = clean_data_metadata_fullsites))
  
  ,tar_target(data_exploring_fullsites, fun_data_exploring_fullsites(data_and_meta_clean_fullsites = clean_data_metadata_fullsites))
  
  ,tar_target(LCBD_fullsites, fun_LCBD_fullsites(data_and_meta_clean_fullsites = clean_data_metadata_fullsites))
  
  )