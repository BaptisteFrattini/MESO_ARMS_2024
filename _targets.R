library(targets)
library(tarchetypes)
targets::tar_source()

list(
  tar_target(raw_data, "data/raw-data/data_MESO_ARMS(2).csv", format = "file")
  
  # Temperature ####
  
  ,tar_target(temp_data_1A, "data/raw-data/Temperature/P50ARMS1A.csv", format = "file")
  
  ,tar_target(temp_data_1B, "data/raw-data/Temperature/P50ARMS1B.csv", format = "file")
  
  ,tar_target(temp_data_2A, "data/raw-data/Temperature/P50ARMS2A.csv", format = "file")
  
  ,tar_target(temp_data_2B, "data/raw-data/Temperature/P50ARMS2B.csv", format = "file")
  
  ,tar_target(temp_data_3A, "data/raw-data/Temperature/P50ARMS3A.csv", format = "file")
  
  ,tar_target(temp_data_3B, "data/raw-data/Temperature/P50ARMS3B.csv", format = "file")
  
  ,tar_target(temp_in_situ, fun_temp(temp_1A = temp_data_1A,
                                     temp_1B = temp_data_1B,
                                     temp_2A = temp_data_2A,
                                     temp_2B = temp_data_2B,
                                     temp_3A = temp_data_3A,
                                     temp_3B = temp_data_3B))
  
  # Cleaning data ####
  
  ,tar_target(clean_data_metadata, fun_data_cleaning(raw_data = raw_data))
  
  # Data Analysis ####
  
  ,tar_target(explore_data, fun_data_exploring(data_and_meta_clean = clean_data_metadata))

  ,tar_target(stack_c_chart, fun_stack_c_chart(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(venn_diag, fun_venn_diag(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(beta_decomp, fun_decomp_betadiv(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(PCoA, fun_PCoA(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(boxplot_categories, fun_boxplot_categories(data_and_meta_clean = clean_data_metadata))
  
  ,tarchetypes::tar_quarto(report, "quarto_file.qmd")
  )