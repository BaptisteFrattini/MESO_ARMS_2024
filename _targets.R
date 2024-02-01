library(targets)
library(tarchetypes)
targets::tar_source()

list(
  tar_target(raw_data, "data/raw-data/data_MESO_ARMS(2).csv", format = "file")
  
  ,tar_target(clean_data_metadata, fun_data_cleaning(raw_data = raw_data))
  
  ,tar_target(explore_data, fun_data_exploring(data_and_meta_clean = clean_data_metadata))

  ,tar_target(stack_c_chart, fun_stack_c_chart(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(venn_diag, fun_venn_diag(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(beta_decomp, fun_decomp_betadiv(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(PCoA, fun_PCoA(data_and_meta_clean = clean_data_metadata))
  
  ,tar_target(boxplot_categories, fun_boxplot_categories(data_and_meta_clean = clean_data_metadata))
  )