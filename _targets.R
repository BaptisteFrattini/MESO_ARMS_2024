library(targets)
library(tarchetypes)
targets::tar_source()

list(
  tar_target(raw_data, "data/raw-data/data_MESO_ARMS.csv", format = "file")
  
)