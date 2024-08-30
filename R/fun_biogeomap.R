#' testing biogeonat
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

biogeo_test <- function(metadata_data_mean){
 
  # data_and_meta_clean = targets::tar_read("clean_data_metadata") 
  library(plyr)
  library(igraph)
  df_mean <- read.csv(data_and_meta_clean["path_data_mean"], row.names = 1)
  meta_mean <- read.csv(data_and_meta_clean["path_meta_mean"], row.names = 1)
  
  df_mean <- df_mean %>% 
    group_by(meta_mean$campain) %>% 
    summarise_all(sum, na.rm = TRUE)
  
  df_mean <- as.data.frame(df_mean)
  
  #### step 1 ####
  rownames(df_mean) = NULL
  data <- reshape2::melt(df_mean)
  nrow(data)

  data = data.frame(arms = as.factor(data$'meta_mean$campain'),
                    msp = as.factor(data$variable),
                    abun = data$value)
  
  class(data$abun)
  
  data <- data[which(data$abun > 0, arr.ind = TRUE),]
  
  g <- graph_from_data_frame(
    d = data, 
    directed = FALSE
  )
 
  V(g)$type <- V(g)$name %in% rownames(data)
  
  plot(g, vertex.label = V(g)$name, vertex.size = 5, 
       vertex.color = ifelse(V(g)$type, "lightblue", "salmon"),
       edge.width = E(g)$weight)
  
  
  write_graph(g, file = "data/derived-data/Graph.gephi.graphml", format = "graphml")
  
  return()
  }  
    