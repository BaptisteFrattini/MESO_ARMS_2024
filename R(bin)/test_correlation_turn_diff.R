data_and_meta_clean_fullsites = targets::tar_read("clean_data_metadata_fullsites") 
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

S <- vegan::specnumber(data_mean)

mat.diff <- outer(S, S, FUN = "-") 
## Diff ##
df.diff <- melt(as.matrix(mat.diff), varnames = c("row", "col"))
df.diff <- subset(df.diff, row != col)

# Convert factors to characters
df.diff$row <- as.character(df.diff$row)
df.diff$col <- as.character(df.diff$col)

# Create a new column with sorted combinations
df.diff$sorted_comparison <- apply(df.diff[, c("row", "col")], 1, function(x) paste(sort(x), collapse="_"))

# Identify and remove redundant comparisons
df_unique <- df.diff %>% distinct(sorted_comparison, .keep_all = TRUE)

# Remove the temporary column
df.diff <- df_unique[, -ncol(df_unique)]

merged_df <- merge(df.diff, meta_mean, by.x = c("row"), by.y = c("arms"))
merged_df2 <- merge(merged_df, meta_mean, by.x = c("col"), by.y = c("arms"))

df.diff.V <- subset(merged_df2,  campain.x == "RODARMS" & campain.y == "RODARMS")
df.diff.V$comparisons <- ifelse(df.diff.V$site.x == df.diff.V$site.y, "Same_site", "Different_site")
df.diff.V <- subset(df.diff.V, df.diff.V$comparisons == "Different_site")

V = data.frame(value = df.diff.V$value, 
               comp = rep("V", nrow(df.diff.V)))

df.diff.W <- subset(merged_df2, campain.x != "P50ARMS" & campain.y != "P50ARMS")
df.diff.W$comparisons <- ifelse(df.diff.W$campain.x == df.diff.W$campain.y, "Same_Island", "Different_Island")
df.diff.W <- subset(df.diff.W, df.diff.W$comparisons == "Different_Island")


W = data.frame(value = df.diff.W$value, 
               comp = rep("W", nrow(df.diff.W)))

df.diff.Z <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
df.diff.Z <- subset(df.diff.Z, campain.x != campain.y)
df.diff.Z$comparisons <- ifelse(df.diff.Z$site.x == df.diff.Z$site.y, "Same_site", "Different_site")
df.diff.Z <- subset(df.diff.Z, df.diff.Z$comparisons == "Same_site")

Z = data.frame(value = df.diff.Z$value, 
               comp = rep("Z", nrow(df.diff.Z)))

df.diff.Y <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
df.diff.Y <- subset(df.diff.Y, campain.x == campain.y)
df.diff.Y$comparisons <- ifelse(df.diff.Y$site.x == df.diff.Y$site.y, "Same_site", "Different_site")
df.diff.Y <- subset(df.diff.Y, df.diff.Y$comparisons == "Different_site" & df.diff.Y$campain.x == "RUNARMS")

Y = data.frame(value = df.diff.Y$value, 
               comp = rep("Y", nrow(df.diff.Y)))

df.diff.X <- subset(merged_df2, campain.x != "RODARMS" & campain.y != "RODARMS")
df.diff.X <- subset(df.diff.X, campain.x == campain.y)
df.diff.X$comparisons <- ifelse(df.diff.X$site.x == df.diff.X$site.y, "Same_site", "Different_site")
df.diff.X <- subset(df.diff.X, df.diff.X$comparisons == "Different_site" & df.diff.X$campain.x == "P50ARMS")

X = data.frame(value = df.diff.X$value, 
               comp = rep("X", nrow(df.diff.X)))

df_diff <- rbind(V,W,X,Y,Z)

df_diff$abs_value <- abs(df_diff$value)

my_comparisons <- list( c("V", "W"), c("V", "X"), c("V", "Y"), c("V","Z"), c("W","X"), c("W", "Y"), c("W","Z"), c("X","Y"), c("X","Z"), c("Y","Z"))

## Turnover ##

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

df_turn <- rbind(V,W,X,Y,Z)
length(df_diff$value)
length(df_turn$value)

plot(df_diff$abs_value, df_turn$value)

my_comparisons <- list( c("V", "W"), c("V", "X"), c("V", "Y"), c("V","Z"), c("W","X"), c("W", "Y"), c("W","Z"), c("X","Y"), c("X","Z"), c("Y","Z"))

# Richness difference #


abs_diff <- abs(as.vector(as.matrix(mat.diff)))

turn <- as.vector(as.matrix(mat.turn))

plot(abs_diff, turn)
