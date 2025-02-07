#' Temperature model copernicus
#'
#' @param data_and_meta_clean the path to the clean data
#' 
#'
#' @return the path to the subseted raw data file
#' @export
#'
fun_temperature_extraction <- function(data_temp_copernicus, 
                                       data_temp_copernicus_2, 
                                       data_temp_copernicus_rodrigues, 
                                       data_temp_copernicus_rodrigues_2,
                                       gps_sites){
   # 
   # data_temp_copernicus = targets::tar_read("temp_copernicus")
   # data_temp_copernicus_2 = targets::tar_read("temp_copernicus_2")
   # data_temp_copernicus_rodrigues = targets::tar_read("temp_copernicus_rodrigues")
   # data_temp_copernicus_rodrigues_2 = targets::tar_read("temp_copernicus_rodrigues_2")
   # gps_sites = targets::tar_read("data_gps_sites")

  data_gps <- read.csv(gps_sites, dec = ",",sep = ";", row.names = 1)
  data_gps_RUNA <- data_gps[c(7:15),] 
  rows <- rownames(data_gps_RUNA)
  data_gps_RUNA <- data.frame(Longitude = data_gps_RUNA$Longitude,
                              Latitude = data_gps_RUNA$Latitude)
  rownames(data_gps_RUNA) <- rows
  
  
  library(ncdf4)
  #### import raster ####
  dat.ncfd <- ncdf4::nc_open(data_temp_copernicus)
  
  attributes(dat.ncfd$var)
  attributes(dat.ncfd$dim)
  lat <- ncdf4::ncvar_get(dat.ncfd, "latitude")
  long <- ncdf4::ncvar_get(dat.ncfd, "longitude")
  depth <- ncdf4::ncvar_get(dat.ncfd, "depth")
  time <- ncdf4::ncvar_get(dat.ncfd, "time")
  head(time)
  tunits <- ncatt_get(dat.ncfd, "time", "units") #check units
  nt <- dim(time)
  lswt_array <- ncvar_get(dat.ncfd, "thetao")
  dim(lswt_array)
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="GMT")
  range(time_obs)
  length(time_obs)
  
  
  lonlattime <- as.matrix(expand.grid(long,lat,time_obs))
  lswt_vec_long <- as.vector(lswt_array)
  lswt_obs <- data.frame(cbind(lonlattime, lswt_vec_long))
  colnames(lswt_obs) <- c("x","y","z","l")
  lswt_obs$l <- as.numeric(lswt_obs$l)
  lswt_obs_clean <- na.omit(lswt_obs)
  
  
  n <- length(levels(as.factor(lswt_obs_clean$z)))

  
  ####curve RUNA1####
  
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~x+y, sub, mean)   
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RUNA[1,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]  
    
  }
  
  T.RUNA1 <- as.data.frame(cbind(temp.moy, timeline))
  
  ####curve RUNA5####
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~x+y, sub, mean)   
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RUNA[5,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]  
    
  }
  
  T.RUNA5 <- as.data.frame(cbind(temp.moy, timeline))
  
  ####curve RUNA9####
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  
  data_gps_RUNA[9,2] <- -21.38
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~x+y, sub, mean)   
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RUNA[9,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]  
    
  }
  
  T.RUNA9 <- as.data.frame(cbind(temp.moy, timeline))
  
  data_copernicus_temperature_shallow_2019_20 <- as.data.frame(cbind(T.RUNA1$timeline, T.RUNA1$temp.moy, T.RUNA5$temp.moy, T.RUNA9$temp.moy))
  colnames(data_copernicus_temperature_shallow_2019_20) <- c("date", "RUNA1", "RUNA5", "RUNA9")
  
  
  data_copernicus_temperature_shallow_2019_20$RUNA1 <- zoo::rollmean(as.numeric(data_copernicus_temperature_shallow_2019_20$RUNA1), k = 7, align = "left", fill = NA)
  data_copernicus_temperature_shallow_2019_20$RUNA5 <- zoo::rollmean(as.numeric(data_copernicus_temperature_shallow_2019_20$RUNA5), k = 7, align = "left", fill = NA)
  data_copernicus_temperature_shallow_2019_20$RUNA9 <- zoo::rollmean(as.numeric(data_copernicus_temperature_shallow_2019_20$RUNA9), k = 7, align = "left", fill = NA)
  
  path_to_data_copernicus_shallow_2019_20 <- "data/derived-data/data_copernicus_temperature_shallow_2019_20.csv"
  
  write.csv(data_copernicus_temperature_shallow_2019_20, path_to_data_copernicus_shallow_2019_20)

  
  
  #### import raster 2 ####
  dat.ncfd <- ncdf4::nc_open(data_temp_copernicus_2)
  
  attributes(dat.ncfd$var)
  attributes(dat.ncfd$dim)
  lat <- ncdf4::ncvar_get(dat.ncfd, "latitude")
  long <- ncdf4::ncvar_get(dat.ncfd, "longitude")
  depth <- ncdf4::ncvar_get(dat.ncfd, "depth")
  time <- ncdf4::ncvar_get(dat.ncfd, "time")
  head(time)
  tunits <- ncatt_get(dat.ncfd, "time", "units") #check units
  nt <- dim(time)
  lswt_array <- ncvar_get(dat.ncfd, "thetao")
  dim(lswt_array)
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="GMT")
  range(time_obs)
  length(time_obs)
  
  
  lonlattime <- as.matrix(expand.grid(long,lat,time_obs))
  lswt_vec_long <- as.vector(lswt_array)
  lswt_obs <- data.frame(cbind(lonlattime, lswt_vec_long))
  colnames(lswt_obs) <- c("x","y","z","l")
  lswt_obs$l <- as.numeric(lswt_obs$l)
  lswt_obs_clean <- na.omit(lswt_obs)
  
  
  n <- length(levels(as.factor(lswt_obs_clean$z)))
  
  ####curve RUNA1####
  
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~x+y, sub, mean)   
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RUNA[1,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]  
    
  }
  
  T.RUNA1 <- as.data.frame(cbind(temp.moy, timeline))
  
  ####curve RUNA5####
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~x+y, sub, mean)   
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RUNA[5,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]  
    
  }
  
  T.RUNA5 <- as.data.frame(cbind(temp.moy, timeline))
  
  ####curve RUNA9####
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  
  data_gps_RUNA[9,2] <- -21.38
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~x+y, sub, mean)   
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RUNA[9,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]  
    
  }
  
  T.RUNA9 <- as.data.frame(cbind(temp.moy, timeline))
  
  data_copernicus_temperature_shallow_2022_23 <- as.data.frame(cbind(T.RUNA1$timeline, T.RUNA1$temp.moy, T.RUNA5$temp.moy, T.RUNA9$temp.moy))
  colnames(data_copernicus_temperature_shallow_2022_23) <- c("date", "RUNA1", "RUNA5", "RUNA9")
  
  
  data_copernicus_temperature_shallow_2022_23$RUNA1 <- zoo::rollmean(as.numeric(data_copernicus_temperature_shallow_2022_23$RUNA1), k = 7, align = "left", fill = NA)
  data_copernicus_temperature_shallow_2022_23$RUNA5 <- zoo::rollmean(as.numeric(data_copernicus_temperature_shallow_2022_23$RUNA5), k = 7, align = "left", fill = NA)
  data_copernicus_temperature_shallow_2022_23$RUNA9 <- zoo::rollmean(as.numeric(data_copernicus_temperature_shallow_2022_23$RUNA9), k = 7, align = "left", fill = NA)
  
  path_to_data_copernicus_shallow_2022_23 <- "data/derived-data/data_copernicus_temperature_shallow_2022_23"
  
  write.csv(data_copernicus_temperature_shallow_2022_23, path_to_data_copernicus_shallow_2022_23)
  
  
  #### import raster 3 ####
  
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE, row.names = 1)
  data_gps_RODA <- data_gps[c(4:6),] 
  # data_gps_RODA[1,1] <- -19.662
  # data_gps_RODA[1,2] <- 63.5
  data_gps_RODA[2,1] <- -19.620
  data_gps_RODA[2,2] <- 63.420
  data_gps_RODA[3,1] <- -19.593
  data_gps_RODA[3,2] <- 63.533
  
  
  dat.ncfd <- ncdf4::nc_open(data_temp_copernicus_rodrigues)

  attributes(dat.ncfd$var)
  attributes(dat.ncfd$dim)
  lat <- ncdf4::ncvar_get(dat.ncfd, "latitude")
  long <- ncdf4::ncvar_get(dat.ncfd, "longitude")
  depth <- ncdf4::ncvar_get(dat.ncfd, "depth")
  time <- ncdf4::ncvar_get(dat.ncfd, "time")
  head(time)
  tunits <- ncatt_get(dat.ncfd, "time", "units") #check units
  nt <- dim(time)
  lswt_array <- ncvar_get(dat.ncfd, "thetao")
  dim(lswt_array)
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="GMT")
  range(time_obs)
  length(time_obs)


  lonlattime <- as.matrix(expand.grid(long,lat,time_obs))
  lswt_vec_long <- as.vector(lswt_array)
  lswt_obs <- data.frame(cbind(lonlattime, lswt_vec_long))
  colnames(lswt_obs) <- c("x","y","z","l")
  lswt_obs$l <- as.numeric(lswt_obs$l)
  lswt_obs_clean <- na.omit(lswt_obs)


  n <- length(levels(as.factor(lswt_obs_clean$z)))

  ####curve RUNA1####


  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  # i =1
  
  for (i in 1:n) {

    sub <- subset(lswt_obs_clean, z == DATE[i])

    moy <- aggregate(l~y+x, sub, mean)

    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)

    rast.temp <- raster::extract(rast.moy, data_gps_RODA[1,], method = 'simple')

    temp.moy[i] <- rast.temp[1]

    timeline[i] <- DATE[i]

  }

  T.RUNA1 <- as.data.frame(cbind(temp.moy, timeline))

  ####curve RUNA5####
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL

  for (i in 1:n) {

    sub <- subset(lswt_obs_clean, z == DATE[i])

    moy <- aggregate(l~y+x, sub, mean)

    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)

    rast.temp <- raster::extract(rast.moy, data_gps_RODA[2,], method = 'simple')

    temp.moy[i] <- rast.temp[1]

    timeline[i] <- DATE[i]

  }

  T.RUNA5 <- as.data.frame(cbind(temp.moy, timeline))

  ####curve RUNA9####
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL

  for (i in 1:n) {

    sub <- subset(lswt_obs_clean, z == DATE[i])

    moy <- aggregate(l~y+x, sub, mean)

    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)

    rast.temp <- raster::extract(rast.moy, data_gps_RODA[3,], method = 'simple')

    temp.moy[i] <- rast.temp[1]

    timeline[i] <- DATE[i]

  }

  T.RUNA9 <- as.data.frame(cbind(temp.moy, timeline))

  data_copernicus_temperature_rodrigues_2016_2017 <- as.data.frame(cbind(T.RUNA1$timeline, T.RUNA1$temp.moy, T.RUNA5$temp.moy, T.RUNA9$temp.moy))
  colnames(data_copernicus_temperature_rodrigues_2016_2017) <- c("date", "RODA1", "RODA2", "RODA3")
  
  data_copernicus_temperature_rodrigues_2016_2017$RODA1 <- zoo::rollmean(as.numeric(data_copernicus_temperature_rodrigues_2016_2017$RODA1), k = 7, align = "left", fill = NA)
  data_copernicus_temperature_rodrigues_2016_2017$RODA2 <- zoo::rollmean(as.numeric(data_copernicus_temperature_rodrigues_2016_2017$RODA2), k = 7, align = "left", fill = NA)
  data_copernicus_temperature_rodrigues_2016_2017$RODA3 <- zoo::rollmean(as.numeric(data_copernicus_temperature_rodrigues_2016_2017$RODA3), k = 7, align = "left", fill = NA)
  
  path_to_data_copernicus_rodrigues_2016_2017 <- "data/derived-data/data_copernicus_temperature_rodrigues_2016_2017"
 
  
  
  write.csv(data_copernicus_temperature_rodrigues_2016_2017, path_to_data_copernicus_rodrigues_2016_2017)
  
  
  #### import raster 4 ####
  
  data_gps <- read.csv(gps_sites, sep = ";", dec = ",", header = TRUE, row.names = 1)
  data_gps_RODA <- data_gps[c(4:6),] 
  
  dat.ncfd <- ncdf4::nc_open(data_temp_copernicus_rodrigues_2)
  
  attributes(dat.ncfd$var)
  attributes(dat.ncfd$dim)
  lat <- ncdf4::ncvar_get(dat.ncfd, "latitude")
  long <- ncdf4::ncvar_get(dat.ncfd, "longitude")
  depth <- ncdf4::ncvar_get(dat.ncfd, "depth")
  time <- ncdf4::ncvar_get(dat.ncfd, "time")
  head(time)
  tunits <- ncatt_get(dat.ncfd, "time", "units") #check units
  nt <- dim(time)
  lswt_array <- ncvar_get(dat.ncfd, "thetao")
  dim(lswt_array)
  time_obs <- as.POSIXct(time, origin = "1970-01-01", tz="GMT")
  range(time_obs)
  length(time_obs)
  
  
  lonlattime <- as.matrix(expand.grid(long,lat,time_obs))
  lswt_vec_long <- as.vector(lswt_array)
  lswt_obs <- data.frame(cbind(lonlattime, lswt_vec_long))
  colnames(lswt_obs) <- c("x","y","z","l")
  lswt_obs$l <- as.numeric(lswt_obs$l)
  lswt_obs_clean <- na.omit(lswt_obs)
  
  
  n <- length(levels(as.factor(lswt_obs_clean$z)))
  
  ####curve RUNA1####
  # data_gps_RODA[1,1] <- -19.662
  # data_gps_RODA[1,2] <- 63.5
  data_gps_RODA[2,1] <- -19.620
  data_gps_RODA[2,2] <- 63.420
  data_gps_RODA[3,1] <- -19.593
  data_gps_RODA[3,2] <- 63.533
  
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  # i =1
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~y+x, sub, mean)
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RODA[1,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]
    
  }
  
  T.RUNA1 <- as.data.frame(cbind(temp.moy, timeline))
  
  ####curve RUNA5####
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~y+x, sub, mean)
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RODA[2,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]
    
  }
  
  T.RUNA5 <- as.data.frame(cbind(temp.moy, timeline))
  
  ####curve RUNA9####
  DATE <- levels(as.factor(lswt_obs_clean$z))
  temp.moy <- NULL
  timeline <- NULL
  
  for (i in 1:n) {
    
    sub <- subset(lswt_obs_clean, z == DATE[i])
    
    moy <- aggregate(l~y+x, sub, mean)
    
    rast.moy <- raster::rasterFromXYZ(moy, digits = 3)
    
    rast.temp <- raster::extract(rast.moy, data_gps_RODA[3,], method = 'simple')
    
    temp.moy[i] <- rast.temp[1]
    
    timeline[i] <- DATE[i]
    
  }
  
  T.RUNA9 <- as.data.frame(cbind(temp.moy, timeline))
  
  data_copernicus_temperature_rodrigues_2022_2023 <- as.data.frame(cbind(T.RUNA1$timeline, T.RUNA1$temp.moy, T.RUNA5$temp.moy, T.RUNA9$temp.moy))
  colnames(data_copernicus_temperature_rodrigues_2022_2023) <- c("date", "RODA1", "RODA2", "RODA3")
  
  data_copernicus_temperature_rodrigues_2022_2023$RODA1 <- zoo::rollmean(as.numeric(data_copernicus_temperature_rodrigues_2022_2023$RODA1), k = 7, align = "left", fill = NA)
  data_copernicus_temperature_rodrigues_2022_2023$RODA2 <- zoo::rollmean(as.numeric(data_copernicus_temperature_rodrigues_2022_2023$RODA2), k = 7, align = "left", fill = NA)
  data_copernicus_temperature_rodrigues_2022_2023$RODA3 <- zoo::rollmean(as.numeric(data_copernicus_temperature_rodrigues_2022_2023$RODA3), k = 7, align = "left", fill = NA)
  
  path_to_data_copernicus_rodrigues_2022_2023 <- "data/derived-data/data_copernicus_temperature_rodrigues_2022_2023"
  
  write.csv(data_copernicus_temperature_rodrigues_2022_2023, path_to_data_copernicus_rodrigues_2022_2023)
  
  path_to_data_copernicus <- c(path_to_data_copernicus_shallow_2019_20, path_to_data_copernicus_shallow_2022_23, path_to_data_copernicus_rodrigues_2016_2017, path_to_data_copernicus_rodrigues_2022_2023)
  names(path_to_data_copernicus) <- c("data_copernicus_shallow_2019_20", "data_copernicus_shallow_2022_23", "data_copernicus_rodrigues_2016_2017", "data_copernicus_rodrigues_2022_2023")
  
  
return(path_to_data_copernicus)

}