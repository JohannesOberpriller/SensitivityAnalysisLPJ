library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(raster)
set.seed(123)

## Defining helper functions to determine the distance to the sea and the shortest distance the next grid

shortest_distances <- function(x, y){
  distances = vector(mode = "numeric", length = length(x))
  for(i in 1:length(x)){
    dists = c()
    for(j in (i+1):(length(x)-1)){
      dists[j] = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))
    }
    distances[i] = min(dists, na.rm = T)
  }
  return(distances)
}


distance_to_sea <- function(lat,lon){
  by_sea = vector(mode = "logical", length = length(lat))
  for(i in 1:length(lat)){
    Data_full_lat = Data_full[which(abs(Data_full$lat - lat[i]) <0.4),]
    Data_full_lon_lat = Data_full_lat[which(abs(Data_full_lat$lon - lon[i])< 0.4),]
    neighbors = nrow(Data_full_lon_lat)
    if(neighbors == 9) {
      by_sea[i] = T
    }
    else{
      by_sea[i] = F
    }
  }
  return(by_sea)
}

# loading the data and assigning the longitudinal and latitudinal coordinates

ncin_qq <-ncdf4::nc_open("/Site_selection_data/qq_ens_mean_0.25deg_reg_2011-2019_v20.0e.nc")

lon <- ncvar_get(ncin_qq,"longitude")
lat <- ncvar_get(ncin_qq,"latitude")

# calculating the sum over the time

qq <- ncvar_get(ncin_qq, "qq")
qq_summed = apply(FUN = mean, X = qq, MARGIN = c(1,2),na.rm = T)

## removing data form RAM issues

rm(qq)
# loading the data

ncin_tg = ncdf4::nc_open("/Site_selection_data/tg_ens_mean_0.25deg_reg_2011-2019_v20.0e.nc")

# calculating the sum over the time

tg <- ncvar_get(ncin_tg, "tg")
tg_summed = apply(FUN = mean, X = tg, MARGIN = c(1,2),na.rm = T)

## removing data form RAM issues

rm(tg)

# loading data
ncin_rr <- nc_open("/Site_selection_data/rr_ens_mean_0.25deg_reg_2011-2019_v20.0e.nc")

# calculating the sum over the time

rr <- ncvar_get(ncin_rr,"rr")
rr_summed = apply(FUN = mean, X = rr, MARGIN = c(1,2), na.rm = T)
rr_summed = 365.25*rr_summed

## removing data form RAM issues
rm(rr)

##generating data frame to work with

lonlat <- as.matrix(expand.grid(lon,lat))

Data <- as.data.frame(cbind(lonlat,as.vector(qq_summed),as.vector(tg_summed),as.vector(rr_summed)))

colnames(Data) = c("lon","lat","radiation","temperature","precipitation")

## excluding low latitudes and high longitudinal coordinates

Data_full = Data[complete.cases(Data),]
Data_full = Data_full[-which(Data_full$lat<37),]
Data_full = Data_full[-which(Data_full$lon >24.5)]

## excluding sites near the sea

sea_sites = distance_to_sea(lat = Data_full$lat, lon = Data_full$lon)

Data_full_without_sea = Data_full[which(sea_sites == T),]

library(splitstackshape)

library(Hmisc)
library(fields)
library(colorRamps)

## Applying random stratified sampling to the data

groupstemp <- cut2(Data_full_without_sea$temperature, g=5)
groupsprec <- cut2(Data_full_without_sea$precipitation, g=5)
grouprad <- cut2(Data_full_without_sea$radiation, g=5)
groupedlon <- cut(Data_full_without_sea$lon, breaks =20)
groupedlat <-  cut(Data_full_without_sea$lat, breaks =20)

Data_full_grouped = as.data.frame(cbind(Data_full_without_sea,groupstemp,groupsprec, grouprad, groupedlat, groupedlon))



#sample = stratified(Data_full_grouped,c("groupstemp","groupsprec","grouprad", "groupedlat", "groupedlon"), size = 2)
sample = stratified(Data_full_grouped,c("groupstemp","groupsprec","grouprad", "groupedlat", "groupedlon"), size = 5 )
sample2 = sample(1:nrow(sample), size = 500)
samples = sample[sample2,]


## calulating distance between sites and removing if, they are to close to each other

shorest_distances <-shortest_distances(samples$lon, samples$lat)

big_enough = which(shorest_distances > 1.)

samples = samples[big_enough,]

## Making plots


colors = colorRampPalette(c("blue","lightblue","orange","red"), bias = 1.5)(64)

colors2 = colorRampPalette(c("red","orange","lightblue","blue"), bias = 2.)(64)

image.plot(lon,lat,qq_summed, xlim = c(-11,25), ylim = c(36,72), col = colors,
           xaxt ='n', xlab = '', yaxt = 'n', ylab = '')
mtext("a) Mean daily solar radiation [W/m^2]", side =3, cex = 1.)
points(x=samples$lon, y=samples$lat,add = T, pch = 19)

p2 <- image.plot(lon,lat,tg_summed,xlim = c(-11,25), ylim = c(36,72), col = colors,
                 xaxt ='n', xlab = '', yaxt = 'n', ylab = '')
mtext("b) Mean annual temperature [C°]", side =3, cex = 1)
points(x=samples$lon, y=samples$lat,add = T, pch = 19)

p3 <- image.plot(lon,lat,rr_summed,xlim = c(-11,25), ylim = c(36,72), col = colors,
                 xaxt ='n', xlab = '', yaxt = 'n', ylab = '')
mtext("c) Mean yearly precipitation [l/m^2]", cex = 1.)
points(x=samples$lon, y=samples$lat,add = T, pch = 19)


## Saving the stuff
library(NOOA)


save_data_frame = samples
save_data_frame$groupedlat <- NULL
save_data_frame$groupedlon <- NULL
save_data_frame$grouprad <- NULL
save_data_frame$groupsprec <- NULL
save_data_frame$groupstemp <- NULL

save_data_frame$countries = countries


for(i in 1:length(save_data_frame$lon)){
  if((save_data_frame$lon[i]-floor(save_data_frame$lon[i]))<0.5){
    save_data_frame$lon[i] = floor(save_data_frame$lon[i])+0.25
  }
  else{
    save_data_frame$lon[i] = floor(save_data_frame$lon[i])+0.75
  }
  if((save_data_frame$lat[i]-floor(save_data_frame$lat[i]))<0.5){
    save_data_frame$lat[i] = floor(save_data_frame$lat[i])+0.25
  }
  else{
    save_data_frame$lat[i] = floor(save_data_frame$lat[i])+0.75
  }
  save_data_frame$radiation[i] = round(save_data_frame$radiation[i],1)
  save_data_frame$temperature[i] = round(save_data_frame$temperature[i],1)
  save_data_frame$precipitation[i] = round(save_data_frame$precipitation[i],1)
}


save_data_frame$countries = coords2country(save_data_frame$lat,save_data_frame$lon)

library(fields)
par(mfrow = c(1,3))

image.plot(lon,lat,qq_summed, xlim = c(-11,25), ylim = c(36,72), col = colors,
           xaxt ='n', xlab = '', yaxt = 'n', ylab = '')
mtext("a) Mean daily solar radiation [W/m^2]", side =3, cex = 1.)
points(x=save_data_frame$lon, y=save_data_frame$lat,add = T, pch = 19)

p2 <- image.plot(lon,lat,tg_summed,xlim = c(-11,25), ylim = c(36,72), col = colors,
                 xaxt ='n', xlab = '', yaxt = 'n', ylab = '')
mtext("b) Mean annual temperature [C°]", side =3, cex = 1)
points(x=save_data_frame$lon, y=save_data_frame$lat,add = T, pch = 19)

p3 <- image.plot(lon,lat,rr_summed,xlim = c(-11,25), ylim = c(36,72), col = colors2,
                 xaxt ='n', xlab = '', yaxt = 'n', ylab = '')
mtext("c) Mean yearly precipitation [l/m^2]", cex = 1.)
points(x=save_data_frame$lon, y=save_data_frame$lat,add = T, pch = 19)


save_data_frame$Site = 1:nrow(save_data_frame)
save_data_frame = save_data_frame[,c(7,6,1,2,3,4,5)]
colnames(save_data_frame) = c("Site","Country", "Longitudinal", "Latitudinal","Radiation [W/m^2]", "Temperature [C°]", "Precipitation [l/m^2]")
saveRDS(save_data_frame, file = "EnvironmentalData/sites_data.rds")





