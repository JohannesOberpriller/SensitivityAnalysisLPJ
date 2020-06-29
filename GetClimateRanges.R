library(ncdf4)

ncpath <- "/home/johannes/crop_ncep/data/Johannes/Data/"
ncname <- "rsds_bced_1960_1999_ipsl-cm5a-lr_hist_rcp8p5_1901-2099_noleap_monmean_chunked_rew"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
gc()
# ncin <- nc_open(ncfname)
# print(ncin)
# variable <- ncvar_get(ncin,"rsdsAdjust")
# dim(variable)
# variable$lon
# 
# lon <- ncvar_get(ncin,"lon")
# lat <- ncvar_get(ncin,"lat")
# time <- ncvar_get(ncin,"time")
# tunits<-ncatt_get(ncin,"time",attname="units")
# tustr<-strsplit(tunits$value, " ")
# dates<-as.Date(time,origin=unlist(tustr)[3])
# 
# which(substring(dates,1,4)>2010)
# which(lon == "22.25")
# which(lat == "22.25")
# 
# max(variable[which(lon == "22.25"),which(lat == "22.25"),which(substring(dates,1,4)>2015)])
# min(variable[which(lon == "22.25"),which(lat == "22.25"),which(substring(dates,1,4)>2015)])

sites <- readRDS("/home/johannes/Documents/PhD/Code/LPJ_Guess_KIT/rLPJGUESS/rLPJGUESS/Testcode/sites_data.rds")


ncvariable = "rsdsAdjust"
get_climate_extremes <- function(sites,ncfname,ncvariable, time ){
  
  ncin <- nc_open(ncfname)
  variable <- ncvar_get(ncin, ncvariable)
  lon <- ncvar_get(ncin,"lon")
  lat <- ncvar_get(ncin,"lat")
  time <- ncvar_get(ncin,"time")
  tunits<-ncatt_get(ncin,"time",attname="units")
  tustr<-strsplit(tunits$value, " ")
  dates<-as.Date(time,origin=unlist(tustr)[3])
  
  maximum = vector(length = nrow(sites))
  minimum = vector(length = nrow(sites))
  for(i in 1:nrow(sites)){
    maximum[i] = max(variable[which(lon == sites$Longitudinal[i]),which(lat == sites$Latitudinal[i]),which(substring(dates,1,4)>time)])
    minimum[i] = min(variable[which(lon == sites$Longitudinal[i]),which(lat == sites$Latitudinal[i]),which(substring(dates,1,4)>time)])
  }
  gc()
  return(list("maximum" = maximum, "minimum"= minimum))
}

hans = get_climate_extremes(sites, ncfname, ncvariable = ncvariable, time = 2010) 

read.filename <- file("/home/johannes/crop_ncep/data/Johannes/Data/GlobalNitrogenDepositionRCP85.bin", "rb")

readBin("/home/johannes/crop_ncep/data/Johannes/Data/GlobalNitrogenDepositionRCP85.bin", double(), n = 200)
readBin(read.filename, character(),  n = 2, endian = "little")
ncin <- nc_open(ncfname)
variable <- ncvar_get(ncin, ncvariable)
