sites <- readRDS("sites_data.rds")

# select all countries which have 5 or more sites

big_enough_countries <- names(which(table(sites$Country) >=5))

rescaling <- function(x) {(x-min(x))/(max(x)-min(x))}

rescaling_abs <- function(x) {abs(x)/max(abs(x))}

## getting the ordered parameter names ##

parameters = readRDS("Parameter_list.rds")
drivernames  = paste0("run_",c("co2","ndep","insol","temp","ph","prec"),"_change")
parameters_pic_abi = c(as.character(parameters$Pic_abi$NameRLPJ), drivernames)
parameters_fag_syl = c(as.character(parameters$Fag_syl$NameRLPJ), drivernames)
parameters_pin_syl = c(as.character(parameters$Pin_syl$NameRLPJ), drivernames)

results =  readRDS(paste0("./LPJrunTest/Results/Pic_abi_0.25_42.25.rds"))
parameternames_ordered_pic_abi = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_pic_abi)])

results =  readRDS(paste0("./LPJrunTest/Results/Fag_syl_0.25_42.25.rds"))
parameternames_ordered_fag_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_fag_syl)])

results =  readRDS(paste0("./LPJrunTest/Results/Pin_syl_0.25_42.25.rds"))
parameternames_ordered_pin_syl = names(results$parameters[[1]][which(names(results$parameters[[1]]) %in% parameters_pin_syl)])

results_mixed = readRDS("./LPJrunTest/Results/mapping_mixed.rds")
parameternames_ordered_mixed2 = results_mixed$parameternames

parameternames_ordered_pic_abi2 = substring(parameternames_ordered_pic_abi,regexpr("_", parameternames_ordered_pic_abi) + 1)
parameternames_ordered_fag_syl2 = substring(parameternames_ordered_fag_syl,regexpr("_", parameternames_ordered_fag_syl) + 1)
parameternames_ordered_pin_syl2 = substring(parameternames_ordered_pin_syl,regexpr("_", parameternames_ordered_pin_syl) + 1)


noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  parameternames_ordered_pic_abi2 = gsub(noisy_things[i],"",parameternames_ordered_pic_abi2)
  parameternames_ordered_fag_syl2 = gsub(noisy_things[i],"",parameternames_ordered_fag_syl2)
  parameternames_ordered_pin_syl2 = gsub(noisy_things[i],"",parameternames_ordered_pin_syl2)
  parameternames_ordered_mixed2 = gsub(noisy_things[i],"",parameternames_ordered_mixed2)
}

ordering_fag_syl = match(parameternames_ordered_pic_abi2,parameternames_ordered_fag_syl2)
ordering_pin_syl = match(parameternames_ordered_pic_abi2,parameternames_ordered_pin_syl2)
ordering_mixed = match(parameternames_ordered_pic_abi2, parameternames_ordered_mixed2)

parameters = readRDS("Parameter_list.rds")
drivers = c("ph","co2","prec","temp","insol","ndep")

variablenames = c(as.character(parameters$Pic_abi[,"NameRLPJ"]),drivers)
variablenames2 = substring(variablenames,regexpr("_", variablenames) + 1)
noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  variablenames2 = gsub(noisy_things[i],"",variablenames2)
}

grouping = read.csv("./Grouping_Fag_syl.csv", header = T, sep = ";")
grouping = c(as.character(grouping[,1]), rep("Drivers",6))

order_grouping = match(parameternames_ordered_pic_abi2, variablenames2)
grouping2 = grouping[order_grouping]


### Loading the results of the linear regressions ###

effects_Pic_abi = readRDS("LPJrunTest/Results/Pic_abi_effects.rds")

effects_Fag_syl = readRDS("LPJrunTest/Results/Fag_syl_effects.rds")

effects_Pin_syl = readRDS("LPJrunTest/Results/Pin_syl_effects.rds")

effects_mixed = readRDS("LPJrunTest/Results/Mixed_effects.rds")


###

library(rgeos)
library(rworldmap)
library(sp)
#library(sf)
library(tidyverse)
library(ggplot2)
library(rgdal)
library(scatterpie)


# get world map
wmap <- getMap(resolution="high")

# get centroids
centroids <- gCentroid(wmap, byid=TRUE)
aoi_boundary_HARV <- readOGR("./EnSv8/ens_v8.shp")
# get a data.frame with centroids
orig_df = as.data.frame(centroids)
df <- as.data.frame(centroids)
coordinates(df) <- c("x", "y")
proj4string(df) <- CRS("+init=epsg:4326")

df_new <- spTransform(df, proj4string(aoi_boundary_HARV))

countries = c("Finland","France","Germany","Italy","Norway","Poland","Romania","Russia","Spain","Sweden","United Kingdom","Ukraine")

big_countries_data_set = data.frame("country" = countries, "lon" = df_new[match(countries,rownames(orig_df)),]$x,
                                    "lat" = df_new[match(countries,rownames(orig_df)),]$y)



result_matrix = as.data.frame(matrix(nrow = 12, ncol = 10))
rownames(result_matrix) = countries



counter = 1
for(country in big_enough_countries){
  print(country)
  sites_per_country = which(sites$Country == country)
  mean_effects_mixed_lai = effects_mixed[["lai"]][["complete"]][sites_per_country,]
  mapped_effects_mixed_lai =  matrix(nrow = nrow(mean_effects_mixed_lai), ncol = 38)
  for(parameter in 1:length(results_mixed$position_mapping)){
    for(site in 1:nrow(mean_effects_mixed_lai)){
      mapped_effects_mixed_lai[site,parameter] = sum(mean_effects_mixed_lai[site,results_mixed$position_mapping[[parameter]]]*
        results_mixed$weight_mapping[[parameter]])
    }
  }
  mapped_effects_mixed_lai = mapped_effects_mixed_lai[complete.cases(mapped_effects_mixed_lai),]
  mapped_effects_mixed_lai = t(apply(mapped_effects_mixed_lai,1,rescaling_abs))

  main_effects_mixed_lai = mapped_effects_mixed_lai[ordering_mixed]

  auxiliary_df = data.frame("group" = grouping2, "effect" = main_effects_mixed_lai)

  results = aggregate(effect ~ group, auxiliary_df, sum)
  result_matrix[counter,] = results$effect
  counter = counter +1
}

colnames(result_matrix) = results$group

mixed_lai_df <- cbind(big_countries_data_set,result_matrix)
mixed_lai_df = mixed_lai_df[-8,]



aoi_boundary_HARV <- st_read("./EnSv8/ens_v8.shp")

library(viridis)
library(ggnewscale)
library(RColorBrewer)
myColors = magma(10)


community_composition_plot <- ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = magma(13)) +
  new_scale("fill") +
  coord_sf() +
  geom_scatterpie(aes(x=lon, y=lat, group = country, r = 100000),
                  data = mixed_lai_df, cols = colnames(mixed_lai_df[4:ncol(mixed_lai_df)]),
                  color = NA)  + scale_fill_brewer(palette = "Paired") + ggtitle("Community composition")



