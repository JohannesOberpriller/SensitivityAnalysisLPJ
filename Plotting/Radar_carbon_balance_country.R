#load the site data
sites <- readRDS("EnvironmentalData/sites_data.rds")

# function to remap the parameters of mixed results
remap_paramters <- function(parameters, weighting_scheme, position_scheme){
  mapped_parameters = matrix(ncol = length(weighting_scheme),
                             nrow = nrow(parameters))
  for(site in 1:nrow(parameters)){
    for(parameter in 1:length(position_scheme)){
      mapped_parameters[site,parameter] = sum(abs(parameters[site,position_scheme[[parameter]]])*
                                                weighting_scheme[[parameter]])
    }
  }
  return(mapped_parameters)
}


# select all countries which have 5 or more sites

big_enough_countries <- names(which(table(sites$Country) >=5))

rescaling <- function(x) {(x-min(x))/(max(x)-min(x))}

rescaling_abs <- function(x) {abs(x)/max(abs(x))}

## getting the ordered parameter names ##

results_mixed = readRDS("./LPJrunTest/Results/mapping_mixed.rds")

# load the the parameters and clear the names to be usable
parameters = readRDS("ParameterMetaData/Parameter_list.rds")
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

parameternames_ordered_pic_abi2 = substring(parameternames_ordered_pic_abi,regexpr("_", parameternames_ordered_pic_abi) + 1)
parameternames_ordered_fag_syl2 = substring(parameternames_ordered_fag_syl,regexpr("_", parameternames_ordered_fag_syl) + 1)
parameternames_ordered_pin_syl2 = substring(parameternames_ordered_pin_syl,regexpr("_", parameternames_ordered_pin_syl) + 1)
parameternames_ordered_mixed2 = results_mixed$parameternames

noisy_things = c("intolerant","tolerant","tree","abi","change","_","syl","leaved")
for(i in 1:length(noisy_things)){
  parameternames_ordered_pic_abi2 = gsub(noisy_things[i],"",parameternames_ordered_pic_abi2)
  parameternames_ordered_fag_syl2 = gsub(noisy_things[i],"",parameternames_ordered_fag_syl2)
  parameternames_ordered_pin_syl2 = gsub(noisy_things[i],"",parameternames_ordered_pin_syl2)
  parameternames_ordered_mixed2 = gsub(noisy_things[i],"",parameternames_ordered_mixed2)
}

ordering_fag_syl = match(parameternames_ordered_pic_abi2,parameternames_ordered_fag_syl2)
ordering_pin_syl = match(parameternames_ordered_pic_abi2,parameternames_ordered_pin_syl2)
ordering_mixed = match(parameternames_ordered_pic_abi2,parameternames_ordered_mixed2)

## load the list of parameters

parameters = readRDS("ParameterMetaData/Parameter_list.rds")
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

effects_mixed = readRDS("LPJrunTest/Results/mixed_effects.rds")


###

library(rgeos)
library(rworldmap)

library(sp)
library(sf)
library(tidyverse)
library(ggplot2)
library(rgdal)


# get world map
wmap <- getMap(resolution="low")
aoi_boundary_HARV <- readOGR("./EnSv8/ens_v8.shp")
# get centroids
centroids <- gCentroid(wmap, byid=TRUE)

# get a data.frame with centroids
orig_df = as.data.frame(centroids)
df <- as.data.frame(centroids)
coordinates(df) <- c("x", "y")
proj4string(df) <- CRS("+init=epsg:4326")
df_new <- spTransform(df, proj4string(aoi_boundary_HARV))

countries = c("Finland","France","Germany","Italy","Norway","Poland","Romania","Russia","Spain","Sweden","United Kingdom","Ukraine")

big_countries_data_set = data.frame("country" = countries, "lon" = df_new[match(countries,rownames(orig_df)),]$x,
                                    "lat" = df_new[match(countries,rownames(orig_df)),]$y)



#Going over the settings and doing the same for each setting


## Pic abi and cflux

# set up matrix to store stuff
result_matrix = as.data.frame(matrix(nrow = 12, ncol = 7))
rownames(result_matrix) = countries

counter = 1

# get the data results and only use sites where the tree grows
mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)
sites_pic_abi = sites[growth_sites_Pic_abi,]

# looping over each country and aggregating effects per process
for(country in big_enough_countries){
  sites_per_country = which(sites_pic_abi$Country == country)
  if(length(sites_per_country) >= 5){
    # calculating effects on the correct site
    mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux"]][["complete"]][growth_sites_Pic_abi,][sites_per_country,]
    # building the average per site
    main_effects_Pic_abi_cflux = apply(abs(mean_effects_Pic_abi_cflux[,1:38]),2,mean)
    # aggregate data
    main_effects_Pic_abi_cflux_weight = main_effects_Pic_abi_cflux/sum(abs(main_effects_Pic_abi_cflux))
    # dataframe to store data
    auxiliary_df = data.frame("group" = grouping2, "effect" = main_effects_Pic_abi_cflux_weight)

    results = aggregate(effect ~ group, auxiliary_df, sum)
    result_matrix[counter,] = results$effect
  }
  counter = counter +1
}

colnames(result_matrix) = results$group


pic_abi_cflux_df <- cbind(big_countries_data_set,result_matrix)
pic_abi_cflux_df = pic_abi_cflux_df[-which(pic_abi_cflux_df$country == "Russia"),]



result_matrix = as.data.frame(matrix(nrow = 12, ncol = 7))
rownames(result_matrix) = countries

counter = 1

mean_effects_Fag_syl_cpool = effects_Fag_syl[["cpool"]][["complete"]]
growth_sites_Fag_syl = which(mean_effects_Fag_syl_cpool[,39]/200 > 2)
sites_Fag_syl = sites[growth_sites_Fag_syl,]

for(country in big_enough_countries){
  sites_per_country = which(sites_Fag_syl$Country == country)
  if(length(sites_per_country) >= 5){
    mean_effects_Fag_syl_cflux = effects_Fag_syl[["cflux"]][["complete"]][growth_sites_Fag_syl,][sites_per_country,]

    main_effects_Fag_syl_cflux = apply(abs(mean_effects_Fag_syl_cflux[,1:38]),2,mean)[ordering_fag_syl]

    main_effects_Fag_syl_cflux_weight = main_effects_Fag_syl_cflux/sum(abs(main_effects_Fag_syl_cflux))

    auxiliary_df = data.frame("group" = grouping2, "effect" = main_effects_Fag_syl_cflux_weight)

    results = aggregate(effect ~ group, auxiliary_df, sum)
    result_matrix[counter,] = results$effect
  }
  counter = counter +1
}

colnames(result_matrix) = results$group


fag_syl_cflux_df <- cbind(big_countries_data_set,result_matrix)
fag_syl_cflux_df = fag_syl_cflux_df[-which(fag_syl_cflux_df$country == "Russia"),]



result_matrix = as.data.frame(matrix(nrow = 12, ncol = 7))
rownames(result_matrix) = countries

counter = 1

mean_effects_Pin_syl_cpool = effects_Pin_syl[["cpool"]][["complete"]]
growth_sites_Pin_syl = which(mean_effects_Pin_syl_cpool[,39]/200 > 2)
sites_Pin_syl = sites[growth_sites_Pin_syl,]

for(country in big_enough_countries){
  sites_per_country = which(sites_Pin_syl$Country == country)
  if(length(sites_per_country) >= 5){
    mean_effects_Pin_syl_cflux = effects_Pin_syl[["cflux"]][["complete"]][growth_sites_Pin_syl,][sites_per_country,]

    main_effects_Pin_syl_cflux = apply(abs(mean_effects_Pin_syl_cflux[,1:38]),2,mean)[ordering_pin_syl]

    main_effects_Pin_syl_cflux_weight = main_effects_Pin_syl_cflux/sum(abs(main_effects_Pin_syl_cflux))

    auxiliary_df = data.frame("group" = grouping2, "effect" = main_effects_Pin_syl_cflux_weight)

    results = aggregate(effect ~ group, auxiliary_df, sum)
    result_matrix[counter,] = results$effect
  }
  counter = counter +1
}

colnames(result_matrix) = results$group


pin_syl_cflux_df <- cbind(big_countries_data_set,result_matrix)
pin_syl_cflux_df = pin_syl_cflux_df[-which(pin_syl_cflux_df$country == "Russia"),]



result_matrix = as.data.frame(matrix(nrow = 12, ncol = 7))
rownames(result_matrix) = countries

counter = 1

mean_effects_mixed_cpool = effects_mixed[["cpool"]][["complete"]]
growth_sites_mixed = which(mean_effects_mixed_cpool[,74]/200 > 2)
sites_mixed = sites[growth_sites_mixed,]

for(country in big_enough_countries){
  sites_per_country = which(sites_mixed$Country == country)
  if(length(sites_per_country) >= 5){
    mean_effects_mixed_cflux = effects_mixed[["cflux"]][["complete"]][growth_sites_mixed,][sites_per_country,]


    remaped_effects_cpool_abs = remap_paramters(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)

    main_effects_mixed_cflux = apply(abs(remaped_effects_cpool_abs),2,mean)[ordering_mixed]

    main_effects_mixed_cflux_weight = main_effects_mixed_cflux/sum(abs(main_effects_mixed_cflux))

    auxiliary_df = data.frame("group" = grouping2, "effect" = main_effects_mixed_cflux_weight)

    results = aggregate(effect ~ group, auxiliary_df, sum)
    result_matrix[counter,] = results$effect
  }
  counter = counter +1
}

colnames(result_matrix) = results$group


mixed_cflux_df <- cbind(big_countries_data_set,result_matrix)
mixed_cflux_df = mixed_cflux_df[-which(mixed_cflux_df$country == "Russia"),]

result_matrix = as.data.frame(matrix(nrow = 12, ncol = 7))
rownames(result_matrix) = countries



counter = 1
for(country in big_enough_countries){
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
  main_effects_mixed_lai = apply(abs(mapped_effects_mixed_lai),2,mean)[ordering_mixed]

  main_effects_mixed_lai = main_effects_mixed_lai/sum(abs(main_effects_mixed_lai))

  auxiliary_df = data.frame("group" = grouping2, "effect" = main_effects_mixed_lai)

  results = aggregate(effect ~ group, auxiliary_df, sum)
  result_matrix[counter,] = results$effect
  counter = counter +1
}

colnames(result_matrix) = results$group

mixed_lai_df <- cbind(big_countries_data_set,result_matrix)
mixed_lai_df = mixed_lai_df[-8,]



colors_me = c("gold","burlywood1","red2","darkolivegreen2","darkgreen","chocolate4","blue2")
aoi_boundary_HARV <- st_read("./EnSv8/ens_v8.shp")

library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(scatterpie)
library(ggplot2)
library(fmsb)
library(ggradar)

colors_env = c("darkblue","blue3","chartreuse2","darkgreen","green3",
               "deepskyblue","green1","gold",'yellow',
               "yellow2","darkorange","cyan3","darkolivegreen3")

max(mixed_cflux_df[,4:10])
radar.testplot_list <-
  lapply(1:nrow(pic_abi_cflux_df), function(i) {

    gt_plot <- ggplotGrob(
      ggradar(rbind(cbind("Name" = "Pic. abi.",pic_abi_cflux_df[i,4:10]),
                    cbind("Name" = "Fag. syl.",fag_syl_cflux_df[i,4:10]),
                    cbind("Name" = "Pin. syl.",pin_syl_cflux_df[i,4:10]),
                    cbind("Name" = "Mixed",mixed_cflux_df[i,4:10])),
              plot.title = rownames(pic_abi_cflux_df)[i], grid.line.width = 0.5,
              plot.legend = F, grid.max = 0.4, grid.mid = 0.2,
              values.radar = c("0%", "20%", "40%"), group.point.size =3)
    )




  })

pdf("./Figures/radar_no_map.pdf", width = 30, height = 25)
grid.arrange(arrangeGrob(radar.testplot_list[[1]],radar.testplot_list[[2]],radar.testplot_list[[3]], ncol =3, nrow =1),
             arrangeGrob(radar.testplot_list[[4]],radar.testplot_list[[5]],radar.testplot_list[[6]], ncol =3, nrow =1),
             arrangeGrob(radar.testplot_list[[7]],radar.testplot_list[[8]],radar.testplot_list[[9]], ncol =3, nrow =1),
             arrangeGrob(ggplot()+theme_void(),radar.testplot_list[[10]],radar.testplot_list[[11]],ggplot()+theme_void(), ncol =4, nrow =1,
                         widths = c(0.5,1,1,0.5)),nrow = 4)
dev.off()





# ##### Carbon balance driver analysis #####
#
# ## Pic abi and cflux
#
# result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
# rownames(result_matrix) = countries
#
# counter = 1
#
# mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
# growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)
# sites_pic_abi = sites[growth_sites_Pic_abi,]
#
# for(country in big_enough_countries){
#   sites_per_country = which(sites_pic_abi$Country == country)
#   mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux"]][["complete"]][growth_sites_Pic_abi,][sites_per_country,]
#
#   main_effects_Pic_abi_cflux_driver = apply(abs(mean_effects_Pic_abi_cflux[,1:38]),2,mean)[which(grouping2 == "Drivers")]
#
#   main_effects_Pic_abi_cflux_weight_driver = main_effects_Pic_abi_cflux_driver/sum(abs(main_effects_Pic_abi_cflux_driver))
#
#   result_matrix[counter,] = main_effects_Pic_abi_cflux_weight_driver
#
#   counter = counter +1
# }
#
# colnames(result_matrix) = parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]
#
#
# pic_abi_cflux_driver_df <- cbind(big_countries_data_set,result_matrix)
# pic_abi_cflux_driver_df = pic_abi_cflux_driver_df[-which(pic_abi_cflux_driver_df$country == "Russia"),]
#
#
#
# result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
# rownames(result_matrix) = countries
#
# counter = 1
#
# mean_effects_fag_syl_cpool = effects_Fag_syl[["cpool"]][["complete"]]
# growth_sites_fag_syl = which(mean_effects_fag_syl_cpool[,39]/200 > 2)
# sites_fag_syl = sites[growth_sites_fag_syl,]
#
# for(country in big_enough_countries){
#   sites_per_country = which(sites_fag_syl$Country == country)
#   if(length(sites_per_country) >= 5){
#     mean_effects_fag_syl_cflux = effects_Fag_syl[["cflux"]][["complete"]][growth_sites_fag_syl,][sites_per_country,]
#
#     main_effects_fag_syl_cflux_driver = apply(abs(mean_effects_fag_syl_cflux[,1:38]),2,mean)[ordering_fag_syl][which(grouping2 == "Drivers")]
#
#     main_effects_fag_syl_cflux_weight_driver = main_effects_fag_syl_cflux_driver/sum(abs(main_effects_fag_syl_cflux_driver))
#
#     result_matrix[counter,] = main_effects_fag_syl_cflux_weight_driver
#   }
#   counter = counter +1
# }
#
# colnames(result_matrix) = parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]
#
#
# fag_syl_cflux_driver_df <- cbind(big_countries_data_set,result_matrix)
# fag_syl_cflux_driver_df = fag_syl_cflux_driver_df[-which(fag_syl_cflux_driver_df$country == "Russia"),]
#
#
#
# result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
# rownames(result_matrix) = countries
#
# counter = 1
#
# mean_effects_pin_syl_cpool = effects_Pin_syl[["cpool"]][["complete"]]
# growth_sites_pin_syl = which(mean_effects_pin_syl_cpool[,39]/200 > 2)
# sites_pin_syl = sites[growth_sites_pin_syl,]
#
# for(country in big_enough_countries){
#   sites_per_country = which(sites_pin_syl$Country == country)
#   if(length(sites_per_country) >= 5){
#     mean_effects_pin_syl_cflux = effects_Pin_syl[["cflux"]][["complete"]][growth_sites_pin_syl,][sites_per_country,]
#
#     main_effects_pin_syl_cflux_driver = apply(abs(mean_effects_pin_syl_cflux[,1:38]),2,mean)[ordering_pin_syl][which(grouping2 == "Drivers")]
#
#     main_effects_pin_syl_cflux_weight_driver = main_effects_pin_syl_cflux_driver/sum(abs(main_effects_pin_syl_cflux_driver))
#
#     result_matrix[counter,] = main_effects_pin_syl_cflux_weight_driver
#   }
#   counter = counter +1
# }
#
# colnames(result_matrix) = parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]
#
#
# pin_syl_cflux_driver_df <- cbind(big_countries_data_set,result_matrix)
# pin_syl_cflux_driver_df = pin_syl_cflux_driver_df[-which(pin_syl_cflux_driver_df$country == "Russia"),]
#
#
#
# result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
# rownames(result_matrix) = countries
#
# counter = 1
#
# mean_effects_mixed_cpool = effects_mixed[["cpool"]][["complete"]]
# growth_sites_mixed = which(mean_effects_mixed_cpool[,74]/200 > 2)
# sites_mixed = sites[growth_sites_mixed,]
#
# for(country in big_enough_countries){
#   print(sites_per_country)
#   sites_per_country = which(sites_mixed$Country == country)
#   if(length(sites_per_country) >= 5){
#     mean_effects_mixed_cflux = effects_mixed[["cflux"]][["complete"]][growth_sites_mixed,][sites_per_country,]
#
#     remaped_effects_cpool_abs = remap_paramters(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)
#
#     main_effects_mixed_cflux_driver = apply(abs(remaped_effects_cpool_abs),2,mean)[ordering_mixed][which(grouping2 == "Drivers")]
#
#     main_effects_mixed_cflux_driver_weight = main_effects_mixed_cflux_driver/sum(abs(main_effects_mixed_cflux_driver))
#
#     result_matrix[counter,] = main_effects_mixed_cflux_driver_weight
#   }
#   counter = counter +1
# }
#
# colnames(result_matrix) =  parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]
#
#
# mixed_cflux_driver_df <- cbind(big_countries_data_set,result_matrix)
# mixed_cflux_driver_df = mixed_cflux_driver_df[-which(mixed_cflux_driver_df$country == "Russia"),]
#
# result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
# rownames(result_matrix) = countries
#
#
#
# counter = 1
# for(country in big_enough_countries){
#   sites_per_country = which(sites$Country == country)
#   mean_effects_mixed_lai = effects_mixed[["lai"]][["complete"]][sites_per_country,]
#   mapped_effects_mixed_lai =  matrix(nrow = nrow(mean_effects_mixed_lai), ncol = 38)
#   for(parameter in 1:length(results_mixed$position_mapping)){
#     for(site in 1:nrow(mean_effects_mixed_lai)){
#       mapped_effects_mixed_lai[site,parameter] = sum(mean_effects_mixed_lai[site,results_mixed$position_mapping[[parameter]]]*
#                                                        results_mixed$weight_mapping[[parameter]])
#     }
#   }
#   mapped_effects_mixed_lai = mapped_effects_mixed_lai[complete.cases(mapped_effects_mixed_lai),]
#   main_effects_mixed_lai_driver = apply(abs(mapped_effects_mixed_lai),2,mean)[ordering_mixed][which(grouping2 == "Drivers")]
#
#   main_effects_mixed_lai_driver = main_effects_mixed_lai_driver/sum(abs(main_effects_mixed_lai_driver))
#
#   result_matrix[counter,] = main_effects_mixed_lai_driver
#   counter = counter +1
# }
#
# colnames(result_matrix) = parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]
#
# mixed_lai_driver_df <- cbind(big_countries_data_set,result_matrix)
# mixed_lai_driver_df = mixed_lai_driver_df[-8,]
#
#
#
# radar.testplot_list <-
#   lapply(1:nrow(pic_abi_cflux_df), function(i) {
#
#     gt_plot <- ggplotGrob(
#       ggradar(rbind(cbind("Name" = "Pic. abi.",pic_abi_cflux_driver_df[i,4:9]),
#                     cbind("Name" = "Fag. syl.",fag_syl_cflux_driver_df[i,4:9]),
#                     cbind("Name" = "Pin. syl.",pin_syl_cflux_driver_df[i,4:9]),
#                     cbind("Name" = "Mixed",mixed_cflux_driver_df[i,4:9])),
#               plot.title = rownames(pic_abi_cflux_driver_df)[i], gridline.max.linetype = 0.4)
#     )
#
#
#
#   })
#
# pdf("./Figures/radar_driver_no_map.pdf", width = 30, height = 25)
# grid.arrange(arrangeGrob(radar.testplot_list[[1]],radar.testplot_list[[2]],radar.testplot_list[[3]], ncol =3, nrow =1),
#              arrangeGrob(radar.testplot_list[[4]],radar.testplot_list[[5]],radar.testplot_list[[6]], ncol =3, nrow =1),
#              arrangeGrob(radar.testplot_list[[7]],radar.testplot_list[[8]],radar.testplot_list[[9]], ncol =3, nrow =1),
#              arrangeGrob(ggplot()+theme_void(),radar.testplot_list[[10]],radar.testplot_list[[11]],ggplot()+theme_void(), ncol =4, nrow =1,
#                          widths = c(0.5,1,1,0.5)),nrow = 4)
# dev.off()
#
#
