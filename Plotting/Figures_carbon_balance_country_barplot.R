sites <- readRDS("sites_data.rds")

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





## Pic abi and cflux

result_matrix = as.data.frame(matrix(nrow = 12, ncol = 7))
rownames(result_matrix) = countries

counter = 1

mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)
sites_pic_abi = sites[growth_sites_Pic_abi,]

for(country in big_enough_countries){
  sites_per_country = which(sites_pic_abi$Country == country)
  if(length(sites_per_country) >= 5){
    mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux"]][["complete"]][growth_sites_Pic_abi,][sites_per_country,]

    main_effects_Pic_abi_cflux = apply(abs(mean_effects_Pic_abi_cflux[,1:38]),2,mean)

    main_effects_Pic_abi_cflux_weight = main_effects_Pic_abi_cflux/sum(abs(main_effects_Pic_abi_cflux))

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
myColors = magma(10)

colors_env = c("darkblue","blue3","chartreuse2","darkgreen","green3",
               "deepskyblue","green1","gold",'yellow',
               "yellow2","darkorange","cyan3","darkolivegreen3")


Pic_abi_carbon_balance <- ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = viridis(10)[10:1]), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env) +
  new_scale("fill") + ggtitle("a) Pic. abi.") + theme(legend.position = "none") +
  coord_sf()

bar.testplot_list <-
  lapply(1:nrow(pic_abi_cflux_df), function(i) {
    new_data = data.frame("processes" = colnames(pic_abi_cflux_df[i,4:10]),
                          "values" = as.numeric(pic_abi_cflux_df[i,4:10]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_me) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(pic_abi_cflux_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = pic_abi_cflux_df[i,"lon"] - 4e5,
                    xmax = pic_abi_cflux_df[i,"lon"] + 4e5,
                    ymin = pic_abi_cflux_df[i,"lat"] - 3e5,
                    ymax = pic_abi_cflux_df[i,"lat"] + 3e5))

Pic_abi_carbon_balance = Reduce("+", bar_annotation_list,Pic_abi_carbon_balance)


Pin_syl_carbon_balance <-ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env) +
  new_scale("fill") + theme(legend.position = "none") +
  coord_sf() + ggtitle("b) Pin. syl.")

bar.testplot_list <-
  lapply(1:nrow(pin_syl_cflux_df), function(i) {
    new_data = data.frame("processes" = colnames(pin_syl_cflux_df[i,4:10]),
                          "values" = as.numeric(pin_syl_cflux_df[i,4:10]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_me) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(pin_syl_cflux_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = pin_syl_cflux_df[i,"lon"] - 4e5,
                    xmax = pin_syl_cflux_df[i,"lon"] + 4e5,
                    ymin = pin_syl_cflux_df[i,"lat"] - 3e5,
                    ymax = pin_syl_cflux_df[i,"lat"] + 3e5))

Pin_syl_carbon_balance = Reduce("+", bar_annotation_list,Pin_syl_carbon_balance)


Fag_syl_carbon_balance <-ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env, name = "Environmental zones") +
  new_scale("fill") +
  coord_sf() +  ggtitle("c) Fag. syl.")

bar.testplot_list <-
  lapply(1:nrow(fag_syl_cflux_df), function(i) {
    new_data = data.frame("processes" = colnames(fag_syl_cflux_df[i,4:10]),
                          "values" = as.numeric(fag_syl_cflux_df[i,4:10]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_me) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(fag_syl_cflux_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = fag_syl_cflux_df[i,"lon"] - 4e5,
                    xmax = fag_syl_cflux_df[i,"lon"] + 4e5,
                    ymin = fag_syl_cflux_df[i,"lat"] - 3e5,
                    ymax = fag_syl_cflux_df[i,"lat"] + 3e5))

Fag_syl_carbon_balance = Reduce("+", bar_annotation_list,Fag_syl_carbon_balance)

mixed_cflux <-ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env, name = "Environmental zones") +
  new_scale("fill") + theme(legend.position = "none") +
  coord_sf() +  ggtitle("d) Mixed stand")


bar.testplot_list <-
  lapply(1:nrow(mixed_cflux_df), function(i) {
    new_data = data.frame("processes" = colnames(mixed_cflux_df[i,4:10]),
                          "values" = as.numeric(mixed_cflux_df[i,4:10]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_me) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(mixed_cflux_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = mixed_cflux_df[i,"lon"] - 4e5,
                    xmax = mixed_cflux_df[i,"lon"] + 4e5,
                    ymin = mixed_cflux_df[i,"lat"] - 3e5,
                    ymax = mixed_cflux_df[i,"lat"] + 3e5))

mixed_cflux = Reduce("+", bar_annotation_list,mixed_cflux)

mixed_lai <-ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env, name = "Environmental zones") +
  new_scale("fill") + theme(legend.position = "none") +
  coord_sf() +  ggtitle("e) Mixed stand")


bar.testplot_list <-
  lapply(1:nrow(mixed_lai_df), function(i) {
    new_data = data.frame("processes" = colnames(mixed_lai_df[i,4:10]),
                          "values" = as.numeric(mixed_lai_df[i,4:10]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_me) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(mixed_lai_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = mixed_lai_df[i,"lon"] - 4e5,
                    xmax = mixed_lai_df[i,"lon"] + 4e5,
                    ymin = mixed_lai_df[i,"lat"] - 3e5,
                    ymax = mixed_lai_df[i,"lat"] + 3e5))

mixed_lai = Reduce("+", bar_annotation_list,mixed_lai)

library(gridExtra)
library(lemon)

legend = g_legend(Fag_syl_carbon_balance)
legdn1 = legend$grobs[["99_58e3232fe5a07e0795cee4636c073554"]]


Fag_syl_carbon_balance = Fag_syl_carbon_balance + theme(legend.position = "none")

library(ggpubr)

title1=text_grob("Carbon Balance", size = 15, face = "bold")

title2=text_grob("Community Composition", size = 15, face = "bold")

pdf("./Figures/maps_lin_bar.pdf", width =20, height =20)

grid.arrange(arrangeGrob(Pic_abi_carbon_balance,Pin_syl_carbon_balance, ncol = 2, nrow = 1, top = title1),
             arrangeGrob(Fag_syl_carbon_balance, mixed_cflux, ncol = 2, nrow = 1),
             arrangeGrob(arrangeGrob(ggplot()+ theme_void(), mixed_lai, widths = c(1,2), ncol =2, nrow =1),arrangeGrob(ggplot()+theme_void(),legdn1, ncol =2),
                         widths = c(3,1),ncol =2, nrow=1, top = title2),
             ncol =1, nrow =3)

dev.off()




##### Carbon balance driver analysis #####

## Pic abi and cflux

result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
rownames(result_matrix) = countries

counter = 1

mean_effects_Pic_abi_cpool = effects_Pic_abi[["cpool"]][["complete"]]
growth_sites_Pic_abi = which(mean_effects_Pic_abi_cpool[,39]/200 > 2)
sites_pic_abi = sites[growth_sites_Pic_abi,]

for(country in big_enough_countries){
  sites_per_country = which(sites_pic_abi$Country == country)
  mean_effects_Pic_abi_cflux = effects_Pic_abi[["cflux"]][["complete"]][growth_sites_Pic_abi,][sites_per_country,]

  main_effects_Pic_abi_cflux_driver = apply(abs(mean_effects_Pic_abi_cflux[,1:38]),2,mean)[which(grouping2 == "Drivers")]

  main_effects_Pic_abi_cflux_weight_driver = main_effects_Pic_abi_cflux_driver/sum(abs(main_effects_Pic_abi_cflux_driver))

  result_matrix[counter,] = main_effects_Pic_abi_cflux_weight_driver

  counter = counter +1
}

colnames(result_matrix) = parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]


pic_abi_cflux_driver_df <- cbind(big_countries_data_set,result_matrix)
pic_abi_cflux_driver_df = pic_abi_cflux_driver_df[-which(pic_abi_cflux_driver_df$country == "Russia"),]



result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
rownames(result_matrix) = countries

counter = 1

mean_effects_fag_syl_cpool = effects_Fag_syl[["cpool"]][["complete"]]
growth_sites_fag_syl = which(mean_effects_fag_syl_cpool[,39]/200 > 2)
sites_fag_syl = sites[growth_sites_fag_syl,]

for(country in big_enough_countries){
  sites_per_country = which(sites_fag_syl$Country == country)
  if(length(sites_per_country) >= 5){
    mean_effects_fag_syl_cflux = effects_Fag_syl[["cflux"]][["complete"]][growth_sites_fag_syl,][sites_per_country,]

    main_effects_fag_syl_cflux_driver = apply(abs(mean_effects_fag_syl_cflux[,1:38]),2,mean)[ordering_fag_syl][which(grouping2 == "Drivers")]

    main_effects_fag_syl_cflux_weight_driver = main_effects_fag_syl_cflux_driver/sum(abs(main_effects_fag_syl_cflux_driver))

    result_matrix[counter,] = main_effects_fag_syl_cflux_weight_driver
  }
  counter = counter +1
}

colnames(result_matrix) = parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]


fag_syl_cflux_driver_df <- cbind(big_countries_data_set,result_matrix)
fag_syl_cflux_driver_df = fag_syl_cflux_driver_df[-which(fag_syl_cflux_driver_df$country == "Russia"),]



result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
rownames(result_matrix) = countries

counter = 1

mean_effects_pin_syl_cpool = effects_Pin_syl[["cpool"]][["complete"]]
growth_sites_pin_syl = which(mean_effects_pin_syl_cpool[,39]/200 > 2)
sites_pin_syl = sites[growth_sites_pin_syl,]

for(country in big_enough_countries){
  sites_per_country = which(sites_pin_syl$Country == country)
  if(length(sites_per_country) >= 5){
    mean_effects_pin_syl_cflux = effects_Pin_syl[["cflux"]][["complete"]][growth_sites_pin_syl,][sites_per_country,]

    main_effects_pin_syl_cflux_driver = apply(abs(mean_effects_pin_syl_cflux[,1:38]),2,mean)[ordering_pin_syl][which(grouping2 == "Drivers")]

    main_effects_pin_syl_cflux_weight_driver = main_effects_pin_syl_cflux_driver/sum(abs(main_effects_pin_syl_cflux_driver))

    result_matrix[counter,] = main_effects_pin_syl_cflux_weight_driver
  }
  counter = counter +1
}

colnames(result_matrix) = parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]


pin_syl_cflux_driver_df <- cbind(big_countries_data_set,result_matrix)
pin_syl_cflux_driver_df = pin_syl_cflux_driver_df[-which(pin_syl_cflux_driver_df$country == "Russia"),]



result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
rownames(result_matrix) = countries

counter = 1

mean_effects_mixed_cpool = effects_mixed[["cpool"]][["complete"]]
growth_sites_mixed = which(mean_effects_mixed_cpool[,74]/200 > 2)
sites_mixed = sites[growth_sites_mixed,]

for(country in big_enough_countries){
  print(sites_per_country)
  sites_per_country = which(sites_mixed$Country == country)
  if(length(sites_per_country) >= 5){
    mean_effects_mixed_cflux = effects_mixed[["cflux"]][["complete"]][growth_sites_mixed,][sites_per_country,]

    remaped_effects_cpool_abs = remap_paramters(mean_effects_mixed_cflux, mixed_results$weight_mapping, mixed_results$position_mapping)

    main_effects_mixed_cflux_driver = apply(abs(remaped_effects_cpool_abs),2,mean)[ordering_mixed][which(grouping2 == "Drivers")]

    main_effects_mixed_cflux_driver_weight = main_effects_mixed_cflux_driver/sum(abs(main_effects_mixed_cflux_driver))

    result_matrix[counter,] = main_effects_mixed_cflux_driver_weight
  }
  counter = counter +1
}

colnames(result_matrix) =  parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]


mixed_cflux_driver_df <- cbind(big_countries_data_set,result_matrix)
mixed_cflux_driver_df = mixed_cflux_driver_df[-which(mixed_cflux_driver_df$country == "Russia"),]

result_matrix = as.data.frame(matrix(nrow = 12, ncol = 6))
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
  main_effects_mixed_lai_driver = apply(abs(mapped_effects_mixed_lai),2,mean)[ordering_mixed][which(grouping2 == "Drivers")]

  main_effects_mixed_lai_driver = main_effects_mixed_lai_driver/sum(abs(main_effects_mixed_lai_driver))

  result_matrix[counter,] = main_effects_mixed_lai_driver
  counter = counter +1
}

colnames(result_matrix) = parameternames_ordered_pic_abi2[which(grouping2 == "Drivers")]

mixed_lai_driver_df <- cbind(big_countries_data_set,result_matrix)
mixed_lai_driver_df = mixed_lai_driver_df[-8,]



colors_driver = c("grey","yellow","green","brown","darkblue","red")
aoi_boundary_HARV <- st_read("./EnSv8/ens_v8.shp")

library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(scatterpie)
myColors = magma(10)

colors_env = c("darkblue","blue3","chartreuse2","darkgreen","green3",
               "deepskyblue","green1","gold",'yellow',
               "yellow2","darkorange","cyan3","darkolivegreen3")


Pic_abi_carbon_balance <- ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = viridis(10)[10:1]), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env) +
  new_scale("fill") + ggtitle("a) Pic. abi.") +
  coord_sf() + theme(legend.position = "none")

pdf("./Figures/blank_map.pdf", height = 5, width = 15)
Pic_abi_carbon_balance
dev.off()

bar.testplot_list <-
  lapply(1:nrow(pic_abi_cflux_driver_df), function(i) {
    new_data = data.frame("processes" = colnames(pic_abi_cflux_driver_df[i,4:9]),
                          "values" = as.numeric(pic_abi_cflux_driver_df[i,4:9]))
    print(str(new_data))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill = colors_driver) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })


bar_annotation_list <- lapply(1:nrow(pic_abi_cflux_driver_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = pic_abi_cflux_driver_df[i,"lon"] - 4e5,
                    xmax = pic_abi_cflux_driver_df[i,"lon"] + 4e5,
                    ymin = pic_abi_cflux_driver_df[i,"lat"] - 3e5,
                    ymax = pic_abi_cflux_driver_df[i,"lat"] + 3e5))

Pic_abi_carbon_balance = Reduce("+", bar_annotation_list,Pic_abi_carbon_balance)

Pin_syl_carbon_balance <-ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env) +
  new_scale("fill") + ggtitle("b) Pin. syl.") +
  coord_sf() +  theme(legend.position = "none")

bar.testplot_list <-
  lapply(1:nrow(pin_syl_cflux_driver_df), function(i) {
    new_data = data.frame("processes" = colnames(pin_syl_cflux_driver_df[i,4:9]),
                          "values" = as.numeric(pin_syl_cflux_driver_df[i,4:9]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_driver) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(pin_syl_cflux_driver_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = pin_syl_cflux_driver_df[i,"lon"] - 4e5,
                    xmax = pin_syl_cflux_driver_df[i,"lon"] + 4e5,
                    ymin = pin_syl_cflux_driver_df[i,"lat"] - 3e5,
                    ymax = pin_syl_cflux_driver_df[i,"lat"] + 3e5))

Pin_syl_carbon_balance = Reduce("+", bar_annotation_list,Pin_syl_carbon_balance)

Fag_syl_carbon_balance <-ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env, name = "Environmental zones") +
  new_scale("fill") +
  coord_sf() +  ggtitle("c) Fag. syl.")


bar.testplot_list <-
  lapply(1:nrow(fag_syl_cflux_driver_df), function(i) {
    new_data = data.frame("processes" = colnames(fag_syl_cflux_driver_df[i,4:9]),
                          "values" = as.numeric(fag_syl_cflux_driver_df[i,4:9]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_driver) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(fag_syl_cflux_driver_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = fag_syl_cflux_driver_df[i,"lon"] - 4e5,
                    xmax = fag_syl_cflux_driver_df[i,"lon"] + 4e5,
                    ymin = fag_syl_cflux_driver_df[i,"lat"] - 3e5,
                    ymax = fag_syl_cflux_driver_df[i,"lat"] + 3e5))

Fag_syl_carbon_balance = Reduce("+", bar_annotation_list,Fag_syl_carbon_balance)

mixed_cflux_driver <-ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env, name = "Environmental zones") +
  new_scale("fill") + theme(legend.position = "none") +
  coord_sf() +  ggtitle("d) Mixed stand")

bar.testplot_list <-
  lapply(1:nrow(mixed_cflux_driver_df), function(i) {
    new_data = data.frame("processes" = colnames(mixed_cflux_driver_df[i,4:9]),
                          "values" = as.numeric(mixed_cflux_driver_df[i,4:9]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_driver) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(mixed_cflux_driver_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = mixed_cflux_driver_df[i,"lon"] - 4e5,
                    xmax = mixed_cflux_driver_df[i,"lon"] + 4e5,
                    ymin = mixed_cflux_driver_df[i,"lat"] - 3e5,
                    ymax = mixed_cflux_driver_df[i,"lat"] + 3e5))

mixed_cflux_driver = Reduce("+", bar_annotation_list,mixed_cflux_driver)


mixed_lai_driver <-ggplot() +
  geom_sf(data = aoi_boundary_HARV, size = 0.01, color = "black",
          aes(fill = EnZ_name, colour = magma(10)), alpha = .4) +
  theme_void()+ scale_fill_manual(values = colors_env, name = "Environmental zones") +
  new_scale("fill") + theme(legend.position = "none") +
  coord_sf() +  ggtitle("d) Mixed stand")

bar.testplot_list <-
  lapply(1:nrow(mixed_lai_driver_df), function(i) {
    new_data = data.frame("processes" = colnames(mixed_lai_driver_df[i,4:9]),
                          "values" = as.numeric(mixed_lai_driver_df[i,4:9]))
    gt_plot <- ggplotGrob(
      ggplot(data = new_data, aes( x = processes, y = values))+ theme_void()+
        geom_bar(stat="identity", color="black", fill=colors_driver) + ylim (0,0.8)
    )
    panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel",]
    gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
  })



bar_annotation_list <- lapply(1:nrow(mixed_lai_driver_df), function(i)
  annotation_custom(bar.testplot_list[[i]],
                    xmin = mixed_lai_driver_df[i,"lon"] - 4e5,
                    xmax = mixed_lai_driver_df[i,"lon"] + 4e5,
                    ymin = mixed_lai_driver_df[i,"lat"] - 3e5,
                    ymax = mixed_lai_driver_df[i,"lat"] + 3e5))

mixed_lai_driver = Reduce("+", bar_annotation_list,mixed_lai_driver)

library(gridExtra)
library(lemon)

legend = g_legend(Fag_syl_carbon_balance)
legdn1 = legend$grobs[["99_58e3232fe5a07e0795cee4636c073554"]]


Fag_syl_carbon_balance = Fag_syl_carbon_balance + theme(legend.position = "none")

library(ggpubr)

title1=text_grob("Carbon Balance", size = 15, face = "bold")

title2=text_grob("Community Composition", size = 15, face = "bold")

pdf("./Figures/maps_lin_driver_bar.pdf", width =20, height =20)

grid.arrange(arrangeGrob(Pic_abi_carbon_balance,Pin_syl_carbon_balance, ncol = 2, nrow = 1, top = title1),
             arrangeGrob(Fag_syl_carbon_balance, mixed_cflux_driver, ncol = 2, nrow = 1),
             arrangeGrob(arrangeGrob(ggplot()+ theme_void(), mixed_lai_driver, widths = c(1,2), ncol =2, nrow =1),
                         arrangeGrob(ggplot() + theme_void(),legdn1, ncol =2),
                         widths = c(3,1),ncol =2, nrow=1, top = title2),
             ncol =1, nrow =3)

dev.off()
